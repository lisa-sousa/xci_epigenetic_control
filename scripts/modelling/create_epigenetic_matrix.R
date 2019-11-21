###################################################################################
#libraries and tools
###################################################################################

library(energy)
library(rtracklayer)
library(GenomicRanges)
library(Biostrings)
library(openxlsx)
library(here)
bedTools = '/home/lisasous/tools/bedtools2/bin/' #set path to bedtools 2 /bin/ folder


###################################################################################
#directories
###################################################################################

feature_matrix_name = 'promoter_matrix_reannotated_normRAdjusted_pro_seq_genes'
output_dir = 'data/modelling/feature_matrix'
tmp_dir = here(output_dir,'tmp')
system(paste('mkdir',tmp_dir))


###################################################################################
#load gene annotation
###################################################################################

##########genes mm9 from gencode
genes_mm9_gencode_file = here('data/annotation_files/gene_annotation','gencode.vM9.annotation.chrX.genes.mm9.bed')
genes_mm9_gencode_table = read.table(genes_mm9_gencode_file,sep = '\t',header=F)
colnames(genes_mm9_gencode_table) = c('chr','start','end','gene_name','score','strand')


##########silencing halftimes
halftimes_file = here('data/silencing_halftimes/fitted_data','halftimes_pro_seq_mm9_reannotated_with_rr.bed')
#halftimes_file = here('data/annotation_files/gene_annotation','gencode.vM9.annotation.chrX.genes.reannotated.with.rr.mm9.bed') #all genes
halftimes_table = read.table(file = halftimes_file, sep = '\t', header = F)
colnames(halftimes_table) = c('chr','start','end','gene_name','halftime','strand')

halftimes_table$TSS[halftimes_table$strand == '+'] = halftimes_table$start[halftimes_table$strand == '+']
halftimes_table$TSS[halftimes_table$strand == '-'] = halftimes_table$end[halftimes_table$strand == '-']


##########gene regions from start to end and 1000 bp around the promoter
gene_regions = GRanges(seqnames = halftimes_table$chr, ranges = IRanges(halftimes_table$start,halftimes_table$end), strand = halftimes_table$strand)
mcols(gene_regions) = data.frame(name = halftimes_table$gene_name, score = halftimes_table$halftime)
#mcols(gene_regions) = data.frame(name = halftimes_table$gene_name, score = 0) #all genes

tmp_gene_regions_file = paste(tmp_dir,'/tmp_gene_regions_file.bed',sep='')
export.bed(gene_regions,tmp_gene_regions_file,format='bed')

gene_regions_1000 = promoters(gene_regions, upstream = 500, downstream = 500)
tmp_gene_regions_1000_file = paste(tmp_dir,'/tmp_gene_regions_1000_file.bed',sep='')
export.bed(gene_regions_1000,tmp_gene_regions_1000_file,format='bed')


###################################################################################
#load chip-seq data
###################################################################################

input_dir = here('data/chip_seq/normalized_counts/reannotated_normRAdjusted_pro_seq_genes/')
#input_dir = here('data/chip_seq/normalized_counts/reannotated_normRAdjusted_all_chrX_genes/') #all genes

#metadata
file_metadata = here('data/chip_seq/metadata','metadata_normalization_regions.txt')
metadata = read.table(file_metadata,sep='\t',header=T,comment.char='')
colnames(metadata)[1] = 'feature'
metadata$feature = as.character(metadata$feature)
metadata$accession_number = as.character(metadata$accession_number)
metadata$experiment_file = as.character(metadata$experiment_file)
metadata$control_file = as.character(metadata$control_file)

chip_seq_table = data.frame(gene_name = halftimes_table$gene_name)

for (i in 1:nrow(metadata)) {
  feature = metadata$feature[i]
  print(feature)
  acNr = metadata$accession_number[i]
  a = metadata$start[i]
  b = metadata$end[i]
  title = paste(feature,acNr,a,b,sep='_')
  bed_file = paste(input_dir,title,'_normalized_normr.bed',sep='')
  
  chip_seq = read.table(bed_file,sep='\t')
  chip_seq = data.frame(gene_name = chip_seq$V4, score = chip_seq$V5)
  colnames(chip_seq) = c('gene_name',title)
  
  chip_seq_table = merge(chip_seq_table,chip_seq,by='gene_name')
}


###################################################################################
#compute genomic features
###################################################################################

##########distance to closest topological domain (TAD) (mapped to mm9)
TADs_file = here('data/annotation_files/TADs','nature11082-s2.xlsx')
TADs = read.xlsx(TADs_file, sheet = 4, colNames = F)
TADs = c(TADs$X2[TADs$X1=="chrX"],TADs$X3[TADs$X1=="chrX"])
all_distances <- abs(outer(halftimes_table$TSS,TADs,FUN="-"))
distance_TAD_table = data.frame(gene_name = halftimes_table$gene_name, distance_to_TAD_border = apply(all_distances,1,min))


##########gene density in 1mb around each gene
tmp_gene_density_file = paste(tmp_dir,'/tmp_gene_density.txt',sep='')
cmd = paste(bedTools,'windowBed -a ',tmp_gene_regions_file,' -b ',genes_mm9_gencode_file,' -w 500000 -c > ',tmp_gene_density_file,sep='') #-w 	Base pairs added upstream and downstream of each entry in A when searching for overlaps in B. Default is 1000 bp.
system(cmd)
gene_density_gencode_table = read.table(tmp_gene_density_file,header=F,sep='\t')[,c(4,7)]
colnames(gene_density_gencode_table) = c('gene_name','gene_density')


##########overlap with Xist initiation site (mapped to mm9)
initiation_sites_file = here('data/annotation_files/xist_initiation_sites','1hr-Xist-Initiation-Sites_mm9.bed')
tmp_overlap_initiation_file = paste(tmp_dir,'/tmp_initiation_sites.txt',sep='')
cmd = paste(bedTools,'intersectBed -a ',tmp_gene_regions_file,' -b ',initiation_sites_file,' -c > ',tmp_overlap_initiation_file,sep='')
system(cmd)
overlap_initiation_table = read.table(tmp_overlap_initiation_file, sep='\t', header=F)
overlap_initiation_table = data.frame(gene_name = overlap_initiation_table$V4, overlap_with_Xist_entry_sites = overlap_initiation_table$V7)
overlap_initiation_table$overlap_with_Xist_entry_sites[overlap_initiation_table$overlap_with_Xist_entry_sites > 0] = 1


##########distance to the Xist locus
TSS_Xist = genes_mm9_gencode_table$end[genes_mm9_gencode_table$gene_name == 'Xist'] #Xist is on the minus strand
distance_xist_table = data.frame(gene_name = halftimes_table$gene_name, distance_to_Xist = abs(halftimes_table$TSS - TSS_Xist))


##########distance (closest) and overlap to Lamina Associated Domain (LADs)
LADs_file = here('data/annotation_files/LADs','1-s2.0-S1097276510003217-mmc2.xlsx')
LADs_table = read.xlsx(LADs_file, sheet = 1)
LADs_table = LADs_table[LADs_table$seqname == 'chrX',]
tmp_LADs_file = paste(tmp_dir,'/tmp_LADs_file.bed',sep='')
export.bed(LADs_table,tmp_LADs_file,format='bed')

tmp_overlap_LADs_file = paste(tmp_dir,'/tmp_overlap_LADs_file.bed',sep='')
cmd = paste(bedTools,'intersectBed -a ',tmp_gene_regions_1000_file,' -b ',tmp_LADs_file,' -c > ',tmp_overlap_LADs_file,sep='')
system(cmd)
overlap_LADs_table = read.table(tmp_overlap_LADs_file, sep='\t', header=F)
overlap_LADs_table = data.frame(gene_name = overlap_LADs_table$V4, overlap_with_LADs = overlap_LADs_table$V7)

LADs = c(LADs_table$start,LADs_table$end)
all_distances <- abs(outer(halftimes_table$TSS,LADs,FUN="-"))
distance_LAD_table = data.frame(gene_name = halftimes_table$gene_name, distance_to_LADs = apply(all_distances,1,min))


##########Hi-C interaction between gene promoter and other regulatory regions (mapped to mm9)
hic_promoter_other_file = here('data/annotation_files/Hi-C_data','ESC_promoter_other_significant_interactions_HIC.txt')
hic_promoter_other_table = read.table(hic_promoter_other_file,sep='\t',header=T)

hic_promoter_promoter_file = here('data/annotation_files/Hi-C_data','ESC_promoter_promoter_significant_interactions_HIC.txt')
hic_promoter_promoter_table = read.table(hic_promoter_promoter_file,sep='\t',header=T)

hic_table = data.frame(gene_name = halftimes_table$gene_name, number_interactions_HiC_all = 0, mean_interaction_strength_HiC_all = 0, number_interactions_HiC_promoter = 0, 
                       mean_interaction_strength_HiC_promoter = 0, mean_interaction_strength_HiC_xist = 0)
for(i in 1:nrow(hic_table)){
  interaction_other = hic_promoter_other_table[grepl(hic_table$gene_name[i],hic_promoter_other_table$Symbol),]
  interaction_promoter = hic_promoter_promoter_table[(grepl(hic_table$gene_name[i],hic_promoter_promoter_table$Symbol.1)) | (grepl(hic_table$gene_name[i],hic_promoter_promoter_table$Symbol)),]
  hic_table$number_interactions_HiC_all[i] = nrow(interaction_other)
  hic_table$mean_interaction_strength_HiC_all[i] = sum(interaction_other$log.observed.expected.)/nrow(interaction_other)
  hic_table$number_interactions_HiC_promoter[i] = nrow(interaction_promoter)
  hic_table$mean_interaction_strength_HiC_promoter[i] = sum(interaction_promoter$log.observed.expected.)/nrow(interaction_promoter)
  hic_table$mean_interaction_strength_HiC_xist[i] = sum(interaction_other$log.observed.expected.[(interaction_other$chr == 'chrX') & (interaction_other$start == '100647076')],
                                          interaction_promoter$log.observed.expected.[(grepl(as.factor('Xist'),interaction_promoter$Symbol.1)) | (grepl(as.factor('Xist'),interaction_promoter$Symbol))])
}
hic_table[is.na(hic_table)] = 0


##########HICap interaction between gene promoter and other regulatory regions (mapped to mm9)
hicap_promoter_enhancer_file = here('data/annotation_files/HiCap_data','Promoter_Enhancer_Interactions.csv')
hicap_promoter_enhancer_table = read.table(hicap_promoter_enhancer_file,sep='\t',header=T)

hicap_promoter_promoter_file = here('data/annotation_files/HiCap_data','Promoter_Promoter_Interactions.csv')
hicap_promoter_promoter_table = read.table(hicap_promoter_promoter_file,sep='\t',header=T)

hicap_table = data.frame(gene_name = halftimes_table$gene_name, number_interactions_HiCap_promoter = 0, number_interactions_HiCap_enhancer = 0, number_interactions_HiCap_all = 0)
for(i in 1:nrow(hicap_table)){
  interactions_promoter = hicap_promoter_promoter_table[as.character(hicap_promoter_promoter_table$gene_name) == as.character(hicap_table$gene_name[i]),]
  interactions_enhancer = hicap_promoter_enhancer_table[as.character(hicap_promoter_enhancer_table$gene_name) == as.character(hicap_table$gene_name[i]),]
  hicap_table$number_interactions_HiCap_promoter[i] = sum(interactions_promoter$read_pairs_in_HiCap_replicate1,interactions_promoter$read_pairs_in_HiCap_replicate2)/(2*nrow(interactions_promoter)) #mean connectivity for two replicates
  hicap_table$number_interactions_HiCap_enhancer[i] = sum(interactions_enhancer$read_pairs_in_HiCap_replicate1,interactions_enhancer$read_pairs_in_HiCap_replicate2)/(2*nrow(interactions_enhancer)) #mean connectivity for two replicates
  hicap_table[is.na(hicap_table)] = 0
  hicap_table$number_interactions_HiCap_all[i] = hicap_table$number_interactions_HiCap_promoter[i] + hicap_table$number_interactions_HiCap_enhancer[i]
}


##########DNA methylation 
dna_methylation_file = here('data/annotation_files/Bisulfite_Seq_data','GSE30202_BisSeq_ES_CpGmeth.tsv.gz')
dna_methylation_table = read.table(dna_methylation_file,sep='\t',header=T)
dna_methylation_table = dna_methylation_table[dna_methylation_table$chr == 'chrX',]
dna_methylation_table = data.frame(chr = dna_methylation_table$chr, start = dna_methylation_table$position, end = dna_methylation_table$position, name = 'meth_base',
                                   score = dna_methylation_table$nMeth/dna_methylation_table$nTot, strand = '.')

tmp_dna_methylation_file = paste(tmp_dir,'/tmp_dna_methylation_file.bed',sep='')
export.bed(dna_methylation_table,tmp_dna_methylation_file,format='bed')

tmp_overlap_dna_methylation_file = paste(tmp_dir,'/tmp_overlap_dna_methylation_file.bed',sep='')
cmd = paste(bedTools,'intersectBed -a ',tmp_gene_regions_1000_file,' -b ',tmp_dna_methylation_file,' -wo > ',tmp_overlap_dna_methylation_file,sep='')
system(cmd)

overlap_dna_methylation_table = read.table(tmp_overlap_dna_methylation_file,sep='\t')
overlap_dna_methylation_table$V4 = as.character(overlap_dna_methylation_table$V4)

dna_methylation_table = data.frame(gene_name = as.character(halftimes_table$gene_name), DNA_methylation_BS_Seq = 0)
for(i in 1:nrow(dna_methylation_table)){
  dna_methylation = overlap_dna_methylation_table[overlap_dna_methylation_table$V4 == dna_methylation_table$gene_name[i],]
  dna_methylation_table$DNA_methylation_BS_Seq[i] = mean(dna_methylation$V11)
}
dna_methylation_table$DNA_methylation_BS_Seq[is.na(dna_methylation_table$DNA_methylation_BS_Seq)] = 0


##########CpG islands
CpG_island_file = here('data/annotation_files/CpG_islands','CGI_mm9_ucsc')
CpG_island_table = read.table(CpG_island_file,sep='\t')
CpG_island_table = data.frame(chr = CpG_island_table$V2, start = CpG_island_table$V3, end = CpG_island_table$V4, name = 'CpG_island', score = 0, strand = '+')

tmp_CpG_island_file = paste(tmp_dir,'/tmp_CpG_island_file.bed',sep='')
export.bed(CpG_island_table,tmp_CpG_island_file,format='bed')

tmp_overlap_CpG_island_file = paste(tmp_dir,'/tmp_overlap_CpG_island_file.bed',sep='')
cmd = paste(bedTools,'intersectBed -a ',tmp_gene_regions_1000_file,' -b ',tmp_CpG_island_file,' -c > ',tmp_overlap_CpG_island_file,sep='')
system(cmd)
overlap_CpG_island_table = read.table(tmp_overlap_CpG_island_file, sep='\t', header=F)
CpG_island_table = data.frame(gene_name = overlap_CpG_island_table$V4, overlap_with_CpG_islands = overlap_CpG_island_table$V7)


##########CpG content
mm9_fasta = here('data/annotation_files/mouse_genome','mm9.fa') ## data not provided, please download mm9 genome in fasta format

tmp_gene_regions_1000_fasta = paste(tmp_dir,'/tmp_gene_regions_1000_file.fa',sep='')
cmd = paste(bedTools,'fastaFromBed -fi ',mm9_fasta,' -bed ',tmp_gene_regions_1000_file,' -fo ',tmp_gene_regions_1000_fasta,sep='')
system(cmd)

sequences = readDNAStringSet(tmp_gene_regions_1000_fasta)
dinucleotide_frequency = dinucleotideFrequency(sequences,step = 1)
frequencies = data.frame(CpG = dinucleotide_frequency[,colnames(dinucleotide_frequency) == 'CG'], G = letterFrequency(sequences,letters = 'G'), C = letterFrequency(sequences,letters = 'C'))
L = width(sequences)
CpG_content_table = data.frame(gene_name = halftimes_table$gene_name, CpG_content = (frequencies$CpG/L) / (((frequencies$G+frequencies$C)/(2*L))^2)) #normalization of CpG content taken from: http://nar.oxfordjournals.org/content/37/19/6305.full.pdf+html


##########distance (closest) to LINE and LINE density in 700kb around each gene
LINEs_file = here('data/annotation_files/LINEs','LINEs_mm9.bed')
LINEs_table = read.table(LINEs_file,sep='\t')
colnames(LINEs_table) = c("chr","start","end","name","score","strand","dispStart","dispEnd","col")
LINEs_table = LINEs_table[LINEs_table$chr == 'chrX',]

tmp_LINEs_file = paste(tmp_dir,'/tmp_LINEs_file.bed',sep='')
export.bed(LINEs_table,tmp_LINEs_file,format='bed')

tmp_LINE_density_file = paste(tmp_dir,'/tmp_LINE_density.txt',sep='')
cmd = paste(bedTools,'windowBed -a ',tmp_gene_regions_file,' -b ',tmp_LINEs_file,' -w 350000 -c > ',tmp_LINE_density_file,sep='') #-w 	Base pairs added upstream and downstream of each entry in A when searching for overlaps in B. Default is 1000 bp.
system(cmd)
LINE_density_table = read.table(tmp_LINE_density_file,header=F,sep='\t')[,c(4,7)]
colnames(LINE_density_table) = c('gene_name','LINE_density')

LINEs = c(LINEs_table$start,LINEs_table$end)
all_distances = abs(outer(halftimes_table$TSS,LINEs,FUN="-"))
distance_LINE_table = data.frame(gene_name = halftimes_table$gene_name, distance_to_LINE = apply(all_distances,1,min))


##########remove tmp files
system(paste('rm -rf',tmp_dir))


###################################################################################
#merge all tables
###################################################################################

merged_tables = merge(halftimes_table[,c(4,5)], chip_seq_table, by = 'gene_name', all.x = T, all.y = F)
merged_tables = merge(merged_tables, distance_TAD_table, by = 'gene_name', all.x = T, all.y = F)
merged_tables = merge(merged_tables, gene_density_gencode_table, by = 'gene_name', all.x = T, all.y = F)
merged_tables = merge(merged_tables, overlap_initiation_table, by = 'gene_name', all.x = T, all.y = F)
merged_tables = merge(merged_tables, distance_xist_table, by = 'gene_name', all.x = T, all.y = F)
merged_tables = merge(merged_tables, overlap_LADs_table, by = 'gene_name', all.x = T, all.y = F)
merged_tables = merge(merged_tables, distance_LAD_table, by = 'gene_name', all.x = T, all.y = F)
merged_tables = merge(merged_tables, hic_table, by = 'gene_name', all.x = T, all.y = F)
merged_tables = merge(merged_tables, hicap_table, by = 'gene_name', all.x = T, all.y = F)
merged_tables = merge(merged_tables, dna_methylation_table, by = 'gene_name', all.x = T, all.y = F)
merged_tables = merge(merged_tables, CpG_island_table, by = 'gene_name', all.x = T, all.y = F)
merged_tables = merge(merged_tables, CpG_content_table, by = 'gene_name', all.x = T, all.y = F)
merged_tables = merge(merged_tables, LINE_density_table, by = 'gene_name', all.x = T, all.y = F)
merged_tables = merge(merged_tables, distance_LINE_table, by = 'gene_name', all.x = T, all.y = F)


###################################################################################
#create matrix and save it
###################################################################################

matrix_table = merged_tables
row.names(matrix_table) = matrix_table$gene_name
matrix_table$gene_name = NULL
print(nrow(matrix_table))

#remove halftime from tablebecause it is used as target
halftime = matrix_table$halftime
matrix_table$halftime = NULL

#make certain features as factors for nomalization
matrix_table$overlap_with_Xist_entry_sites = as.factor(matrix_table$overlap_with_Xist_entry_sites)
matrix_table$overlap_with_CpG_islands = as.factor(matrix_table$overlap_with_CpG_islands)
matrix_table$overlap_with_LADs = as.factor(matrix_table$overlap_with_LADs)

#save data_set and target
data_set = matrix_table
halftime = halftime

save(data_set,halftime,file = here(output_dir,paste(feature_matrix_name,'.RData',sep='')))


