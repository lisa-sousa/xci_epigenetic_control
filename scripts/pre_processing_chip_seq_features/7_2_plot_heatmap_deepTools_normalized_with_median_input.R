#normalize ChIP-Seq counts with median input and save output as bigwig file (25bp bins covering +/- 2200bp around gene TSS)

library(GenomicRanges)
library(rtracklayer)

#directories
deepTools = '/scratch/ngsvin/Chip-seq/deepTools/bin/'
bedTools = '/home/lisasous/programs_libs/programs/bedtools2/bin/'
file_genes = '/project/lncrna/Xist/data/silencing_rates/fitted_data/silencing_rates_pro_seq_mm9.bed'
file_metadata = '/project/lncrna/Xist/data/metadata/metadata_filtered_by_deepTools.txt'
ChIP_dir = '/project/ngs_marsico/Xist/bam/experiment/'
ctrl_dir = '/project/ngs_marsico/Xist/bam/control/'
output_dir_data = '/project/lncrna/Xist/data/ChIP_Seq/deepTools_heatmap/'
output_dir_plot = '/project/lncrna/Xist/scripts/ChIP_seq_analysis/3_heatmap/plots/'
gene_region1 = '/project/lncrna/Xist/data/silencing_rates/fitted_data/silencing_rates_pro_seq_mm9_fastly_silenced.bed'
gene_region2 = '/project/lncrna/Xist/data/silencing_rates/fitted_data/silencing_rates_pro_seq_mm9_slowly_silenced.bed'
gene_region = '/project/lncrna/Xist/data/silencing_rates/fitted_data/silencing_rates_pro_seq_mm9.bed'
number_of_gene_regions = 1

#regions around TSS
a_norm = 1100
b_norm = 1100
a_heatmap = 2000
b_heatmap = 2000

#metadata
metadata = read.table(file_metadata,sep='\t',header=T,comment.char='')
colnames(metadata)[1] = 'feature'
metadata$feature = as.character(metadata$feature)
metadata$accession_number = as.character(metadata$accession_number)
metadata$experiment_file = as.character(metadata$experiment_file)
metadata$control_file = as.character(metadata$control_file)

#gene regions
table_genes = read.table(file_genes,sep='\t',header=T,comment.char = '')
genes = GRanges(seqnames = table_genes$chromosome, ranges = IRanges(table_genes$start,table_genes$end), strand = table_genes$strand)
gene_regions = promoters(genes,upstream = a_norm,downstream = b_norm)

si = Seqinfo(genome='mm9')
gr = GRanges(seqnames = 'chrX',ranges=IRanges(1,seqlengths(si)[seqnames(si) == 'chrX']))
binned_genome_regions = unlist(tile(gr,width = 25)[1])
seqinfo(binned_genome_regions) = Seqinfo(seqnames = 'chrX', seqlength = seqlengths(si)[seqnames(si) == 'chrX'], isCircular=isCircular(si)[seqnames(si)=='chrX'],genome='mm9')

#temporary files
cmd = paste('mkdir ', output_dir_data, 'tmp', sep='')
system(cmd)
file_bed = paste(output_dir_data,'tmp/tmp_gene_regions.bed',sep='')
file_bed_results_ChIP = paste(output_dir_data,'tmp/tmp_coverage_ChIP.bed',sep='')
file_bed_results_Ctrl = paste(output_dir_data,'tmp/tmp_coverage_Ctrl.bed',sep='')

#generate one bigWig and one heatmap per feature
for (i in 1:nrow(metadata)) {
  feature = metadata$feature[i]
  print(feature)
  acNr = metadata$accession_number[i]
  title = paste(feature,acNr,'normalized_median_input',sep='_')
  bamfileChIP = paste(ChIP_dir,metadata$experiment_file[i],sep='')
  bamfileCtrl = paste(ctrl_dir,metadata$control_file[i],sep='')
  outfile_bw = paste(output_dir_data,'big_wig/',title,'.bw',sep='')
  
  if(!file.exists(outfile_bw)){
    #####get global scaling factor
    #export gene regions to bed file
    export.bed(gene_regions,file_bed,format='bed')
    
    #get coverage of ChIP
    #uses bedTools v 2.25.0 -> -a and -b are swapped!!
    cmd = paste(bedTools, 'coverageBed -a ', file_bed, ' -b ', bamfileChIP, ' > ', file_bed_results_ChIP, sep='')
    system(cmd)
    
    bed_results_ChIP = import(file_bed_results_ChIP,extraCols=c('a','b','c','d'))
    read_counts_ChIP = as.integer(mcols(bed_results_ChIP)[,3])
    
    #get coverage of Control
    cmd = paste(bedTools, 'coverageBed -a ', file_bed, ' -b ', bamfileCtrl, ' > ', file_bed_results_Ctrl, sep='')
    system(cmd)
    
    bed_results_Ctrl = import(file_bed_results_Ctrl,extraCols=c('a','b','c','d'))
    read_counts_Ctrl = as.integer(mcols(bed_results_Ctrl)[,3])
    
    #calculate scaling factor
    counts_table = data.frame(read_counts_ChIP = read_counts_ChIP, read_counts_Ctrl = read_counts_Ctrl)
    counts_table$fc = ((counts_table$read_counts_ChIP + 1.0)/(counts_table$read_counts_Ctrl + 1))
    global_scaling_factor = median(counts_table$fc)
    
    #delete tmp files
    cmd = paste('rm -f ', output_dir_data, 'tmp/*.bed', sep='')
    system(cmd)
    
    #####normalize bins on whole X chromosome
    #export binned gene regions as bed file
    export.bed(binned_genome_regions,file_bed,format='bed')
    
    #get coverage of ChIP
    cmd = paste(bedTools, 'coverageBed -a ', file_bed, ' -b ', bamfileChIP, ' > ', file_bed_results_ChIP, sep='')
    system(cmd)
    
    bed_results_ChIP = import(file_bed_results_ChIP,extraCols=c('a','b','c','d'))
    read_counts_ChIP = as.integer(mcols(bed_results_ChIP)[,3])
    
    #get coverage of Control
    cmd = paste(bedTools, 'coverageBed -a ', file_bed, ' -b ', bamfileCtrl, ' > ', file_bed_results_Ctrl, sep='')
    system(cmd)
    
    bed_results_Ctrl = import(file_bed_results_Ctrl,extraCols=c('a','b','c','d'))
    read_counts_Ctrl = as.integer(mcols(bed_results_Ctrl)[,3])
    
    #normalize counts
    counts_table = data.frame(read_counts_ChIP = read_counts_ChIP, read_counts_Ctrl = read_counts_Ctrl)
    normalized_counts = ((counts_table$read_counts_ChIP + 1.0)/(counts_table$read_counts_Ctrl + 1.0))*(1.0/global_scaling_factor)
    
    #append normalized counts for bigWig file
    normalized_binned_genome_regions = binned_genome_regions
    mcols(normalized_binned_genome_regions) = data.frame(score = normalized_counts)
    
    #delete tmp files
    cmd = paste('rm -f ', output_dir_data, 'tmp/*.bed', sep='')
    system(cmd)

    export.bw(normalized_binned_genome_regions, outfile_bw, format='bw')
  }
  
  if(number_of_gene_regions == 1){
    mat_file = paste(output_dir_data,'matrix/',title,'.mat.gz',sep='')
    cmd = paste(deepTools,'computeMatrix reference-point -S ',outfile_bw,' -R ',gene_region,' --referencePoint TSS -a ',a_heatmap,' -b ',b_heatmap,
                ' --binSize 50 --numberOfProcessors max/2 --outFileName ', mat_file, sep='')
    print(cmd)
    system(cmd)
    
    plotFile = paste(output_dir_plot,title,'.pdf',sep='')
    cmd = paste(deepTools,'plotHeatmap -m ',mat_file,' -out ',plotFile,' --sortRegions ascend --zMin 0 --colorList white,darksalmon,brown --colorNumber 20',sep='')
    print(cmd)
    system(cmd)
  }else{
    mat_file = paste(output_dir_data,'matrix/',title,'.mat.gz',sep='')
    cmd = paste(deepTools,'computeMatrix reference-point -S ',outfile_bw,' -R ',gene_region1,' ',gene_region2,' --referencePoint TSS -a ',a_heatmap,' -b ',b_heatmap,
                ' --binSize 50 --numberOfProcessors max/2 --outFileName ', mat_file, sep='')
    print(cmd)
    system(cmd)
    
    plotFile = paste(output_dir_plot,title,'.pdf',sep='')
    cmd = paste(deepTools,'plotHeatmap -m ',mat_file,' -out ',plotFile,' --sortRegions ascend --zMin 0 --colorList white,darksalmon,brown --colorNumber 20 --regionsLabel fast slow',sep='')
    print(cmd)
    system(cmd)
  }

}

cmd = paste('pdfunite ',output_dir_plot,'*.pdf ',output_dir_plot,'heatmaps.pdf',sep='')
system(cmd)

cmd = paste('rm -rf ', output_dir_data, 'tmp', sep='')
system(cmd)

