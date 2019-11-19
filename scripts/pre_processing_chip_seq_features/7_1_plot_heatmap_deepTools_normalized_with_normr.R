#normalize ChIP-Seq counts with normr and save output as bigwig file (x bp bins covering the whole genome)

library(normr)
library(GenomicRanges)
library(rtracklayer)
library(bamsignals)
source("/project/lncrna/Xist/scripts/chip_seq_analysis/getEnrichmentFine.R")

deepTools = '/home/lisasous/tools/deepTools2.0/bin/' 
file_metadata = '/project/lncrna/Xist/data/chip_seq/metadata/metadata_filtered_by_deepTools.txt'
file_mm9_chrom_sizes = '/project/lncrna/Xist/data/annotation_files/mouse_genome/mm9.chrom.sizes'
ChIP_dir = '/project/ngs_marsico/Xist/bam/experiment/'
ctrl_dir = '/project/ngs_marsico/Xist/bam/control/'
output_dir_data = '/project/lncrna/Xist/data/chip_seq/bigwigs/normr/'
output_dir_plot = '/project/lncrna/Xist/plots/chip_seq_analysis/heatmap/data_normalized_normR/'
#gene_regions = c('/project/lncrna/Xist/data/modelling/clustering/silenced_vs_not_silenced_cluster1.txt','/project/lncrna/Xist/data_lisa/modelling/clustering/silenced_vs_not_silenced_cluster2.txt','/project/lncrna/Xist/data_lisa/modelling/clustering/silenced_vs_not_silenced_cluster3.txt','/project/lncrna/Xist/data_lisa/modelling/clustering/silenced_vs_not_silenced_cluster4.txt')
#gene_regions = c('/project/lncrna/Xist/cl1.txt','/project/lncrna/Xist/cl2.txt','/project/lncrna/Xist/cl3.txt')
gene_regions = '/project/lncrna/Xist/data/silencing_halftimes/fitted_data/halftimes_pro_seq_mm9_new_gene_annotation.bed'

a = 2000
b = 2000
binsize.bigwig = 25
binsize.heatmap = 100
scale.pseudo = F

metadata = read.table(file_metadata,sep='\t',header=T,comment.char='')
colnames(metadata)[1] = 'feature'
metadata$feature = as.character(metadata$feature)
metadata$accession_number = as.character(metadata$accession_number)
metadata$experiment_file = as.character(metadata$experiment_file)
metadata$control_file = as.character(metadata$control_file)

#chromosome
mm9_chrom_sizes = read.table(file = file_mm9_chrom_sizes, header=F, sep='\t')
colnames(mm9_chrom_sizes) = c('chr','length')

#create temporary directory
cmd = paste("mkdir ",output_dir_data,"tmp",sep="")
print(cmd)
system(cmd)


for (i in 1:nrow(metadata)) {
  feature = metadata$feature[i]
  print(feature)
  acNr = metadata$accession_number[i]
  title = paste(feature,acNr,'normalized_normr',sep='_')
  bamfileChIP = paste(ChIP_dir,metadata$experiment_file[i],sep='')
  bamfileCtrl = paste(ctrl_dir,metadata$control_file[i],sep='')
  outfile_bw_enrichment = paste(output_dir_data,title,'_enrichment.bw',sep='')
 
  if(!file.exists(outfile_bw_enrichment)){
    #fit the data and compute normalization parameters, default binsize=250 --> is too small --> new binsize=500
    fit <- enrichR(bamfileChIP, bamfileCtrl, genome=mm9_chrom_sizes,countConfigSingleEnd(binsize = 500))
    # ctsChIP <- getCounts(fit)$treatment
    # ctsCtrl <- getCounts(fit)$control
    # post <- getPosteriors(fit)[,1]
    # 
    # binsize.fit = width(getRanges(fit))[2]
    # 
    # pseudoChIP <- sum(post * ctsChIP)/sum(post) / binsize.fit / binsize.vis
    # pseudoCtrl <- sum(post * ctsCtrl)/sum(post) / binsize.fit / binsize.vis
    # regul <- log(pseudoCtrl/pseudoChIP)
    # stdrz <- log(fit@theta[2]/(1-fit@theta[2]) * (1-fit@theta[1])/fit@theta[1])
    # 
    # #divide the genome in 25bp bins and get the read counts for each bin 
    # #genome = data.frame(chr = mm9_chrom_sizes$chr, length=mm9_chrom_sizes$length)
    # #binned_counts = normr:::handleCharCharDf(bamfileChIP, bamfileCtrl, genome, countConfigSingleEnd(binsize = binsize.vis), procs = 4,verbose = F)
    # 
    # coverage = bamCoverage(bamfileChIP, chrX, mapqual = 10)
    # 
    # #compute foldchange of ChIP over control 
    # fc = log((binned_counts[[1]][[1]] + pseudoChIP)/(binned_counts[[1]][[2]] + pseudoCtrl))
    # 
    # #compute the normalized enrichment from the foldchange 
    # normEnrichment <- (fc + regul)/stdrz
    # normEnrichment[normEnrichment < 0] = 0
    # normEnrichment[normEnrichment > 1] = 1
    # normEnrichment[is.na(normEnrichment) | is.infinite(normEnrichment)] = 0
    # normEnrichment = normEnrichment*100 #scale for histone density
    # 
    # bins_enrichment = binned_counts$gr
    # mcols(bins_enrichment) = data.frame(score = normEnrichment)
    bins_enrichment = getEnrichmentFine(bamfileChIP, bamfileCtrl, fit, binsize.bigwig, method="coverage", scale.pseudo = scale.pseudo, procs=20)
    export.bw(bins_enrichment, outfile_bw_enrichment, format='bw')
  }
  #plot a deepTools heatmap of normalized enrichment for given regions
  mat_file_enrichment = paste(output_dir_data,'tmp/',title,'_enrichment.mat.gz',sep='')
  cmd = paste(deepTools,'computeMatrix reference-point -S ',outfile_bw_enrichment,' -R ',paste(gene_regions,collapse = " "),' --referencePoint TSS -a ',a,' -b ',b,
              ' --binSize ',binsize.heatmap,' --numberOfProcessors max/2 --outFileName ',mat_file_enrichment,sep='')
  print(cmd)
  system(cmd)
  
  plotFile = paste(output_dir_plot,title,'_enrichment.pdf',sep='')
  if(length(gene_regions) == 1){
    cmd = paste(deepTools,'plotHeatmap -m ',mat_file_enrichment,' -out ',plotFile,
                ' --averageTypeSummaryPlot mean --sortRegions ascend --colorList white,darksalmon,brown --colorNumber 20 --zMin 0',sep='')
  }else{
    cmd = paste(deepTools,'plotHeatmap -m ',mat_file_enrichment,' -out ',plotFile," --regionsLabel ", paste("cluster",1:length(gene_regions),sep = "",collapse = " "),
                ' --averageTypeSummaryPlot mean --sortRegions ascend --colorList white,darksalmon,brown --colorNumber 20 --zMin 0',sep='')
  }
  print(cmd)
  system(cmd)
}

cmd = paste('pdfunite ',output_dir_plot,'*_enrichment.pdf ',output_dir_plot,'heatmaps_enrichment.pdf',sep='')
print(cmd)
system(cmd)

cmd = paste("rm -rf ",output_dir_data,"tmp",sep="")
print(cmd)
system(cmd)
