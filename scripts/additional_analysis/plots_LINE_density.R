
library(rtracklayer)

output_dir = "/project/lncrna/Xist/plots/additional_analysis/" 
bedTools = '/home/lisasous/tools/bedtools2/bin/'

#create tmp directory
tmp_dir = paste(output_dir,"tmp/",sep="")
cmd = paste("mkdir",tmp_dir)
system(cmd)

#binning chromosome X
file_mm9_chrom_sizes = '/project/lncrna/Xist/data/annotation_files/mouse_genome/mm9.chrom.sizes'
mm9_chrom_sizes = read.table(file = file_mm9_chrom_sizes, header=F, sep='\t')
colnames(mm9_chrom_sizes) = c('chr','length')

bins = seq(2,mm9_chrom_sizes$length[mm9_chrom_sizes$chr=="chrX"],by=700000)
bins_table = data.frame(chr = rep("chrX",(length(bins)-1)),start = bins[1:(length(bins)-1)],end = bins[2:length(bins)]-1, name = bins[1:(length(bins)-1)]-1, score=0, strand=".")

tmp_bins_file = paste(tmp_dir,'/tmp_bins_file.bed',sep='')
export.bed(bins_table,tmp_bins_file,format='bed')

#LINEs
LINEs_file = '/project/lncrna/Xist/data/annotation_files/LINEs/LINEs_mm9.bed'
LINEs_table = read.table(LINEs_file,sep='\t')
colnames(LINEs_table) = c("chr","start","end","name","score","strand","dispStart","dispEnd","col")
LINEs_table = LINEs_table[LINEs_table$chr == 'chrX',]

tmp_LINEs_file = paste(tmp_dir,'/tmp_LINEs_file.bed',sep='')
export.bed(LINEs_table,tmp_LINEs_file,format='bed')

#calculate LINE density in bins of 500.000bp length
tmp_LINE_density_file = paste(tmp_dir,'/tmp_gene_density.txt',sep='')
cmd = paste(bedTools,'windowBed -a ',tmp_bins_file,' -b ',tmp_LINEs_file,' -w 0 -c > ',tmp_LINE_density_file,sep='') #-w 	Base pairs added upstream and downstream of each entry in A when searching for overlaps in B. Default is 1000 bp.
system(cmd)
LINE_density_table = read.table(tmp_LINE_density_file,header=F,sep='\t')[,c(4,7)]
colnames(LINE_density_table) = c('bin','LINE_density')

cairo_pdf(paste(output_dir,'plots_LINE_density.pdf',sep=''),width = 3,height = 3, onefile = TRUE)
ggplot(LINE_density_table, aes(x=LINE_density)) +
  geom_density(alpha=.01, colour="black",lwd=0.3,bw=0.4) +
  ggtitle("Average LINE density on Chr X") +
  scale_x_continuous(name='number of LINEs') +
  scale_y_continuous(name='density (bins)') +
  theme_minimal(base_family = "Source Sans Pro") + 
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),axis.text=element_text(size=6), axis.title=element_text(size=8), plot.title = element_text(size=10)) 
dev.off()

print(paste("Total number of LINEs on Chr X:",nrow(LINEs_table)))
print(paste("Average LINE density on Chr X:",round(mean(LINE_density_table$LINE_density),2)))
print(paste("Max LINE density on Chr X:",round(max(LINE_density_table$LINE_density),2)))

#remove tmp directory
cmd = paste("rm -rf",tmp_dir)
system(cmd)