###################################################################################
#libraries
###################################################################################

library(Cairo)
library(ggplot2)
library(here)
library(ecoflux) #install.packages("remotes");  remotes::install_github("jpshanno/ecoflux")  --> package to obtain scientific notion in 10 ^ format

output_dir_plot = here('plots/additional_analysis/')

###################################################################################
#PRO-seq
###################################################################################

halftimes_file = here('data/silencing_halftimes/fitted_data','halftimes_pro_seq_mm9_reannotated_with_rr.bed')
halftimes_table = read.table(file = halftimes_file, sep = '\t', header = F)
colnames(halftimes_table) = c('chr','start','end','gene_name','halftime','strand')
halftimes_table$TSS[halftimes_table$strand == '+'] = halftimes_table$start[halftimes_table$strand == '+']
halftimes_table$TSS[halftimes_table$strand == '-'] = halftimes_table$end[halftimes_table$strand == '-']
halftimes_table = halftimes_table[order(halftimes_table$TSS),]

file_mm9_chrom_sizes = here('data/annotation_files/mouse_genome','mm9.chrom.sizes')
mm9_chrom_sizes = read.table(file = file_mm9_chrom_sizes, header=F, sep='\t')
colnames(mm9_chrom_sizes) = c('chr','length')

initiation_sites_file = here('data/annotation_files/xist_initiation_sites','1hr-Xist-Initiation-Sites_mm9.bed')
initiation_sites = read.table(initiation_sites_file)

#plot distribution of halftimes
cairo_pdf(paste(output_dir_plot,'plots_halftime_pro_seq.pdf',sep=''),width = 2,height = 2, onefile = TRUE)
ggplot(halftimes_table, aes(x=halftime)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white", breaks = seq(0,3.5,0.1),lwd=0.3) +
  geom_density(alpha=.4,fill="grey", colour="grey",lwd=0.3) + 
  scale_x_continuous(limits=c(-0.5,4),breaks=c(0,1,2,3,3.5), label=c("0","1","2","3",">3.5"), name='half-time [days]') +
  scale_y_continuous(name='# of genes') +
  theme_minimal(base_family = "Source Sans Pro") + 
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),axis.text=element_text(size=6), axis.title=element_text(size=8), plot.title = element_text(size=8)) 
dev.off()

halftimes_silencing_dynamics = data.frame(halftime=halftimes_table$halftime, silencing_class = factor(3, levels=1:3, labels = c("early", "late", "none")), model="silencing dynamics model")
halftimes_silencing_dynamics$silencing_class[halftimes_silencing_dynamics$halftime <= 0.5] = "early"
halftimes_silencing_dynamics$silencing_class[halftimes_silencing_dynamics$halftime >= 0.9 & halftimes_silencing_dynamics$halftime <= 1.3] = "late"

halftimes_xci_escape = data.frame(halftime=halftimes_table$halftime, silencing_class = factor(3, levels=1:3, labels = c("silenced", "not silenced", "none")), model="XCI/escape model")
halftimes_xci_escape$silencing_class[halftimes_xci_escape$halftime <= 0.9] = "silenced"
halftimes_xci_escape$silencing_class[halftimes_xci_escape$halftime >= 1.6] = "not silenced"

halftimes = rbind(halftimes_xci_escape,halftimes_silencing_dynamics)
halftimes$silencing_class = factor(halftimes$silencing_class, levels=c("silenced", "not silenced", "early", "late", "none"))

cairo_pdf(paste(output_dir_plot,'paper_figures_halftime_distribution_silencing_class.pdf',sep=''),width = 6,height = 4, onefile = TRUE)
ggplot(halftimes, aes(x=halftime, color=silencing_class, fill=silencing_class)) +
  geom_histogram(breaks = seq(0,3.5,0.1),position="identity", alpha=0.6) +
  facet_grid(model ~ ., labeller = label_wrap_gen(width = 20, multi_line = TRUE)) +
  theme_minimal(base_family = "Source Sans Pro") + 
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),axis.text=element_text(size=8), axis.title=element_text(size=8), 
        legend.text = element_text(size=8), legend.title = element_text(size=8), legend.position = "bottom") +
  scale_x_continuous(limits=c(-0.5,4),breaks=c(0,1,2,3,3.5), label=c("0","1","2","3",">3.5"), name='half-time [days]') +
  scale_y_continuous(name='# of genes') +
  scale_color_manual("Silencing class",values=c("#87aade", "#2c5aa0", "#ffaaaa", "#a02c2c", "grey"),label=c("silenced (168)", "not silenced (50)", "early (74 genes)","late (40 genes)", "no class")) +
  scale_fill_manual("Silencing class",values=c("#87aade", "#2c5aa0", "#ffaaaa", "#a02c2c", "grey"),label=c("silenced (168)", "not silenced (50)", "early (74 genes)","late (40 genes)", "no class")) +
  guides(fill=guide_legend(nrow=2), col=guide_legend(nrow=2))
dev.off() 

#plot halftime vs genomic position
cairo_pdf(paste(output_dir_plot,'paper_figures_genomic_position.pdf',sep=''),width = 2.5,height = 3, onefile = TRUE)
ggplot(data=halftimes_table,aes(x=TSS,y=halftime)) +
  geom_vline(xintercept = 100678572, colour="grey") + #TSS of Xist on mm9
  geom_text(aes(x=102000000,y=2.5), label=expression(italic(Xist)~'locus'), size=2.5, hjust=0, color="grey") +
  geom_line(aes(x=halftimes_table$TSS,y=smooth(halftimes_table$halftime),colour="running median")) +
  geom_point(aes(colour="gene"),size=0.5) +
  theme_minimal(base_family = "Source Sans Pro") + 
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), axis.ticks.x=element_line(color="grey"),axis.text=element_text(size=8), axis.title=element_text(size=8), 
        axis.text.x = element_text(hjust=c(0,0.5,1)),legend.text = element_text(size=8),legend.position = "top") +
  scale_x_continuous(limits=c(0,max(halftimes_table$TSS)),breaks=c(0,max(halftimes_table$TSS)/2,max(halftimes_table$TSS)),name='genomic position Chr X [x10⁸]',
                     labels = c(0,round((max(halftimes_table$TSS)/2)/1e+08,2),round((max(halftimes_table$TSS))/1e+08,2))) +
  scale_y_continuous(limits=c(0,3.6),breaks=c(0,1,2,3,3.5), label=c("0","1","2","3",">3.5"), name='half-time [days]') +
  scale_colour_manual(name='', values=c('gene'='#666666', 'running median'='#2c5aa0', guide='legend')) +
  guides(fill = guide_legend(override.aes = list(linetype = 1, shape=''),nrow=3), colour = guide_legend(override.aes = list(linetype=c(0,1), shape=c(16,NA)),nrow=3))
dev.off()

###################################################################################
#mRNA-seq undifferentiated
###################################################################################

mrna_seq_undiff_file = here('data/silencing_halftimes/fitted_data','halftimes_mrna_seq_undiff_mm9_reannotated_with_rr.bed')
mrna_seq_undiff_table = read.table(file = mrna_seq_undiff_file, sep = '\t', header = F)
colnames(mrna_seq_undiff_table) = c('chr','start','end','gene_name','halftime','strand')
mrna_seq_undiff_table$TSS[mrna_seq_undiff_table$strand == '+'] = mrna_seq_undiff_table$start[mrna_seq_undiff_table$strand == '+']
mrna_seq_undiff_table$TSS[mrna_seq_undiff_table$strand == '-'] = mrna_seq_undiff_table$end[mrna_seq_undiff_table$strand == '-']
mrna_seq_undiff_table = mrna_seq_undiff_table[order(mrna_seq_undiff_table$TSS),]

#plot distribution of halftimes
cairo_pdf(paste(output_dir_plot,'plots_halftime_mrna_seq_undiff.pdf',sep=''),width = 2,height = 2, onefile = TRUE)
ggplot(mrna_seq_undiff_table, aes(x=halftime)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white", breaks = seq(0,3.5,0.1),lwd=0.3) +
  geom_density(alpha=.4,fill="grey", colour="grey",lwd=0.3) + 
  scale_x_continuous(limits=c(-0.5,4),breaks=c(0,1,2,3,3.5), label=c("0","1","2","3",">3.5"), name='half-time [days]') +
  scale_y_continuous(name='# of genes') +
  theme_minimal(base_family = "Source Sans Pro") + 
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),axis.text=element_text(size=6), axis.title=element_text(size=8), plot.title = element_text(size=8)) 
dev.off()


###################################################################################
#mRNA-seq differentiated
###################################################################################

mrna_seq_diff_file = here('data/silencing_halftimes/fitted_data','halftimes_mrna_seq_diff_mm9_reannotated_with_rr.bed')
mrna_seq_diff_table = read.table(file = mrna_seq_diff_file, sep = '\t', header = F)
colnames(mrna_seq_diff_table) = c('chr','start','end','gene_name','halftime','strand')
mrna_seq_diff_table$TSS[mrna_seq_diff_table$strand == '+'] = mrna_seq_diff_table$start[mrna_seq_diff_table$strand == '+']
mrna_seq_diff_table$TSS[mrna_seq_diff_table$strand == '-'] = mrna_seq_diff_table$end[mrna_seq_diff_table$strand == '-']
mrna_seq_diff_table = mrna_seq_diff_table[order(mrna_seq_diff_table$TSS),]

#plot distribution of halftimes
cairo_pdf(paste(output_dir_plot,'plots_halftime_mrna_seq_diff.pdf',sep=''),width = 2,height = 2, onefile = TRUE)
ggplot(mrna_seq_diff_table, aes(x=halftime)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white", breaks = seq(0,3.5,0.1),lwd=0.3) +
  geom_density(alpha=.4,fill="grey", colour="grey",lwd=0.3) + 
  scale_x_continuous(limits=c(-0.5,4),breaks=c(0,1,2,3,3.5), label=c("0","1","2","3",">3.5"), name='half-time [days]') +
  scale_y_continuous(name='# of genes') +
  theme_minimal(base_family = "Source Sans Pro") + 
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),axis.text=element_text(size=6), axis.title=element_text(size=8), plot.title = element_text(size=8)) 
dev.off()


