library(gridExtra)
library(ggplot2)
library(here)

###################################################################################
#load data
###################################################################################

pro_seq_file = here('data/silencing_halftimes/fitted_data','halftimes_pro_seq_mm10.bed')
pro_seq_table = read.table(pro_seq_file,sep='\t',header=F)
colnames(pro_seq_table) = c("Chromosomes","Start","End","Genes","halftime","Strand")

mrna_seq_undiff_file = here('data/silencing_halftimes/fitted_data','halftimes_mrna_seq_undiff_mm10.bed')
mrna_seq_undiff_table = read.table(mrna_seq_undiff_file,sep='\t',header=T)
colnames(mrna_seq_undiff_table) = c("Chromosomes","Start","End","Genes","halftime","Strand")

mrna_seq_diff_file = here('data/silencing_halftimes/fitted_data','halftimes_mrna_seq_diff_mm10.bed')
mrna_seq_diff_table = read.table(mrna_seq_diff_file,sep='\t',header=F)
colnames(mrna_seq_diff_table) = c("Chromosomes","Start","End","Genes","halftime","Strand")

pro_mrna_undiff_table = merge(pro_seq_table,mrna_seq_undiff_table,by='Genes') 
pro_mrna_diff_table = merge(pro_seq_table,mrna_seq_diff_table,by='Genes') 
mrna_undiff_diff_table = merge(mrna_seq_undiff_table,mrna_seq_diff_table,by='Genes') 

###################################################################################
#plotting function
###################################################################################

scatterplot_dense_colors <- function(x1, x2, xlab, ylab){

  df = data.frame(x1,x2)
  
  ## Use densCols() output to get density at each point
  x = densCols(x1,x2, colramp=colorRampPalette(c("black", "white")))
  df$dens = as.factor(col2rgb(x)[1,] + 1L)
  
  ## Map densities to colors
  cols = colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(256)
  cols = colorRampPalette(c("grey","black"))(256)
  df$col = cols[df$dens]
  
  cor = cor.test(x1,x2)
  
  ## Plot it, reordering rows so that densest points are plotted on top
  ggplot = ggplot(data=df[order(df$dens),]) +
    geom_point(aes(x1,x2,color=dens),size=0.5) +
    scale_color_grey(start=0.7,end=0) + 
    geom_smooth(aes(x1,x2),method = "lm", se = FALSE,color="#ff8080", size=0.5) +
    geom_text(aes(x=0,y=3), label=paste('r=',round(cor$estimate,2),"\np=",signif(cor$p.value,2),sep=""), size=2.5, hjust=0, family = "Source Sans Pro") +
    scale_x_continuous(breaks=c(0,1,2,3,3.5), label=c("0","1","2","3",">3.5"), name=xlab) +
    scale_y_continuous(breaks=c(0,1,2,3,3.5), label=c("0","1","2","3",">3.5"), name=ylab) +
    theme_minimal(base_family = "Source Sans Pro") + 
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),axis.text.x = element_text(size=8), axis.text.y = element_text(size=8), 
          axis.title=element_text(size=8), legend.position = "none") 
  
  return(ggplot)
}

###################################################################################
#generate plots
###################################################################################


pro_vs_undiff = scatterplot_dense_colors(pro_mrna_undiff_table$halftime.x,pro_mrna_undiff_table$halftime.y,"PRO-seq undiff.","mRNA-seq undiff.")
pro_vs_diff = scatterplot_dense_colors(pro_mrna_diff_table$halftime.x,pro_mrna_diff_table$halftime.y,'PRO-seq undiff.','mRNA-seq diff.')
diff_vs_undiff = scatterplot_dense_colors(mrna_undiff_diff_table$halftime.x,mrna_undiff_diff_table$halftime.y,'mRNA-seq undiff.','mRNA-seq diff.')

cairo_pdf(here('plots/additional_analysis','paper_figures_mrna_seq_vs_pro_seq.pdf'),width = 2,height = 4.3, onefile = TRUE)
grid.arrange(diff_vs_undiff,pro_vs_undiff,pro_vs_diff,ncol=1)
dev.off()


