###################################################################################
#libraries
###################################################################################

library(Cairo)
library(ggplot2)
library(here)
library(reshape2)

###################################################################################
#directories and files
###################################################################################

CAGdelta5_file = here("data/annotation_files/xist_mutants_sakata","GSE93031_allelic_expression.txt")
model_matrix_file = here("data/modelling/feature_matrix/promoter_pro_seq_genes_epigenetic.RData")
clustering_matrix = here("data/modelling/model/silencing_dynamics_model/results_thr_0.5_0.9_1.3","best_clustering_data_set_k3.RData")
output_dir = here("plots/additional_analysis/")

###################################################################################
#load CAGdelta5 data
###################################################################################

table = read.table(file = CAGdelta5_file,sep='\t',header=T)
table_CAGdelta5 = table[,c(1,2,7)]
table_CAGdelta5_chrX = table_CAGdelta5[table_CAGdelta5$chrom=="chrX",]
table_CAGdelta5_chrX = na.omit(table_CAGdelta5_chrX)

#threhold for silenced vs not silenced as defined in Sakata et. al.
thr_silencing = 10

genes_CAGdetla5_silenced = table_CAGdelta5_chrX$gene[table_CAGdelta5_chrX$avePat.CAGdelta5 <= thr_silencing]
genes_CAGdetla5_silenced = as.character(genes_CAGdetla5_silenced)

genes_CAGdetla5_not_silenced = table_CAGdelta5_chrX$gene[table_CAGdelta5_chrX$avePat.CAGdelta5 > thr_silencing]
genes_CAGdetla5_not_silenced = as.character(genes_CAGdetla5_not_silenced)

###################################################################################
#load feature matrix
###################################################################################

load(model_matrix_file)
table_halftimes = data.frame(gene = rownames(data_set), halftime = halftime)

table_halftimes_silenced = table_halftimes[table_halftimes$gene %in% genes_CAGdetla5_silenced,]
table_halftimes_not_silenced = table_halftimes[table_halftimes$gene %in% genes_CAGdetla5_not_silenced,]

###################################################################################
#cumulative distribution of repeat A dependent vs independent genes
###################################################################################

n = nrow(table_halftimes_silenced)
m = nrow(table_halftimes_not_silenced)
table_plot = data.frame(halftime = c(sort(table_halftimes_silenced$halftime),sort(table_halftimes_not_silenced$halftime)),
                   class = c(rep("repeatA independent",n),rep("repeatA dependent",m)))

####plot cumulative distribution function
cairo_pdf(paste(output_dir,'paper_figures_CD_repeatA_sakata.pdf',sep=''),width = 2,height = 3, onefile = TRUE)
ggplot(table_plot, aes(halftime, colour = class)) + 
  stat_ecdf() + 
  scale_x_continuous(breaks=c(0,1,2,3,3.5), label=c("0","1","2","3",">3.5"), name='half-time [days]') +
  scale_y_continuous(breaks=c(0,0.5,1), name='Empirical Cumulative Distribution') +
  scale_color_grey(start=0.4,end=0.8) +
  theme_minimal(base_family = "Source Sans Pro") + 
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),axis.text.x = element_text(size=8), axis.text.y = element_text(size=8), 
        axis.title=element_text(size=8),legend.text = element_text(size=8), legend.title = element_text(size=8), legend.position = "top") +
  guides(fill=guide_legend(nrow=3), col=guide_legend(nrow=3))
dev.off()

###################################################################################
#distribution of halftimes and features of repeat A dependened vs independend genes 
###################################################################################

#plot distribution of halftimes
cairo_pdf(paste(output_dir,'analysis_paper_sakata.pdf',sep=''),width = 2.5,height = 3, onefile = TRUE)
ggplot(table_plot, aes(x=halftime,fill=class)) +
  geom_density(alpha=.4, colour="black",lwd=0.3) + 
  scale_fill_grey(start = 0.1,end=0.8) +
  scale_x_continuous(limits=c(-0.5,4),breaks=c(0,1,2,3,3.5), label=c("0","1","2","3",">3.5"), name='half-time [days]') +
  scale_y_continuous(name='# of genes') +
  theme_minimal(base_family = "Source Sans Pro") + 
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),axis.text=element_text(size=6), axis.title=element_text(size=8), plot.title = element_text(size=8),
        legend.text = element_text(size=8), legend.title = element_text(size=8), legend.position = "bottom") +
  guides(fill=guide_legend(nrow=2), col=guide_legend(nrow=2))

###plot boxplots of features
matrix_CAGdelta5_silenced = data_set[rownames(data_set) %in% genes_CAGdetla5_silenced,]
matrix_not_CAGdelta5_silenced = data_set[rownames(data_set) %in% genes_CAGdetla5_not_silenced,]
n = nrow(matrix_CAGdelta5_silenced)
m = nrow(matrix_not_CAGdelta5_silenced)

for(i in 1:ncol(matrix_CAGdelta5_silenced)){
  feature = colnames(matrix_CAGdelta5_silenced)[i]
  if(!is.factor(matrix_CAGdelta5_silenced[,i])){
    wilcox = wilcox.test(matrix_CAGdelta5_silenced[,i],matrix_not_CAGdelta5_silenced[,i])$p.value
    feature_title = gsub("HiC ","Hi-C ",gsub("number interactions","number",gsub("mean interaction ","",
                     gsub("_"," ",gsub("8WG16","unphosphorylated",gsub('_GLIB_*-?[0-9]*_-?[0-9]*','',gsub('_ENCODE_[A-z]{4}_*-?[0-9]*_-?[0-9]*','',
                     gsub('_GSE[0-9]*_-?[0-9,A-z]*_-?[0-9,A-z]*','',feature,perl=T))))))))
    data_plot = data.frame(repeatA = c(rep("independent",n),rep("dependent",m)), feature = c(matrix_CAGdelta5_silenced[,i],matrix_not_CAGdelta5_silenced[,i]))
    gg_box = ggplot(data_plot, aes(x=repeatA,y=feature)) + 
      geom_boxplot(colour = "#4d4d4d",alpha = 0.7,outlier.size=0.1,lwd=0.4) +
      labs(title=feature_title, subtitle= paste("wilcox test:",signif(wilcox,2))) +
      scale_y_continuous(name = "signal",labels = scales::scientific) +
      theme_minimal(base_family = "Source Sans Pro") + 
      theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),axis.text.x = element_text(size=8, angle = 45, hjust=1), axis.text.y = element_text(size=8), 
            axis.title=element_text(size=8),plot.title = element_text(size=9),plot.subtitle = element_text(size=8)) 
    print(gg_box)
  }
}

dev.off()

###################################################################################
#calculate fration of repeat A dependend and independend genes per cluster
###################################################################################

#select genes that are silenced or not silenced in CAGdelta5 mutant
load(clustering_matrix)
data_set_plot$cluster = as.numeric(data_set_plot$cluster)
cluster_silenced = data_set_plot$cluster[rownames(data_set_plot) %in% genes_CAGdetla5_silenced]
cluster_not_silenced = data_set_plot$cluster[rownames(data_set_plot) %in% genes_CAGdetla5_not_silenced]
total_n_of_genes_in_each_cluster = table(data_set_plot$cluster)

table = rbind(table(factor(cluster_not_silenced,levels = 1:max(data_set_plot$cluster)))/total_n_of_genes_in_each_cluster,
              table(factor(cluster_silenced,levels = 1:max(data_set_plot$cluster)))/total_n_of_genes_in_each_cluster)
table = rbind(1 - colSums(table),table)
rownames(table) = c("not covered", "repeat A dependent genes", "repeat A indenpendent genes")

cairo_pdf(paste(output_dir,'paper_figures_cluster_repeatA_sakata.pdf',sep=''),width = 2,height = 3, onefile = TRUE)
ggplot(aes(x=Var2,y=value,fill=Var1),data=melt(table*100)) +
  geom_bar(stat="identity",color="black") +
  scale_x_continuous(breaks=c(1,2,3), label=c("1", "2", "3"), name='cluster') +
  scale_y_continuous(breaks=c(0,25,50,75,100), name='fraction of genes (%)') +
  scale_fill_grey(name="group",start=0.1,end=0.9) +
  theme_minimal(base_family = "Source Sans Pro") + 
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),axis.text.x = element_text(size=8), axis.text.y = element_text(size=8), 
        axis.title=element_text(size=8),legend.text = element_text(size=8), legend.title = element_text(size=8), legend.position = "bottom") +
  guides(fill=guide_legend(nrow=3), col=guide_legend(nrow=3))
dev.off()

####Fisher test
cluster1 = 2
cluster2 = 1
#repeat A dependent
dependent = table(cluster_not_silenced)
#repeat independent + not covered
independent = total_n_of_genes_in_each_cluster - dependent
#fisher test
matrix = matrix(c(dependent[cluster1],dependent[cluster2],independent[cluster1],independent[cluster2]), nr=2)
fisher.test(matrix)


