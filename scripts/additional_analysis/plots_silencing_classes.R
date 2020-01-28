###################################################################################
#libraries
###################################################################################

library(Cairo)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(here)

###################################################################################
#load data
###################################################################################

table_halftimes = read.table(here('data/silencing_halftimes/fitted_data','halftimes_pro_seq_mm9_reannotated_with_rr.bed'))
halftimes = table_halftimes$V5

early = halftimes[halftimes < 0.5]
silenced = halftimes[halftimes < 0.9]
late = halftimes[halftimes > 0.9 & halftimes < 1.3]
not_silenced = halftimes[halftimes > 1.6]

table_pro_seq = data.frame(halftime = early, silencing_class = rep("early",length(early)), model = rep("silencing dynamics model", length(early)))
table_pro_seq = rbind(table_pro_seq, data.frame(halftime = late, silencing_class = rep("late",length(late)), model = rep("silencing dynamics model", length(late))))
table_pro_seq = rbind(table_pro_seq, data.frame(halftime = silenced, silencing_class = rep("silenced",length(silenced)), model = rep("XCI/escape model", length(silenced))))
table_pro_seq = rbind(table_pro_seq, data.frame(halftime = not_silenced, silencing_class = rep("not silenced",length(not_silenced)), model = rep("XCI/escape model", length(not_silenced))))

###################################################################################
#boxplot of different silncing classes
###################################################################################


cairo_pdf(here('plots/additional_analysis','plots_silencing_classes.pdf'),width = 2,height = 3, onefile = TRUE)
ggplot = ggplot(table_pro_seq, aes(x=silencing_class,y=halftime, fill=model)) + 
  geom_boxplot(colour = "#4d4d4d",alpha = 0.7,outlier.size=-1,lwd=0.4) + 
  ggtitle("Silencing classes based \non PRO-seq") + 
  scale_fill_manual("Model",values=c("#2c5aa0", "#a02c2c")) +
  scale_x_discrete(name = "silencing class",labels=c("early","late","silenced","not silenced")) + 
  scale_y_continuous(breaks=c(0,1,2,3,3.5), label=c("0","1","2","3",">3.5"), name='half-time [days]') +
  theme_minimal(base_family = "Source Sans Pro") + 
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),axis.text.x = element_text(size=8, angle = 45, hjust=1, margin = margin(t=0,b=0)), axis.text.y = element_text(size=8), 
        axis.title=element_text(size=8, margin = margin(t=0)),plot.title = element_text(size=8), legend.text = element_text(size=8), legend.title = element_text(size=8), legend.position = "bottom") +
  guides(fill=guide_legend(nrow=2), col=guide_legend(nrow=2))
legend = get_legend(ggplot)
ggplot = ggplot + theme(legend.position="none")
grid.arrange(ggplot,legend,ncol=1,heights=c(2.5,0.5))
dev.off()



###################################################################################
#plots paper 2e-f)
###################################################################################

####load pro-seq data
table_halftimes = read.table(here('data/silencing_halftimes/fitted_data','halftimes_pro_seq_mm9_reannotated_with_rr.bed'))
halftimes = table_halftimes$V5

early = halftimes[halftimes < 0.5]
silenced = halftimes[halftimes < 0.9]
late = halftimes[halftimes > 0.9 & halftimes < 1.3]
not_silenced = halftimes[halftimes > 1.6]

table_pro_seq = data.frame(halftime = early, silencing_class = rep("early",length(early)), model = rep("silencing dynamics model", length(early)))
table_pro_seq = rbind(table_pro_seq, data.frame(halftime = late, silencing_class = rep("late",length(late)), model = rep("silencing dynamics model", length(late))))
table_pro_seq = rbind(table_pro_seq, data.frame(halftime = silenced, silencing_class = rep("silenced",length(silenced)), model = rep("XCI/escape model", length(silenced))))
table_pro_seq = rbind(table_pro_seq, data.frame(halftime = not_silenced, silencing_class = rep("not silenced",length(not_silenced)), model = rep("XCI/escape model", length(not_silenced))))


####load marks data
load(here('data/modelling/feature_matrix','promoter_pro_seq_genes_epigenetic.RData'))
table_halftimes = data.frame(gene = rownames(data_set), halftime = halftime)

table_marks_paper = read.table(file = here('data/annotation_files/silencing_classes','silencing_classes_marks.txt'),sep='\t',header = F)
colnames(table_marks_paper) = c('gene','silencing_class')

table_marks = merge(table_halftimes,table_marks_paper,by='gene')[,2:3]
levels(table_marks$silencing_class) = c("early","escapee","interm.","late")
table_marks$silencing_class = factor(table_marks$silencing_class,levels = c('early','interm.','late','escapee'),ordered = TRUE)
table_marks$model = "none"
table_marks$source = "Differentiating mESCs"


####load borenzstein
load(here('data/modelling/feature_matrix','promoter_pro_seq_genes_epigenetic.RData'))
table_halftimes = data.frame(gene = rownames(data_set), halftime = halftime)

table_NSMB_paper = read.table(file = here('data/annotation_files/silencing_classes','silencing_classes_borensztein.txt'))
colnames(table_NSMB_paper) = c('gene','silencing_class')

table_boren = merge(table_halftimes,table_NSMB_paper,by='gene')
table_boren = table_boren[table_boren$silencing_class != "Bias",2:3]
levels(table_boren$silencing_class) = c("Bias","early","escapee","interm.","late")
table_boren$silencing_class = factor(table_boren$silencing_class,levels = c('early','interm.','late','escapee'),ordered = TRUE)
table_boren$model = "none"
table_boren$source = "Pre-implantation embryos"

####boxplots
table_pro_seq$source = "PRO-seq in undiff. mESC"
table = rbind(table_marks, table_boren, table_pro_seq)
table$source = factor(table$source, levels = c("Differentiating mESCs","Pre-implantation embryos","PRO-seq in undiff. mESC"))

cairo_pdf(here('plots/additional_analysis','paper_figures_silencing_classes.pdf'),width = 4,height = 3.5, onefile = TRUE)
ggplot(table, aes(x=silencing_class,y=halftime, fill=model)) +
  geom_boxplot(colour = "#4d4d4d",alpha = 0.7,outlier.size=0.1,lwd=0.4) +
  facet_grid(. ~ source, labeller = label_wrap_gen(width = 20, multi_line = TRUE),scales = "free_x") +
  scale_fill_manual("Model",values=c("#2c5aa0","white", "#a02c2c")) +
  scale_x_discrete(name = "silencing class") + 
  scale_y_continuous(breaks=c(0,1,2,3,3.5), label=c("0","1","2","3",">3.5"), name='half-time [days]') +
  theme_minimal(base_family = "Source Sans Pro") + 
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),axis.text.x = element_text(size=8, angle = 45, hjust=1, margin = margin(t=0,b=0)), axis.text.y = element_text(size=8), 
        axis.title=element_text(size=8, margin = margin(t=0)),strip.text = element_text(size=8), legend.text = element_text(size=8), legend.title = element_text(size=8), legend.position = "bottom") +
  guides(fill=guide_legend(nrow=2), col=guide_legend(nrow=2))
dev.off() 
