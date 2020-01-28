
########################################
#libraries
########################################

library(randomForest)
library(Cairo)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(reshape2)
library(matrixStats)
library(foreach)
library(doParallel)
library(cluster)
library(fpc)
library(corrplot)
library(scales)
library(gridExtra)
library(here)
library(openxlsx)
library(ecoflux) #install.packages("remotes");  remotes::install_github("jpshanno/ecoflux")  --> package to obtain scientific notion in 10 ^ format
source(here("scripts/modelling","model_functions.R"))


########################################
#directories and files
########################################

output_directory_plots = 'plots/modelling'
output_directory_data = 'data/modelling/model'

file_feature_matrix = here("data/modelling/feature_matrix","promoter_pro_seq_genes_epigenetic.RData")
file_high_confidence_escapees = here("data/annotation_files/escapees","escapees.RData")
file_feature_matrix_predictions = here("data/modelling/feature_matrix","promoter_all_genes_epigenetic.RData")
file_gene_annotation = here("data/annotation_files/gene_annotation","gencode.vM9.annotation.chrX.genes.reannotated.with.rr.bed")
file_SNPs = here("data/annotation_files/SNPs","gencode.vM9.annotation.SNP.count.bed")
file_RPKM = here("data/silencing_halftimes/raw_data","GSE121144_PROseq.xlsx")
file_halftimes = here("data/silencing_halftimes/fitted_data","halftimes_pro_seq_mm10_RSS_initial_ratio.txt")
file_results = file(here(output_directory_data,'log_results.txt'), open = "a")

########################################
#setup parameters
########################################

#threshold combination for "silencing dynamics model"
thr_silencing_lower_seq = c(0.5,0.6,0.7)
thr_silencing_middle_seq = c(0.7,0.8,0.9,1)
thr_silencing_upper_seq = c(1,1.1,1.2,1.3,1.4)
class0_label = "early silenced genes"
class1_label = "late silenced genes"

#threshold combination for "XCI/escape" model
thr_silencing_lower_seq = c(0.9,1.0,1.1,1.2,1.3,1.4)
thr_silencing_middle_seq = "-"
thr_silencing_upper_seq = c(1.4,1.5,1.6,1.7,1.8,1.9,2.0)
class0_label = "silenced genes"
class1_label = "not silenced genes"

ncores = 20 # number of cores for parallelization
test_set_min_size = 10 # minimum number of genes in OOB set
ntree = 1000 # number of trees in the random forest 
thr_class_error = 0.05 # maximal difference in error between class 0 and class 1, e.g. if class 0 has error 0.2 then class 1 has to have an error between 0.15 and 0.25
runs = 500 # number of runs for stability test of random forest model
mtry_runs = 100 # number of runs to calculate optimal mtry
B=300 #number of bootstraps for cluster stability analysis
thr_anova_p_value = 0.05 # p-value threshold for anova test
idx_RPKM = c(10,15) #idx_RPKM = c(5) #idx_RPKM = c(5,7)

########################################
#cluster setup for parallel computing
########################################

comb <- function(...) {
  mapply('cbind', ..., SIMPLIFY=FALSE)
}

my_cluster = makeCluster(ncores)
registerDoParallel(my_cluster)

########################################
#run pipeline on all combination of thresholds
########################################

for(thr_silencing_lower in thr_silencing_lower_seq){
  for(thr_silencing_middle in thr_silencing_middle_seq){
    for(thr_silencing_upper in thr_silencing_upper_seq){
      
      ########################################
      #load data
      ########################################
      
      load(file_feature_matrix_predictions)
      data_set_predictions = data_set
      
      data = load_data(file_feature_matrix, thr_silencing_lower, thr_silencing_middle, thr_silencing_upper, file_high_confidence_escapees)
      
      smaller_set = min(c(sum(data[[2]]==0),sum(data[[2]]==1)))
      training = round(min(c(smaller_set/100*80,smaller_set-test_set_min_size)))
      sampsize=c(training,training)
      
      if(training >= 25){
        ########################################
        #setup for new threshold
        ########################################
        
        print(paste('thresholds:',thr_silencing_lower,thr_silencing_middle,thr_silencing_upper))
        
        #new entry in log file
        cat("#########################\n#########################\n#########################", file = file_results, sep = "\n")
        cat(paste('thresholds:',thr_silencing_lower,thr_silencing_middle,thr_silencing_upper), file = file_results, sep = "\n")
        cat(paste('size class 0:',sum(data[[2]]==0)), file = file_results, sep = "\n")
        cat(paste('size class 1:',sum(data[[2]]==1)), file = file_results, sep = "\n")
        cat(paste('sampsize:',sampsize[1]), file = file_results, sep = "\n")
        
        #create output directories and files
        output_directory_data_thr = here(paste(output_directory_data,"/results_thr_",thr_silencing_lower,"_",thr_silencing_middle,"_",thr_silencing_upper,"/",sep=''))
        cmd = paste("mkdir",output_directory_data_thr)
        system(cmd)
        
        output_directory_plots_thr =  here(paste(output_directory_plots,"/results_thr_",thr_silencing_lower,"_",thr_silencing_middle,"_",thr_silencing_upper,"/",sep=''))
        cmd = paste("mkdir",output_directory_plots_thr)
        system(cmd)
        
        #plot boxplots of class 0 vs class 1
        plot_data_boxpots(output_directory_plots_thr, data, class0_label, class1_label)
        
        ########################################
        #train random forest on all features
        ########################################
        
        print("train model on all features...")
        
        print("optimize mtry parameter...")
        mtry_seq = seq(1,ncol(data[[1]]),by=2) #mtry_seq = c(seq(1,46,by=3),seq(50,100,by=8),seq(110,ncol(data[[1]]),by=15))
        opt_mtry = optimize_mtry(data[[1]],data[[2]],mtry_seq,ntree,sampsize,thr_class_error,mtry_runs)
        
        print("Random Forest stability test...")
        random_forest_model_all_features = stability_test(data[[1]], data[[2]], opt_mtry, ntree, sampsize, runs)
        model_error = random_forest_model_all_features[[1]]
        
        #write to log file
        cat(paste('all features: total error:',round(mean(model_error$OOB),2), 
                  '; class 0 error:',round(mean(model_error$class0),2),
                  '; class 1 error:',round(mean(model_error$class1),2)), file = file_results, sep = "\n")
        
        print(paste('all features: total error:',round(mean(model_error$OOB),2), 
                    '; class 0 error:',round(mean(model_error$class0),2),
                    '; class 1 error:',round(mean(model_error$class1),2)))
        
        #save random forest model
        file_random_forest_model_all_features = paste(output_directory_data_thr,'random_forest_model_all_features.RData',sep='')
        save(random_forest_model_all_features,file=file_random_forest_model_all_features)
        
        #plot the results
        plot_feature_importance(output_directory_plots_thr, random_forest_model_all_features, class0_label, class1_label)
        plot_feature_importance_sorted(output_directory_plots_thr, random_forest_model_all_features, class0_label, class1_label)
        
        ########################################
        #optain optimal number of top feature, then train random forest model on top features and do predictions
        ########################################
        
        print("feature selection...")
        
        print("optimize number of top features...")
        top_feature_table = optimize_top_features(data[[1]],data[[2]],ntree,sampsize,runs,random_forest_model_all_features,thr_class_error,mtry_runs)
        x = which.min(top_feature_table[,2])
        print(paste('number of top features:',x))
        
        #select top features in data set
        selected_features0 = random_forest_model_all_features[[5]]
        selected_features1 = random_forest_model_all_features[[6]]
        selected_features = unique(c(selected_features0[1:x],selected_features1[1:x]))
        data_set_selected_features = data[[1]][,colnames(data[[1]])%in%selected_features]
        data_set_predictions_selected_features = data_set_predictions[,colnames(data_set_predictions)%in%selected_features]
        
        print("optimize mtry parameter...")
        mtry_seq = seq(1,ncol(data_set_selected_features),by=1)
        opt_mtry = optimize_mtry(data_set_selected_features,data[[2]],mtry_seq,ntree,sampsize,thr_class_error,mtry_runs)
        
        print("Random Forest stability test, predictions and proximities...")
        random_forest_model_top_features = stability_predictions_proximities(data_set_selected_features, data[[2]], opt_mtry, ntree, sampsize, runs, data_set_predictions_selected_features)
        model_error = random_forest_model_top_features[[1]]
        
        #write to log file
        cat(paste('top features: total error:',round(mean(model_error$OOB),2), 
                  '; class 0 error:',round(mean(model_error$class0),2),
                  '; class 1 error:',round(mean(model_error$class1),2)), 
            file = file_results, sep = "\n")
        cat(paste('number of top features:',x), file = file_results, sep = "\n")
        cat(paste('optimal mtry:',opt_mtry), file = file_results, sep = "\n")
        
        print(paste('top features: total error:',round(mean(model_error$OOB),2), 
                    '; class 0 error:',round(mean(model_error$class0),2),
                    '; class 1 error:',round(mean(model_error$class1),2)))
        
        #plot optimization of top features and error rates
        plot_optimization_top_features(output_directory_plots_thr,top_feature_table,class0_label, class1_label)
        plot_error_rates(output_directory_plots_thr, random_forest_model_all_features, random_forest_model_top_features,class0_label, class1_label)
        
        #predictions on new genes
        predictions = create_predictions_list(random_forest_model_top_features[[2]],data_set_predictions,file_gene_annotation,file_SNPs,file_RPKM,idx_RPKM,file_feature_matrix,file_halftimes)
        predictions_old_genes = predictions[[1]]
        predictions_new_genes = predictions[[2]]
        
        #save random forest model and predictions
        file_random_forest_model_top_features = paste(output_directory_data_thr,'random_forest_model_top_features.RData',sep='')
        save(random_forest_model_top_features,file=file_random_forest_model_top_features)
        
        file_predictions_old_genes = paste(output_directory_data_thr,'random_forest_predictions_old_genes.txt',sep='')
        write.table(predictions_old_genes, file_predictions_old_genes, sep='\t', col.names = T, row.names = F, quote = F)
        
        file_predictions_new_genes = paste(output_directory_data_thr,'random_forest_predictions_new_genes.txt',sep='')
        write.table(predictions_new_genes, file_predictions_new_genes, sep='\t', col.names = T, row.names = F, quote = F)
        
        ########################################
        #perform proximity clustering
        ########################################
        
        print("Proximity clustering...")
        
        #clustering on proximity averaged over multiple random forests
        cat('Proximity clustering on averaged proximities', file = file_results, sep = "\n")
        proximity = random_forest_model_top_features[[3]]
        plot_proximities(output_directory_plots_thr,"average",proximity)
        
        index = find_optimal_k(proximity,data[[2]],B)
        
        clustering_predictions = predictions_old_genes$class[predictions_old_genes$Genes %in% row.names(proximity)]
        names(clustering_predictions) = row.names(proximity)
        proximity_clustering(output_directory_plots_thr,output_directory_data_thr,"average",proximity,index,B,clustering_predictions,data,thr_anova_p_value)
        
        #clustering on proximities of best random forest
        best_rf = get_best_rf(data_set_selected_features,data[[2]],ntree,opt_mtry,sampsize,thr_class_error,runs)
        cat('Proximity clustering on proximities of best Random Forest', file = file_results, sep = "\n")
        proximity = best_rf[[2]]
        plot_proximities(output_directory_plots_thr,"best",proximity)
        
        index = find_optimal_k(proximity,data[[2]],B)
        
        clustering_predictions = best_rf[[1]]
        names(clustering_predictions) = row.names(proximity)
        proximity_clustering(output_directory_plots_thr,output_directory_data_thr,"best",proximity,index,B,clustering_predictions,data,thr_anova_p_value)
      }
    }
  }
}

########################################
#close cluster and file
########################################

close(file_results)
stopCluster(my_cluster)

