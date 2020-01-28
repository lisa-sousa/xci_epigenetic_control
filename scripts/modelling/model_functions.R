########################################
#this script contains all functions used in model.R
#functions for model building, predictions and proximity clustering
#author: Lisa Barros de Andarde e Sousa
#email: lisa.barros.andrade.sousa@gmail.com
########################################
#Data: Load feature matrix and halftimes
########################################

load_data <- function(feature_matrix_file, thr_silencing_lower, thr_silencing_middle, thr_silencing_upper, high_confidence_escapees_file){
  
  print("load data...")
  
  load(feature_matrix_file)
  
  if(thr_silencing_middle == '-'){
    #XCI/escappe model
    data_set = data_set[halftime < thr_silencing_lower | halftime > thr_silencing_upper,]
    halftime = halftime[halftime < thr_silencing_lower | halftime > thr_silencing_upper]
    target = halftime
    target[target < thr_silencing_lower] = 0
    target[target > thr_silencing_upper] = 1
    target = as.factor(target)
  }else{
    #silencing dynamics model
    data_set = data_set[halftime < thr_silencing_lower | (halftime > thr_silencing_middle & halftime < thr_silencing_upper),]
    halftime = halftime[halftime < thr_silencing_lower | (halftime > thr_silencing_middle & halftime < thr_silencing_upper)]
    target = halftime
    target[target < thr_silencing_lower] = 0 #fastly silenced
    target[target > thr_silencing_middle] = 1 #slowly silenced
    target = as.factor(target)
    
    #exclude high confidence escapees 
    load(high_confidence_escapees_file)
    halftime = halftime[!rownames(data_set) %in% escapees]
    target = target[!rownames(data_set) %in% escapees]
    data_set = data_set[!rownames(data_set) %in% escapees,]
  }
  return(list(data_set, target, halftime))
}

########################################
#Data: Boxplots for classes
########################################

plot_data_boxpots <- function(output_directory_plots_thr, data, class0_label, class1_label){
  
  cairo_pdf(file = paste(output_directory_plots_thr,'data_boxplots.pdf',sep=''), width = 4, height = 4,onefile = T)
  hic_inteactions = c("mean_interaction_strength_HiC_all","mean_interaction_strength_HiC_promoter","mean_interaction_strength_HiC_xist")

  data_set = data[[1]]
  target = data[[2]]
  halftime = data[[3]]
  
  for(i in 1:ncol(data_set)){
    
    column = data_set[,i]
    feature = colnames(data_set)[i]
    
    if(is.factor(column)){
      column_target = cbind.data.frame(column,halftime)
      feature0 = column_target[column_target[,1]==0,] # the feature
      feature1 = column_target[column_target[,1]==1,] # the feature
      
      wilcox = scientific_10x(wilcox.test(feature0[,2],feature1[,2])$p.value)
      
      feature_title = gsub("HiC ","Hi-C ",gsub("number interactions","number",gsub("mean interaction ","",
                           gsub("_"," ",gsub("8WG16","unphosphorylated",gsub('_GLIB_*-?[0-9]*_-?[0-9]*','',gsub('_ENCODE_[A-z]{4}_*-?[0-9]*_-?[0-9]*','',
                           gsub('_GSE[0-9]*_-?[0-9,A-z]*_-?[0-9,A-z]*','',feature,perl=T))))))))
      
      gg_box = ggplot(column_target, aes(x=column,y=halftime)) + 
        geom_boxplot(colour = "#4d4d4d",fill="lightgrey",alpha = 0.7,outlier.size=0.1,lwd=0.4) +
        labs(title=feature_title, subtitle=parse(text=paste0('"p-value = "~',wilcox))) +
        scale_x_discrete(name = "overlap",breaks=c(0,1),labels=c("no","yes")) + 
        scale_y_continuous(breaks=c(0,1,2,3,3.5), label=c("0","1","2","3",">3.5"), name='half-time [days]') + 
        theme_minimal(base_family = "Source Sans Pro") + 
        theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),axis.text.x = element_text(size=8, angle = 45, hjust=1), axis.text.y = element_text(size=8), 
              axis.title=element_text(size=8),plot.title = element_text(size=9),plot.subtitle = element_text(size=8)) 
      print(gg_box)
    }else{
      column_target = cbind.data.frame(column,halftime,target)
      class0 = column_target[column_target$target==0,]
      class1 = column_target[column_target$target==1,]
      
      pears_corr = scientific_10x(cor(column_target[,1],column_target[,2]))
      pears_pValue = scientific_10x(cor.test(column_target[,1],column_target[,2])$p.value)

      wilcox = scientific_10x(wilcox.test(class0[,1],class1[,1])$p.value)

      feature_title = gsub("HiC ","Hi-C ",gsub("number interactions","number",gsub("mean interaction ","",
                                 gsub("_"," ",gsub("8WG16","unphosphorylated",gsub('_GLIB_*-?[0-9]*_-?[0-9]*','',gsub('_ENCODE_[A-z]{4}_*-?[0-9]*_-?[0-9]*','',
                                                            gsub('_GSE[0-9]*_-?[0-9,A-z]*_-?[0-9,A-z]*','',feature,perl=T))))))))
      
      gg_box = ggplot(column_target, aes(x=target,y=column)) + 
        geom_boxplot(colour = "#4d4d4d",fill="lightgrey",alpha = 0.7,outlier.size=0.1,lwd=0.4) +
        labs(title=feature_title, subtitle=parse(text=paste0('"p-value = "~',wilcox))) +
        scale_x_discrete(name = "class",breaks=c(0,1),labels=c(class0_label,class1_label)) + 
        scale_y_continuous(name = "signal",labels = scientific_10x, breaks = scales::pretty_breaks(n=2), limits = c(0,max(1,max(column_target$column)))) +
        theme_minimal(base_family = "Source Sans Pro") + 
        theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),axis.text.x = element_text(size=8, angle = 45, hjust=1), axis.text.y = element_text(size=8), 
              axis.title=element_text(size=8),plot.title = element_text(size=9),plot.subtitle = element_text(size=8)) 
      print(gg_box)
    }
    if(feature %in% hic_inteactions){
      column_target = cbind.data.frame(column,halftime)
      column_target$column[column_target$column > 0] = 1
      column_target$column = as.factor(column_target$column)
      
      feature0 = column_target[column_target[,1]==0,] # the feature
      feature1 = column_target[column_target[,1]==1,] # the feature
      
      wilcox = scientific_10x(wilcox.test(feature0[,2],feature1[,2])$p.value)
      
      feature_title = gsub("HiC ","Hi-C ",gsub("number interactions","number",gsub("mean interaction ","",
                           gsub("_"," ",gsub("8WG16","unphosphorylated",gsub('_GLIB_*-?[0-9]*_-?[0-9]*','',gsub('_ENCODE_[A-z]{4}_*-?[0-9]*_-?[0-9]*','',
                           gsub('_GSE[0-9]*_-?[0-9,A-z]*_-?[0-9,A-z]*','',feature,perl=T))))))))
      
      gg_box = ggplot(column_target, aes(x=column,y=halftime)) + 
        geom_boxplot(colour = "#4d4d4d",fill="lightgrey",alpha = 0.7,outlier.size=0.1,lwd=0.4) +
        labs(title=feature_title, subtitle=parse(text=paste0('"p-value = "~',wilcox))) +
        scale_x_discrete(name = "overlap",breaks=c(0,1),labels=c("no","yes")) + 
        scale_y_continuous(breaks=c(0,1,2,3,3.5), label=c("0","1","2","3",">3.5"), name='half-time [days]') + 
        theme_minimal(base_family = "Source Sans Pro") + 
        theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),axis.text.x = element_text(size=8, angle = 45, hjust=1), axis.text.y = element_text(size=8), 
              axis.title=element_text(size=8),plot.title = element_text(size=9),plot.subtitle = element_text(size=8)) 
      print(gg_box)
      }
  }
  dev.off()
}


########################################
#Random Forest: mtry parameter optimization
########################################

optimize_mtry <- function(data_set,target,mtry_seq,ntree,sampsize,thr_class_error,mtry_runs){
  
  #calculate average error rate (mtry_runs) for each mtry value 
  error_rate_mtry = data.frame(mtry=mtry_seq,oob=0,class1=0,class2=0)
  for(i in 1:length(mtry_seq)){
    error_rate = foreach(j=1:mtry_runs,.combine='comb', .multicombine=TRUE, .packages='randomForest') %dopar% {
      rf = randomForest(data_set, target, type="classification", ntree=ntree, mtry=mtry_seq[i], replace=FALSE, sampsize=sampsize)
      error_rate = rf$err.rate
      return(list(error_rate[ntree,1],error_rate[ntree,2],error_rate[ntree,3]))
    }
    error_rate_mtry[i,2:4] = colMeans(sapply(error_rate, unlist)) #add the mean error rate over "mtry_runs" runs for mtry value
  }
  
  #eliminate all combination where difference of errors between the groups is greater than threshold, butif distance > threshold adapt threshold
  difference_class_error = abs(error_rate_mtry[,3] - error_rate_mtry[,4])
  if((thr_class_error - min(difference_class_error)) < 0){thr_class_error = min(difference_class_error) + 0.01}
  
  error_rate_mtry = error_rate_mtry[difference_class_error < thr_class_error,]#only use mtry values where the differnece between class 1 and class2 error lies within the given threshold
  opt_mtry = error_rate_mtry$mtry[which.min(error_rate_mtry$oob)]
  return(opt_mtry)
}

########################################
#Random Forest: stability test (average error rate and feature importance over x runs)
########################################

sort_features_by_importance <- function(feature_importance, features){
  
  row.names(feature_importance) = features
  medians = rowMedians(as.matrix(feature_importance))
  names(medians) = features
  medians = sort(medians, decreasing=TRUE)
  feature_importance = feature_importance[c(names(sort(medians, decreasing=TRUE))),]
  return(feature_importance)
}

select_features_by_importance <- function(feature_importance, features, importance_threshold){
  medians = rowMedians(as.matrix(feature_importance))
  names(medians) = features
  medians = sort(medians, decreasing=TRUE)
  selected_features = names(medians[medians > importance_threshold])
  return(selected_features)
}

stability_test <- function(data_set, target, mtry, ntree, sampsize, runs){
  
  #build x Random Forests on optimal parameters
  results = foreach(i=1:runs,.combine='comb', .multicombine=TRUE, .packages='randomForest') %dopar% {
    rf = randomForest(data_set, target, type="classification", ntree=ntree, mtry=mtry, replace=FALSE, sampsize=sampsize, importance=TRUE)
    return(c(as.list(as.data.frame(t(rf$err.rate[ntree,]))),as.list(as.data.frame(importance(rf)))))
  }
  
  model_error = data.frame(OOB = round(t(results[[1]])*100,2), class0 = round(t(results[[2]])*100,2), class1 = round(t(results[[3]])*100,2))
  meanImpClass0 = sort_features_by_importance(results[[4]], colnames(data_set))
  meanImpClass1 = sort_features_by_importance(results[[5]], colnames(data_set))
  meanDecGini = sort_features_by_importance(results[[7]], colnames(data_set))
  selected_features0 = select_features_by_importance(results[[4]], colnames(data_set),0)
  selected_features1 = select_features_by_importance(results[[5]], colnames(data_set),0)
 
  return(list(model_error,meanImpClass0,meanImpClass1,meanDecGini,selected_features0,selected_features1))
}

########################################
#Random Forest: top feature optimization
########################################

optimize_top_features <- function(data_set,target,ntree,sampsize,runs,random_forest_model,thr_class_error,mtry_runs){
  
  #calculate error rate for each # of top features (x), maximim of x=20
  selected_features0 = random_forest_model[[5]]
  selected_features1 = random_forest_model[[6]]
  error_table = NULL
  for(x in 2:min(c(length(selected_features0),length(selected_features1),30))){
    data_set_top_x = data_set
    selected_features = unique(c(selected_features0[1:x],selected_features1[1:x]))
    data_set_top_x = data_set_top_x[,colnames(data_set_top_x)%in%selected_features]
    
    mtry_seq = seq(1,ncol(data_set_top_x),by=1)
    opt_mtry = optimize_mtry(data_set_top_x,target,mtry_seq,ntree,sampsize,thr_class_error,mtry_runs)
    model_error = stability_test(data_set_top_x, target, opt_mtry, ntree, sampsize, runs)[[1]]
    
    error_table = rbind(error_table,c(x,round(colMeans(model_error),2)))
  }
  return(error_table)
}

########################################
#Random Forest: stability test (average error rate and feature importance over x runs) with predictions and proximities
########################################

stability_predictions_proximities <- function(data_set, target, mtry, ntree, sampsize, runs, data_set_predictions){
  
  #build x Random Forests on optimal parameters
  results = foreach(i=1:runs,.combine='comb', .multicombine=TRUE, .packages=c('randomForest','flock')) %dopar% {
    #Random Forest model
    rf = randomForest(data_set, target, type="classification", ntree=ntree, mtry=mtry, replace=FALSE, sampsize=sampsize, importance=TRUE, proximity = TRUE)
    #prediction on new genes
    rf_predictions_class = as.numeric(predict(rf,data_set_predictions,type = 'class'))-1
    rf_predictions_vote = predict(rf,data_set_predictions,type = 'vote')[,1] #only save probability for class 0
    predictions = data.frame(class = rf_predictions_class, vote = rf_predictions_vote)
    return(c(as.list(as.data.frame(t(rf$err.rate[ntree,]))),as.list(predictions),rf$proximity))
  }
  
  model_error = data.frame(OOB = round(t(results[[1]])*100,2), class0 = round(t(results[[2]])*100,2), class1 = round(t(results[[3]])*100,2))
  predictions = data.frame(class = rowMeans(results[[4]]),vote = rowMeans(results[[5]]))
  rownames(predictions) = rownames(data_set_predictions)
  proximity = results[6:length(results)]
  proximity = matrix(rowMeans(do.call(rbind, proximity)),nrow=nrow(data_set))
  row.names(proximity) = row.names(data_set)
  colnames(proximity) = row.names(data_set)

  return(list(model_error,predictions,proximity))
}

########################################
#Random Forest: create predictions list
########################################

create_predictions_list <- function(predictions,data_set_predictions,file_gene_annotation,file_SNPs,file_RPKM,idx_RPKM,file_feature_matrix,file_halftimes){
  
  gene_annotation = read.table(file=file_gene_annotation,header=F,sep="\t")[,c(4,2,3,6,5)]
  colnames(gene_annotation) = c("Genes","start","end","strand","biotype")
  gene_annotation$gene_length = gene_annotation$end - gene_annotation$start
  
  predictions$Genes = row.names(predictions)
  
  table_SNPs = read.table(file=file_SNPs, header=F ,sep="\t")[c(4,5)]
  colnames(table_SNPs) = c('Genes','number_of_SPNs')
  
  table_RPKM = read.xlsx(file_RPKM,sheet = 1,na.strings=c("NA","nd"))
  table_RPKM = table_RPKM[table_RPKM$Chromosomes == "chrX",c(1,idx_RPKM)]
  
  #data_set_predictions$Genes = row.names(data_set_predictions)
  
  predictions_list = merge(gene_annotation,predictions,by="Genes")
  predictions_list = merge(predictions_list,table_SNPs,by="Genes")
  predictions_list = merge(predictions_list,table_RPKM,by="Genes")
  #predictions_list = merge(predictions_list,data_set_predictions,by="Genes")
  
  load(file_feature_matrix)
  
  predictions_list = predictions_list[order(predictions_list$vote,decreasing = T),]
  predictions_old_genes = predictions_list[predictions_list$Genes %in% row.names(data_set),]
  predictions_new_genes = predictions_list[!(predictions_list$Genes %in% row.names(data_set)),]
  
  halftimes = read.table(file_halftimes,header = T)
  predictions_old_genes = merge(halftimes,predictions_old_genes,by="Genes")
  predictions_old_genes = predictions_old_genes[,c(1,2,3,4,6,5,8,9,7,14,15,16)]
  return(list(predictions_old_genes,predictions_new_genes))
  #return(predictions_list)
}

########################################
#Random Forest: get best Random Forest for given parameters
########################################

get_best_rf <- function(data_set,target,ntree,mtry,sampsize,thr_class_error,runs){
  
  proximity = NULL
  predictions_class = NULL
  best_error_rate = 1
  
  while (is.null(predictions_class)) {
    for(i in 1:runs){
      rf = randomForest(data_set, target, type="classification", ntree=ntree, mtry=mtry, replace=FALSE, sampsize=sampsize, importance=TRUE, proximity = TRUE)
      error_rate = rf$err.rate[ntree,1]
      class_0_error = rf$err.rate[ntree,2]
      class_1_error = rf$err.rate[ntree,3]
      difference = abs(class_0_error - class_1_error)
      if(error_rate < best_error_rate & difference < thr_class_error){
        best_error_rate = error_rate
        proximity = rf$proximity
        predictions_class = as.numeric(predict(rf,data_set,type = 'class'))-1
      }
    }
    thr_class_error = thr_class_error + 0.05
  }
  names(predictions_class) = rownames(data_set)
  print(paste("Error rate of best random forest: ",round(error_rate,2),round(class_0_error,2),round(class_1_error,2)))
  return(list(predictions_class,proximity))
}

########################################
#Random Forest: plotting functions
########################################

#bxoplot with average feature importance over x Random Forests
plot_feature_importance <- function(output_directory_plots_thr, random_forest_model, class0_label, class1_label){
  
  cairo_pdf(file = paste(output_directory_plots_thr,'feature_importance.pdf',sep=''), width = 5, height = 4,onefile = T)
  
  fill = 'lightgrey'
  line = 'black'
  par(mfrow=c(1,1),mar=c(15,10,10,5),oma=c(5,5,5,5))
  
  ####plot feature importance for class 0
  meanImpClass0 = random_forest_model[[2]]
  selected_features0 = random_forest_model[[5]]
  meanImpClass0_plot = meanImpClass0[selected_features0,]
  
  #remove GEO ID from feature name
  rownames(meanImpClass0_plot) = gsub("HiC ","Hi-C ",gsub("number interactions","number",gsub("mean interaction ","",gsub("_"," ",gsub("8WG16","unphosphorylated",gsub('_GLIB_*-?[0-9]*_-?[0-9]*','',
                                      gsub('_ENCODE_[A-z]{4}_*-?[0-9]*_-?[0-9]*','',gsub('_GSE[0-9]*_-?[0-9,A-z]*_-?[0-9,A-z]*','',rownames(meanImpClass0_plot),perl=T))))))))
  
  n = min(nrow(meanImpClass0_plot),50) #more than 50 features is hard to plot
  gg_box0 = ggplot(melt(t(meanImpClass0_plot)[,n:1]), aes(x=Var2,y=value)) + 
    geom_boxplot(fill="white",colour = "#4d4d4d",alpha = 0.7,outlier.size=0.2,lwd=0.4) +
    ggtitle(paste("feature importance",class0_label)) +
    scale_y_continuous(name = "mean decrease in accuracy",limits=c(-5, 10)) +
    scale_x_discrete(name="") + coord_flip() +
    theme_minimal(base_family = "Source Sans Pro") +
    theme(panel.grid.minor = element_blank(),axis.text.x = element_text(size=8), axis.text.y = element_text(size=6), 
          axis.title=element_text(size=8),plot.title = element_text(size=8,hjust=0.5)) 

  
  ####plot feature importance for class 0
  meanImpClass1 = random_forest_model[[3]]
  selected_features1 = random_forest_model[[6]]
  meanImpClass1_plot = meanImpClass1[selected_features1,]
  
  #remove GEO ID from feature name
  rownames(meanImpClass1_plot) = gsub("HiC ","Hi-C ",gsub("number interactions","number",gsub("mean interaction ","",gsub("_"," ",gsub("8WG16","unphosphorylated",gsub('_GLIB_*-?[0-9]*_-?[0-9]*','',
                                      gsub('_ENCODE_[A-z]{4}_*-?[0-9]*_-?[0-9]*','',gsub('_GSE[0-9]*_-?[0-9,A-z]*_-?[0-9,A-z]*','',rownames(meanImpClass1_plot),perl=T))))))))
  
  n = min(nrow(meanImpClass1_plot),50) #more than 50 features is hard to plot
  gg_box1 = ggplot(melt(t(meanImpClass1_plot)[,n:1]), aes(x=Var2,y=value)) + 
    geom_boxplot(fill="white",colour = "#4d4d4d",alpha = 0.7,outlier.size=0.2,lwd=0.4) +
    ggtitle(paste("feature importance",class1_label)) +
    scale_y_continuous(name = "mean decrease in accuracy",limits=c(-5, 10)) +
    scale_x_discrete(name="") + coord_flip() +
    theme_minimal(base_family = "Source Sans Pro") +
    theme(panel.grid.minor = element_blank(),axis.text.x = element_text(size=8), axis.text.y = element_text(size=6), 
          axis.title=element_text(size=8),plot.title = element_text(size=8,hjust=0.5)) 
  
  grid.arrange(gg_box0,gg_box1,ncol=2,widths=c(4,4))
  
  dev.off()
}

plot_feature_importance_sorted <- function(output_directory_plots_thr, random_forest_model, class0_label, class1_label){
  
  cairo_pdf(file = paste(output_directory_plots_thr,'feature_importance_sorted.pdf',sep=''), width = 9, height = 4,onefile = T)
  
  meanImpClass0 = data.frame(feature=row.names(random_forest_model_all_features[[2]]),meanImpClass0 = rowMedians(random_forest_model_all_features[[2]]))
  meanImpClass0$meanImpClass0[meanImpClass0$meanImpClass0 < 0] = 0
  
  meanImpClass1 = data.frame(feature=row.names(random_forest_model_all_features[[3]]),meanImpClass1 = rowMedians(random_forest_model_all_features[[3]]))
  meanImpClass1$meanImpClass1[meanImpClass1$meanImpClass1 < 0] = 0
  
  importance = merge(meanImpClass0,meanImpClass1,by="feature")
  importance$mean = rowMeans(importance[,2:3])
  importance = importance[order(importance$mean,decreasing = T),]
  row.names(importance) = importance$feature
  importance = as.matrix(importance[1:100,2:3]) #more than 100 features is hard to plot
  
  rownames(importance) = gsub("HiC ","Hi-C ",gsub("number interactions","number",gsub("mean interaction ","",gsub("_"," ",gsub("8WG16","unphosphorylated",gsub('_GLIB_*-?[0-9]*_-?[0-9]*','',
                              gsub('_ENCODE_[A-z]{4}_*-?[0-9]*_-?[0-9]*','',gsub('_GSE[0-9]*_-?[0-9,A-z]*_-?[0-9,A-z]*','',rownames(importance),perl=T))))))))
  colnames(importance) = c(class0_label,class1_label)
  matrix = melt(t(importance))
  matrix$Var1 = factor(matrix$Var1,levels = c(class1_label,class0_label),ordered = TRUE)
  
  # Create the heatmap
  breaks=sort(c(min(importance),0.01,0.2,0.3,0.4,0.5,1,2,3,4,5,10,max(importance)))
  mycol = c("white","#E2E2E2","#DBDCE0","#CFD0DA","#BEC1D4","#ABB0CC","#959CC3","#7D87B9","#6371AF","#4359A7","#023FA5","darkblue")
  
  
  ggplot = ggplot(matrix, aes(Var2, Var1, fill = value)) +
    geom_tile(color="#4d4d4d",size=0.3) +
    theme_minimal(base_family = "Source Sans Pro") +
    theme(axis.text.x = element_text(size = 6,angle=45, hjust=1),axis.text.y = element_text(size = 8), legend.position="bottom", 
          legend.text = element_text(size=7), legend.title = element_text(size=8), axis.title=element_text(size=8),plot.background=element_blank(),panel.border=element_blank()) +
    coord_fixed(ratio=1.5) + 
    scale_fill_gradientn(name="Mean Decrease in accuracy",colours=mycol, values=rescale(breaks), guide="colorbar", limits = c(0,max(matrix$value))) +
    scale_y_discrete(name="",position = "top",expand = c(0, 0)) + 
    scale_x_discrete(name="",expand = c(0, 0)) 
  print(ggplot)
  dev.off()
}

#plot error rate vs number of selected top features
plot_optimization_top_features <- function(output_directory_plots_thr,error_table,class0_label,class1_label){
  
  min = which.min(error_table[,2])
  
  cairo_pdf(file = paste(output_directory_plots_thr,'optimization_top_features.pdf',sep=''), width = 6, height = 4)
  ggplot = ggplot(data=melt(error_table[,2:4]), aes(x=Var1, y=value, shape=Var2, colour=Var2)) +
    geom_point() +
    geom_line() +
    scale_color_brewer(name="Error",labels=c("total",class0_label,class1_label)) +
    scale_shape_ordinal(name="Error",labels=c("total",class0_label,class1_label)) +
    geom_vline(xintercept = min,color="red") +
    scale_y_continuous(name = "error rate (%)",limits = c(0,50)) + 
    scale_x_continuous(name = "number of features",breaks = sort(c(0,5,10,15,20,25,30,min))) +
    theme_minimal(base_family = "Source Sans Pro") +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8),
          legend.text = element_text(size=8), legend.title = element_text(size=8), axis.title=element_text(size=8),plot.background=element_blank(),panel.border=element_blank()) 
  print(ggplot)
  dev.off() 
}

#plot error rate of full model and top features model
plot_error_rates <- function(output_directory_plots_thr, random_forest_model_all_features, random_forest_model_top_features,class0_label,class1_label){
  
  model_error_all = random_forest_model_all_features[[1]]
  model_error_top = random_forest_model_top_features[[1]]
  
  model_error_all$feature_set = 'all_features'
  model_error_top$feature_set = 'top_features'
  model_error = rbind(model_error_all,model_error_top)
  plot_df = melt(model_error,id.vars = c("feature_set"))
  
  
  cairo_pdf(file = paste(output_directory_plots_thr,'error_rate.pdf',sep=''), width = 3, height = 3)
  gg_box = ggplot(plot_df, aes(x=variable,y=value,fill = feature_set)) + 
    geom_boxplot(alpha = 0.7,outlier.size=0.1,lwd=0.4) +
    scale_fill_manual("feature set",values=c("#cccccc", "#87aade"),labels=c("all features","top features")) +
    
    scale_y_continuous(name = "error rate (%)",limits = c(0,50)) + 
    scale_x_discrete(name = "",labels=c("total",class0_label,class1_label)) +

    theme_minimal(base_family = "Source Sans Pro") + 
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),axis.text.x = element_text(size=8, angle = 45, hjust=1), axis.text.y = element_text(size=8), 
          axis.title=element_text(size=8),legend.position = c(0.85, 0.9),legend.text = element_text(size=8),legend.title = element_text(size=8)) 
  print(gg_box)
  dev.off()
}

########################################
#Proximity clustering: convert proximity to distance
########################################

get_distance <- function(proximity){
  distance = as.dist(1 - proximity)
  return(distance)
}

########################################
#Proximity clustering: find optimal k
########################################

clutserboot_pam <- function(distance, B, k) {
  out <- tryCatch(
    {
      # Just to highlight: if you want to use more than one 
      # R expression in the "try" part then you'll have to 
      # use curly brackets.
      # 'tryCatch()' will return the last evaluated expression 
      # in case the "try" part was completed successfully
      
      message("Calculate JI.")
      
      boot_pam = clusterboot(dist(distance), B=B, bootmethod =c("boot"),multipleboot = T, clustermethod = pamkCBI, krange = k, seed = 15555, count=F)
      boot_pam$bootmean
      # The return value of `readLines()` is the actual value 
      # that will be returned in case there is no condition 
      # (e.g. warning or error). 
      # You don't need to state the return value via `return()` as code 
      # in the "try" part is not wrapped insided a function (unlike that
      # for the condition handlers for warnings and error below)
    },
    error=function(cond) {
      message("Error in clusterboot function.")
      message("Here's the original error message:")
      message(cond)
      # Choose a return value in case of error
      return(0)
    },
    warning=function(cond) {
      message("Clusterboot caused warning.")
      message("Here's the original warning message:")
      message(cond)
      # Choose a return value in case of warning
      return(NULL)
    },
    finally={
      # NOTE:
      # Here goes everything that should be executed at the end,
      # regardless of success or error.
      # If you want more than one expression to be executed, then you 
      # need to wrap them in curly brackets ({...}); otherwise you could
      # just have written 'finally=<expression>' 
      message(paste("\nProceed with k=", k+1))
    }
  )    
  return(out)
}

find_optimal_k <- function(proximity,target,B){

  #optimize k
  index = rep(1,length(target)-1)
  n_0 = sum(target==0)#genes with class 0
  n_1 = sum(target==1)#genes with class 1
  diff = abs(n_0-n_1)#difference between class sizes
  
  #convert proximities to distances
  distance = get_distance(proximity)
  
  for(k in 2:(length(target)-1)){
    pam = pam(distance,k,diss=T)
    
    data_set_k = data.frame(cluster = pam$clustering,target=target)
    
    #calculate purity index for each cluster
    p = rep(0,k)
    for(i in 1:k){
      cluster_target_df_i = data_set_k[data_set_k$cluster == i,]
      #upscale smaller class
      if(n_0 < n_1){
        x_0 = sum(cluster_target_df_i$target == 0)
        upscaled_x_0 = x_0 + x_0/n_0 * diff
        x_1 = sum(cluster_target_df_i$target == 1)
        n_i = upscaled_x_0 + x_1
        p[i] = (upscaled_x_0/n_i) * (x_1/n_i)
      }else{
        x_0 = sum(cluster_target_df_i$target == 0)
        x_1 = sum(cluster_target_df_i$target == 1)
        upscaled_x_1 = x_1 + x_1/n_1 * diff
        n_i = x_0 + upscaled_x_1
        p[i] = (x_0/n_i) * (upscaled_x_1/n_i)
      }
    }
    purity = sum(p)/k
    
    #filter for stable clusterings
    JI = clutserboot_pam(distance, B, k)
    
    if(sum(JI < 0.6) == 0){index[k] = 4*purity}else{index[k] = 1}
    
    if(k > 4){
      if(sum(index[(k-3):k])==4){break}
    }
  }
  return(index)
}

########################################
#Proximity clustering: plot clustering
########################################

proximity_clustering <- function(output_directory_plots_thr,output_directory_data_thr,title,proximity,index,B,predictions,data,thr_anova_p_value){

  #merge clustering with data set, target and predictions
  data_set = cbind(predictions,target=data[[2]],halftime=data[[3]],data[[1]])
  for(i in 1:ncol(data_set)){
    if(is.factor(data_set[,i])){data_set[,i] = as.numeric(data_set[,i])-1}
  }
  
  #convert proximities to distances
  distance = get_distance(proximity)
  
  cat(paste("tree: ",title), file = file_results, sep = "\n")
  
  for(k in 2:6){ #max number of clusters for plotting: 6, exclude k=1
    if(index[k] != 1 | (index[k] == 1 & k==2) | (index[k] == 1 & k==3)){
      boot_pam = clusterboot(dist(distance), B=B, bootmethod =c("boot"),multipleboot = T, clustermethod = pamkCBI, krange = k, seed = 15555, count=F)
      
      #sort clusters by number of class 0/1 genes in cluster
      sort_cluster = data.frame(cluster = boot_pam$result$partition,prediction = data_set$predictions)
      class1_genes = data.frame(cluster = 1:k,class1_genes = rep(0,k))
      for(i in 1:k){class1_genes$class1_genes[i]=sum(sort_cluster$prediction[sort_cluster$cluster==i])}
      class1_genes = class1_genes[order(class1_genes$class1_genes),]
      class1_genes$new_cluster = 1:k
      for(i in 1:k){boot_pam$result$partition[boot_pam$result$partition==i] = class1_genes$new_cluster[class1_genes$cluster==i]+k}
      boot_pam$result$partition = boot_pam$result$partition-k
      
      #cluster stability
      plot_clustering_stability(output_directory_plots_thr,title,k,boot_pam,class1_genes)
      
      cat(paste("k: ",k), file = file_results, sep = "\n")
      cat(paste("index: ",round(index[k],2)), file = file_results, sep = "\n")
      cat(paste(paste("cluster",1:k),round(boot_pam$bootmean[class1_genes$cluster],4), sep=": "), file = file_results, sep = "\n")
      
      data_set_plot = cbind(cluster=boot_pam$result$partition,data_set)
      data_set_plot = data_set_plot[order(data_set_plot$cluster,data_set_plot$predictions,data_set_plot$halftime),]
      data_set_plot$cluster = as.factor(data_set_plot$cluster)
      
      #save data set (for futher analysis)
      file_cluster_k = paste(output_directory_data_thr,title,'_clustering_data_set_k',k,'.RData',sep='')
      save(data_set_plot,file=file_cluster_k)
      
      #select significant features according to anova or wilcox (in case of two clusters)
      p_value_table = data.frame(feature = colnames(data_set_plot)[4:ncol(data_set_plot)], p_value = rep(-1,(ncol(data_set_plot)-3)))
      for(i in 4:ncol(data_set_plot)){
        data_anova = data_set_plot[,c(1,i)]
        colnames(data_anova) = c('cluster','feature')
        if(length(unique(data_set_plot$cluster))==2){
          #wilcox test
          p_value = wilcox.test(data_anova$feature[data_anova$cluster==1],data_anova$feature[data_anova$cluster==2])$p.value
        }else{
          #anova analysis
          feature.aov = aov(feature ~ cluster, data = data_anova)
          p_value = unlist(summary(feature.aov))["Pr(>F)1"]
        }
        p_value_table$p_value[i-3] = p_value
      }
      
      p_value_table = p_value_table[order(p_value_table$p_value),]
      data_set_plot_sorted = data_set_plot[as.character(p_value_table$feature)]#order features by significance
      
      #plot boxplots
      plot_clustering_anova(output_directory_plots_thr,title,k,data_set_plot,data_set_plot_sorted,p_value_table)
      
      #plot heatmap
      data_set_norm = cbind(data_set_plot[1:3],data_set_plot_sorted)
      title_heatmap = paste("clustering_heatmap_k",k,sep="")
      plot_clustering_heatmap(output_directory_plots_thr,k,title,data_set_norm,title_heatmap)
      
      ###plot only significant features in heatmap
      significant_features = p_value_table$feature[p_value_table$p_value < thr_anova_p_value][1:50] #more than 50 features is hard to plot
      data_set_norm = cbind(data_set_plot[1:3],data_set_plot_sorted[colnames(data_set_plot_sorted) %in% significant_features])
      title_heatmap = paste("clustering_heatmap_significant_features_k",k,sep="")
      plot_clustering_heatmap(output_directory_plots_thr,k,title,data_set_norm,title_heatmap)
    }
  } 
}

########################################
#Proximity clustering: plotting fucntions
########################################

#heatmap with proximities
plot_proximities <- function(output_directory_plots_thr,title,proximity){
 
  cairo_pdf(file = paste(output_directory_plots_thr,title,'_proximity.pdf',sep=''), width = 10, height = 10)
  par(mar=c(1,4,1,0),oma=c(0,2,0,0)) 
  
  corrplot(proximity, order='hclust', hclust.method='ward.D2', col=colorRampPalette(c('blue','white','red'))(200), diag=F, cl.cex=1, tl.cex=0.3, tl.col = "black",is.corr = FALSE,
           bg='white', addgrid.col=NA)
  
  dev.off()
}

#plot with stability of each cluster
plot_clustering_stability <- function(output_directory_plots_thr,title,k,boot_pam,class1_genes){

  cairo_pdf(file = paste(output_directory_plots_thr,title,'_clustering_stability_k',k,'.pdf',sep=''), width = 4, height = k*1.5)
  
  boot_pam$bootresult = boot_pam$bootresult[class1_genes$cluster,]
  boot_pam$bootmean = boot_pam$bootmean[class1_genes$cluster]
  boot_plot = melt(boot_pam$bootresult)
  boot_plot$Var1 = paste("cluster", boot_plot$Var1)
  
  text <- data.frame(label = paste("JS =",round(boot_pam$bootmean,3)), Var1 = unique(boot_plot$Var1))
  ggplot = ggplot(boot_plot, aes(x=value)) +
    geom_histogram(fill = "grey",color="darkgrey",breaks = seq(0,1,0.01),position="identity", alpha=0.6) +
    facet_grid(Var1 ~ ., labeller = label_wrap_gen(width = 20, multi_line = TRUE)) +
    theme_minimal(base_family = "Source Sans Pro") + 
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),axis.text=element_text(size=8), axis.title=element_text(size=8),
          strip.text.y = element_text(size = 8),plot.title = element_text(size=10, hjust = 0.5)) +
    scale_x_continuous(limits=c(0,1),breaks=c(0,0.5,1), name='Jaccard Similarity') +
    scale_y_continuous(name='Frequency') +
    ggtitle("cluster stability") +
    geom_text(data = text, mapping = aes(x = -Inf, y = Inf, label = label),size=3,hjust=-0.5,vjust=3)
  print(ggplot)
  dev.off()
  
}

#boxplots of feature distribution per cluster sorted by anova p-value
plot_clustering_anova <- function(output_directory_plots_thr,title,k,data_set_plot,data_set_plot_sorted,p_value_table){

  cairo_pdf(file = paste(output_directory_plots_thr,title,'_clustering_boxplots_k',k,'.pdf',sep=''), width = 4, height = 4, onefile = T)
  cluster = as.factor(data_set_plot$cluster)
  
  for(i in 1:ncol(data_set_plot_sorted)){
    feature = colnames(data_set_plot_sorted)[i]
    data_plot = data.frame(cluster = cluster,feature = data_set_plot_sorted[,i])
    
    feature_title = gsub("HiC ","Hi-C ",gsub("number interactions","number",gsub("mean interaction ","",
                         gsub("_"," ",gsub("8WG16","unphosphorylated",gsub('_GLIB_*-?[0-9]*_-?[0-9]*','',gsub('_ENCODE_[A-z]{4}_*-?[0-9]*_-?[0-9]*','',
                         gsub('_GSE[0-9]*_-?[0-9,A-z]*_-?[0-9,A-z]*','',feature,perl=T))))))))
    
    gg_box = ggplot(data_plot, aes(x=cluster,y=feature)) + 
      geom_boxplot(colour = "#4d4d4d",fill="lightgrey",alpha = 0.7,outlier.size=0.1,lwd=0.4) +
      labs(title=feature_title, subtitle= parse(text=paste0('"p-value = "~',scientific_10x(p_value_table$p_value[i])))) +
      scale_x_discrete(name = "cluster") + 
      scale_y_continuous(name = "signal",labels = scientific_10x, breaks = scales::pretty_breaks(n=2)) +
      theme_minimal(base_family = "Source Sans Pro") + 
      theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),axis.text.x = element_text(size=8, hjust=1), axis.text.y = element_text(size=8), 
            axis.title=element_text(size=8),plot.title = element_text(size=9),plot.subtitle = element_text(size=8)) 
    print(gg_box)
  }
  dev.off() 
}

#plot of feature heatmap, features sorted by anova p-value and cluster
plot_clustering_heatmap <- function(output_directory_plots_thr,k,title,data_set_norm,title_heatmap){
  
  colnames(data_set_norm) = gsub("HiC ","Hi-C ",gsub("number interactions","number",gsub("mean interaction ","",
                       gsub("_"," ",gsub("8WG16","unphosphorylated",gsub('_GLIB_*-?[0-9]*_-?[0-9]*','',gsub('_ENCODE_[A-z]{4}_*-?[0-9]*_-?[0-9]*','',
                       gsub('_GSE[0-9]*_-?[0-9,A-z]*_-?[0-9,A-z]*','',colnames(data_set_norm),perl=T))))))))
  
  #z-score transformation of data set
  for(i in 4:(ncol(data_set_norm))){
    if(!is.factor(data_set_norm[,i])){
      data_set_norm[,i] = (data_set_norm[,i] - mean(data_set_norm[,i])) / sd(data_set_norm[,i])
    }
  }
  data_set_norm = replace(data_set_norm, data_set_norm > 3, 3)
  data_set_norm = replace(data_set_norm, data_set_norm < -3, 3)
  
  data_set_norm$predictions[data_set_norm$predictions == 0] = -3
  data_set_norm$predictions[data_set_norm$predictions == 1] = 3
  
  data_set_norm$target[data_set_norm$target == 0] = -3
  data_set_norm$target[data_set_norm$target == 1] = 3
  data_set_norm$cluster = as.numeric(data_set_norm$cluster)
  
  matrix = melt(t(data_set_norm))
  cairo_pdf(file = paste(output_directory_plots_thr,title,"_",title_heatmap,'.pdf',sep=''), width = 8, height = 8, onefile = T)
  # Create the heatmap
  ggplot = ggplot(matrix, aes(Var2, Var1, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-3,3)) + # Change gradient color
    theme_minimal(base_family = "Source Sans Pro",base_size=8) +
    theme(axis.text.x = element_blank(),axis.text.y = element_text(size = 8),legend.text = element_text(size=8), 
          legend.title = element_text(size=8), axis.title=element_text(size=8),plot.background=element_blank(),panel.border=element_blank()) +
    coord_fixed(ratio=4) + 
    scale_y_discrete(name="feature",limits = rev(levels(matrix$Var1))) +
    scale_x_discrete(name="gene promoters") 
  
  print(ggplot)  
  dev.off()
}

