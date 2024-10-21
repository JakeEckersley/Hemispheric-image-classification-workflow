# model tuning
library(randomForest)
library(future.apply)
library(parallelly)
library(parallel)


agreement<-function(results){
  accuracy<-sum(results$observed==results$predicted)/nrow(results)
  res_tab<-with(results, table(observed, predicted))
  if("Sun" %in% levels(results$observed) & "Sun" %in% levels(results$predicted)==FALSE){
    res_tab[which(rownames(res_tab)=='Sky'),]<-res_tab[which(rownames(res_tab)=='Sun'),]+res_tab[which(rownames(res_tab)=='Sky'),]
    res_tab<-res_tab[-which(rownames(res_tab)=='Sun'),]
  }
  expected<-list()
  for (i in 1:nrow(res_tab)){
    expected[[i]]<-(as.double(sum(res_tab[i,]))*as.double(sum(res_tab[,i])))/nrow(results)
  }
  expected<-sum(unlist(expected))/nrow(results)
  kappa<-(accuracy-expected)/(1-expected)
  return(as.data.frame(cbind(accuracy,kappa)))
}

### file paths #####################################################################
end_path<-"F:/Training tiffs/"
sample_dataframe<-readRDS(paste0(end_path,'Pixel_value_training_dataset.rds'))
img_details<-read.csv(paste0(end_path,"img_details.csv"))
metadata<-read.csv(paste0(end_path,'img_ID_metadata.csv'))[,c(1,9,12)]


# do this step if you change sample_dataframe
#sample_dataframe<-merge(sample_dataframe,metadata,by='Code',all.x=T,all.y = F)
#
# for (i in 1:nrow(sample_dataframe)){
#   if(sample_dataframe$Sky[i]==''){
#   sample_dataframe$Sky[i]<-img_details$Sky[which(img_details$Code==sample_dataframe$Code[i])]
#   }
# }
# saveRDS(sample_dataframe,paste0(end_path,'Pixel_value_training_dataset.rds'))


### editing sample_dataframe  #####################################################################
sample_dataframe$Class<-as.factor(sample_dataframe$Class)
sample_dataframe$Code<-as.factor(sample_dataframe$Code)
SDF <- split(sample_dataframe, sample_dataframe$Code)
rm(sample_dataframe)


# specify variables to use in the model
all_pixel_vars<-which(colnames(SDF[[1]])%in%c("red","green","blue", "Hue_sRGB","Saturation_sRGB","Value_sRGB","L","A","B","Chroma","Hue","Red_PP","Green_PP","Blue_PP","X","Y","Z"))

# specify variables to use in each colour model
RGB_cols<-which(colnames(SDF[[1]])%in%c("red","green","blue"))
HSV_cols<-which(colnames(SDF[[1]])%in%c( "Hue_sRGB","Saturation_sRGB","Value_sRGB"))
LABCH_cols<-which(colnames(SDF[[1]])%in%c("L","A","B","Chroma","Hue"))
ppRGB_cols<-which(colnames(SDF[[1]])%in%c("Red_PP","Green_PP","Blue_PP"))
XYZ_cols<-which(colnames(SDF[[1]])%in%c("X","Y","Z"))


alt_col_mods<-list(
  all_cols = all_pixel_vars,
  RGB_cols = RGB_cols,
  HSV_cols = HSV_cols,
  LABCH_cols = LABCH_cols,
  ppRGB_cols = ppRGB_cols,
  XYZ_cols = XYZ_cols
)



RF_train<-function(sdf_i,alt_col_mods,node_preds,ntrees){
  #check if the image has data
  if(!nrow(sdf_i)>0){
    return(NULL)
  }
  
  sdf_i$Class<-droplevels(sdf_i$Class)
  Code<-as.character(sdf_i$Code[1])
  
  # check there is a leaf sky and wood category in each image
  if(!length(levels(sdf_i$Class))==3){
    print(paste0('Level missing in ',sdf_i[1,1]))
    print(levels(sdf_i$Class))
    #total_a[[i]]<-1
    return(NULL)
  }
  
  
  ####################################################################
  ######    running test and train procedure for the RF model   ######
  ####################################################################
  
  ###### create groupings for test/train
  set.seed(42)
  
  # number of folds
  k <- 5
  
  # randomly assign test or train
  j<-sample(rep(1:k,each=round(nrow(sdf_i))/k))
  
  
  #splitting test/train dataset
  train <- sdf_i[j!= k, ]
  test <- sdf_i[j == k, ]
  
  results<-list()
  for (j in 1:length(alt_col_mods)){
    col_mod<-alt_col_mods[[j]]
    
    # check the number of node preds doesn't exceed the number of vars
    if(node_preds>length(col_mod)-1){
      results[[j]]<-NULL
    }else{
    
    # train the classifier
    classifier_RF = randomForest(x = train[,col_mod], 
                                 y = train$Class, mtry = node_preds,
                                 ntree = ntrees,importance = T)
    
    #var_importance<-classifier_RF$importance[order(classifier_RF$importance[,4],decreasing=T),drop=F,]
    
    # predict categories of test values
    predicted <-predict(classifier_RF, newdata = test[,col_mod])
    
    # combine test/train in a data frame
    x_rf <- cbind(predicted,test)
    colnames(x_rf)[3]<-'observed'
    col_vars<-names(alt_col_mods)[j]
    sky_condition<-train$Sky[1]
    results[[j]]<-cbind(agreement(x_rf),Code,ntrees,node_preds,col_vars,sky_condition)
    }
  }# col mod loop
  
  # add overall classification accuracy to a list
  return(do.call(rbind,results))
} # image loop

##################################
### w/ parallel
##################################
#trees<-c(seq(50,200,by=50),seq(300,1000,by=100))
trees<-400
#node_vars<-seq(2,length(all_pixel_vars)/2)
node_vars<-4
ntree_acc<-list()
for(p in 1:length(trees)){
  ntrees<-trees[p]
  # loop through different node sizes
  node_pred_acc<-list()
  
  for(kk in 1:length(node_vars)){
    node_preds<-node_vars[kk]

    # Set up plan for parallel processing
    plan(multisession, workers = availableCores(omit = 1))
    
    # Use future_lapply for parallel processing
    mult_acc_results <- future_lapply(SDF, RF_train, alt_col_mods, node_preds, trees[p], future.seed = 42L)
    
    node_pred_acc[[kk]]<-do.call(rbind,acc_results)
    
    print(paste('n_trees =',ntrees,'node_preds =',node_preds))
    
  } #n_preds_node loop
  
  print(paste('n_trees =',ntrees))
  ntree_acc[[p]]<-do.call(rbind,node_pred_acc)
  
} #ntrees loop
optimisation_results<-do.call(rbind,ntree_acc)
saveRDS(optimisation_results,paste0(end_path,'RF_model_optimisation_results.rds'))


optimisation_results<-readRDS(paste0(end_path,'RF_model_optimisation_results.rds'))
#all images
optimisation_results_mean<-aggregate(optimisation_results$accuracy,by=list(optimisation_results$ntrees,optimisation_results$node_preds,optimisation_results$col_vars),FUN=mean)
optimisation_results_mean[order(optimisation_results_mean$x,decreasing = T),]

optimisation_results_mean<-aggregate(optimisation_results$accuracy,by=list(optimisation_results$ntrees,optimisation_results$node_preds,optimisation_results$col_vars,optimisation_results$sky_condition),FUN=mean)







########################################################
### work out whether we can drop a few colour models
########################################################

# Store the names and column indexes in a list
colour_models <- list(
  RGB_cols = RGB_cols,
  HSV_cols = HSV_cols,
  LABCH_cols = LABCH_cols,
  ppRGB_cols = ppRGB_cols,
  XYZ_cols = XYZ_cols
)

# create name from colour model names
get_comb_name <- function(models) {
  # Remove "_cols" from names
  paste0(gsub("_cols", "", comb), collapse = "_")
}

# Get combinations of 2 to 4 colour models
combinations_list <- list()

for (n in 2:5) {
  combs <- combn(names(colour_models), n, simplify = FALSE)
  for (comb in combs) {
    # combine index
    combined_indexes <- unlist(lapply(comb, function(model) colour_models[[model]]))
    # name
    comb_name <- get_comb_name(comb)
    
    # add to list
    combinations_list[[comb_name]] <- combined_indexes
  }
}


alt_col_mods<-c(colour_models,combinations_list)

# Set up plan for parallel processing
plan(multisession, workers = detectCores()-1)


# Use future_lapply for parallel processing
# only use a sub-sample of SDF. computationally expensive.
acc_results_opt <- future_lapply(SDF[sample(0:length(SDF),50)], RF_train, alt_col_mods, 2, 400, future.seed = 42L)

acc_results_mean<-aggregate(acc_results_opt$accuracy,by=list(acc_results_opt$ntrees,acc_results_opt$node_preds,acc_results_opt$col_vars),FUN=mean)
acc_results_mean<-acc_results_mean[order(acc_results_mean$x,decreasing = T),]

acc_results_mean_sky<-aggregate(acc_results_opt$accuracy,by=list(acc_results_opt$ntrees,acc_results_opt$node_preds,acc_results_opt$col_vars,acc_results_opt$sky_condition),FUN=mean)
acc_results_mean_sky<-acc_results_mean_sky[order(acc_results_mean_sky$x,decreasing = T),]
acc_results_mean_sky[which(acc_results_mean_sky$sky_condition=='clear'),]
acc_results_mean_sky[which(acc_results_mean_sky$sky_condition=='overcast'),]
acc_results_mean_sky[which(acc_results_mean_sky$sky_condition=='partly cloudy'),]
