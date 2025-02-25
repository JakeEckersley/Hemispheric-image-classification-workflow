
# imagewise RF model training and application
### Libraries + functions #####################################################################
library(randomForest)
library(terra)
library(scales)
library(dismo)
library(gtools)
library(sf)
library("segmented")
source('LAI to PAI functions.r')

### file paths #####################################################################

end_path<-"F:/Training tiffs/"
sample_dataframe<-readRDS(paste0(end_path,'Pixel_value_training_dataset_v4.rds'))
img_details<-read.csv(paste0(end_path,"img_details.csv"))
metadata<-read.csv(paste0(end_path,'img_ID_metadata.csv'))[,c(1,9,12)]
all_colourmodels<-c(Sys.glob(file.path("F:/PhotoproRGB_LAB_imgs",'*.tif')),Sys.glob(file.path("F:/PhotoproRGB_LAB_imgs",'*.tiff')))
updated_shapes="F:/Training tiffs/edited_shapes/"
shapefiles<-Sys.glob(file.path(updated_shapes, "*.shp"))

# re-do this step if you change sample_dataframe
# sample_dataframe<-merge(sample_dataframe,metadata,by='Code',all.x=T,all.y = F)
# 
# for (i in 1:nrow(sample_dataframe)){
#   if(sample_dataframe$Sky[i]==''){
#     sample_dataframe$Sky[i]<-img_details$Sky[which(img_details$Code==sample_dataframe$Code[i])]
#   }
# }
# saveRDS(sample_dataframe,paste0(end_path,'Pixel_value_training_dataset.rds'))


### editing sample_dataframe  #####################################################################
sample_dataframe$Class<-as.factor(sample_dataframe$Class)
sample_dataframe$Code<-as.factor(sample_dataframe$Code)
SDF <- split(sample_dataframe, sample_dataframe$Code)
rm(sample_dataframe)


###################################################################################################
### individual image RF classifier ################################################################
###################################################################################################
# acc_results<-list()
# tables<-list()
# results<-list()
# 
# # these are brand new.
# class_by_seg<-list()
# img_a<-list()
# class_by_ring<-list()
# class_by_seg_4LAI<-list()
# img_a_4LAI<-list()
# class_by_seg_4LAI2200<-list()
# class_by_seg_4LAI60deg<-list()
# class_by_seg_57deg<-list()
# var_importance<-list()

###################################################################################################
### first half: model training     ################################################################
###################################################################################################
#classify_img<-function(SDF_i,end_path,all_colourmodels,Code,classifier_RF){}

#tables<-readRDS(paste0(end_path,'missclass_tables_v4.rds'))
#acc_results<-readRDS(paste0(end_path,'individual_RF_global_accuracy_results_v4.rds'))
#acc_results<-split(acc_results,acc_results$Code)
#results<-readRDS(paste0(end_path,'individual_RF_sample_prediction_accuracy_results_v4.rds'))
samp_list<-Sys.glob(file.path(paste0(end_path,"Sample locations"),'*.shp'))


# redo=T
# re-do everything below 72 nothing in 'best' list for 106 down.
# need to do the above
# best<-list()

# currently running 1 to 106
# running everything >0.85 tomorrow
# class_by_ring<-readRDS(paste0(end_path,'segmented_byzenith__45zen_noazm_maxlikelihood_v5.rds'))
# class_by_ring<-split(class_by_ring,class_by_ring$Code)
# 
# img_a<-as.data.frame(readRDS(paste0(end_path,'Total_a_by_img_zen0to90_v5.rds')))
# img_a<-split(img_a,img_a$Code)
# 
# img_a_4LAI<-as.data.frame(readRDS(paste0(end_path,'Total_a_by_img_zen0to70_v5.rds')))
# img_a_4LAI<-split(img_a_4LAI,img_a_4LAI$Code)
# 
# class_by_seg_4LAI<-readRDS(paste0(end_path,'segmented_7zen_24azm_maxlikelihood_v5.rds'))
# class_by_seg_4LAI<-split(class_by_seg_4LAI,class_by_seg_4LAI$Code)
# 
# class_by_seg_4LAI2200<-readRDS(paste0(end_path,'segmented_5zen_8azm_maxlikelihood_v5.rds'))
# class_by_seg_4LAI2200<-split(class_by_seg_4LAI2200,class_by_seg_4LAI2200$Code)
# 
# class_by_seg_57deg<-readRDS(paste0(end_path,'segmented_zen57deg_8azm_maxlikelihood_v5.rds'))
# class_by_seg_57deg<-split(class_by_seg_57deg,class_by_seg_57deg$Code)
# 
# acc_results<-readRDS(paste0(end_path,'individual_RF_global_accuracy_results__v5.rds'))
# #acc_results<-split(acc_results,acc_results$Code)
# 
# tables<-readRDS(paste0(end_path,'missclass_tables_v5.rds'))
# #tables<-split(tables,tables$Code)
# 
# results<-readRDS(paste0(end_path,'individual_RF_sample_prediction_accuracy_results_v5.rds'))
# names(results)<-lapply(results,function(x)unique(x[,which(colnames(x)=='Code')]))
# #results<-split(results,results$Code)
# 
# var_importance<-readRDS(paste0(end_path,'colour_model_band_importance_v5.rds'))
# var_importance<-split(var_importance,var_importance$Code)
# 
# best<-readRDS(paste0(end_path,'best_model_v5.rds'))
# #best<-split(best,best$Code)
# 
# class_by_seg<-readRDS(paste0(end_path,'segmented_45zen_8azm_maxlikelihood_v5.rds'))
# class_by_seg<-split(class_by_seg,class_by_seg$Code)
# 




class_by_ring<-readRDS(paste0(end_path,'segmented_byzenith__45zen_noazm_maxlikelihood_v6.rds'))
img_a<-readRDS(paste0(end_path,'Total_a_by_img_zen0to90_v6.rds'))
img_a_4LAI<-readRDS(paste0(end_path,'Total_a_by_img_zen0to70_v6.rds'))
class_by_seg_4LAI<-readRDS(paste0(end_path,'segmented_7zen_24azm_maxlikelihood_v6.rds'))
class_by_seg_4LAI2200<-readRDS(paste0(end_path,'segmented_5zen_8azm_maxlikelihood_v6.rds'))
class_by_seg_57deg<-readRDS(paste0(end_path,'segmented_zen57deg_8azm_maxlikelihood_v6.rds'))
acc_results<-readRDS(paste0(end_path,'individual_RF_global_accuracy_results_v6.rds'))
tables<-readRDS(paste0(end_path,'missclass_tables_v6.rds'))
results<-readRDS(paste0(end_path,'individual_RF_sample_prediction_accuracy_results_v6.rds'))
var_importance<-readRDS(paste0(end_path,'colour_model_band_importance_v6.rds'))
best<-readRDS(paste0(end_path,'best_model_v6.rds'))
class_by_seg<-readRDS(paste0(end_path,'segmented_45zen_8azm_maxlikelihood_v6.rds'))

#double check if there's anything to be done with the absolute worst ones.
acc_order<-do.call(rbind,acc_results)
acc_order<-acc_order[order(acc_order$accuracy),]
acc_order<-acc_order[1:(nrow(acc_order)/20),]
redos<-as.vector(acc_order$Code)

# check if their image hasn't been classified
# classified.tifs<-Sys.glob(file.path(paste0(end_path,'classified_imgs_17Nov24/', "*.tif")))
# ims<-gsub('_individ_seg_classified','',tools::file_path_sans_ext(basename(classified.tifs)))
# redos<-names(acc_results)[which(!names(acc_results)%in%ims)]

#check which ones are good but haven't been segmented
#redos<-names(acc_results)[which(!names(acc_results)[names(acc_results)%in%unlist(lapply(acc_results,function(x)if(x[1,1]>0.85){return(x[1,2])}else{return(NULL)}))]%in%names(class_by_seg_4LAI))]
# last thing is to re-do the whichmax image for whatever is the worst acc_result(to recreate the figure....)
wtf<-list()
for(i in 1:length(tables)){
  if(sum(tables[[i]])<200){wtf[[i]]<-i}else{wtf[[i]]<-NULL}
}
redos<-names(SDF)[which(names(SDF)%in%names(tables)[do.call(rbind,wtf)])]
best[which(names(SDF)%in%redos)]

for(i in 46:length(SDF)){
  start<-Sys.time()
  if(i%%10==0){
    saveRDS(class_by_ring,paste0(end_path,'segmented_byzenith__45zen_noazm_maxlikelihood_v7.rds'))
    saveRDS(img_a,paste0(end_path,'Total_a_by_img_zen0to90_v7.rds'))
    saveRDS(img_a_4LAI,paste0(end_path,'Total_a_by_img_zen0to70_v7.rds'))
    saveRDS(class_by_seg_4LAI,paste0(end_path,'segmented_7zen_24azm_maxlikelihood_v7.rds'))
    saveRDS(class_by_seg_4LAI2200,paste0(end_path,'segmented_5zen_8azm_maxlikelihood_v7.rds'))
    saveRDS(class_by_seg_57deg,paste0(end_path,'segmented_zen57deg_8azm_maxlikelihood_v7.rds'))
    saveRDS(acc_results,paste0(end_path,'individual_RF_global_accuracy_results_v7.rds'))
    saveRDS(tables,paste0(end_path,'missclass_tables_v7.rds'))
    saveRDS(results,paste0(end_path,'individual_RF_sample_prediction_accuracy_results_v7.rds'))
    saveRDS(var_importance,paste0(end_path,'colour_model_band_importance_v7.rds'))
    saveRDS(best,paste0(end_path,'best_model_v7.rds'))
    saveRDS(class_by_seg,paste0(end_path,'segmented_45zen_8azm_maxlikelihood_v7.rds'))
  }
  ######################### list of checks
  #########################################
  
  
  if(!names(SDF)[i]%in%redos){next}#else{best[[i]]<-0}
  
  #if(Code=='DSC04736'){next}
  
  if(nrow(SDF[[i]])>0){
    #split off each image dataset
    sdf_i<-SDF[[i]]
    sdf_i$Class<-droplevels(sdf_i$Class)
    Code<-names(SDF)[i]
  }else{
    next
  }
  
  acidx<-which(names(acc_results)==Code)
  
  if(length(acidx)==0){
    print(paste0(Code,' not in accc'))
    next
    }
  
  # if(acc_results[[acidx]][1,1]>0.85&(!names(acc_results)[acidx]%in%redos)){
  #   print(paste0(Code,' accuracy too high = ',acc_results[[acidx]][1,1]))
  #   next
  # }else{
    samps_all<-terra::vect(samp_list[which(gsub('_SampleLocations','',tools::file_path_sans_ext(basename(samp_list)))==Code)])
  #}
  
  if(i>length(best)){
    best[[i]]<-0
  }
  
  if(is.null(best[[i]])){
    best[[i]]<-0
    best_start<-0
  }else{
    best_start<-best[[i]]
    if(best[[i]]==5){next}
  }
  #if(!file.exists(paste0(end_path,"georef_classified_imgs_12Nov24/",Code,"_individ_seg_classified.tif"))){
    
  # check there is a leaf sky and wood category in each image
  if(!length(levels(sdf_i$Class))==3){
    print(paste0('Level missing in ',sdf_i[1,1]))
    print(levels(sdf_i$Class))
    #total_a[[i]]<-1
    next
  }

  ####################################################################
  ######    train first RF model   ######
  ####################################################################
  
  ###### create groupings for test/train
  # number of folds
  k <- 4
  
  # randomly assign test or train
  set.seed(42)
  j<-sample(rep(1:k,each=round(nrow(sdf_i))/k))
  
  # specify variables to use in the model
  pixel_vars<-which(colnames(sdf_i)%in%c('L','A','B','Hue','Chroma',"red","green","blue", "Hue_sRGB","Saturation_sRGB","Value_sRGB"))
  sdf_i$Class<-class_to_numeric(sdf_i$Class)
  
  
  #splitting test/train dataset
  train <- sdf_i[j != k, ]
  test <- sdf_i[j == k, ]
  
  # train the classifier
  classifier_RF_1 = randomForest(x = train[,pixel_vars], 
                               y = train$Class, 
                               ntree = 200,mtry=2,importance = T) 
  
  ## add something that captures 'importance'
  idx<-which(names(var_importance)==Code)
  if(length(idx)==0){idx<-length(var_importance)+1}
  var_importance[[idx]]<-as.data.frame(classifier_RF_1$importance[order(classifier_RF_1$importance[,4],decreasing=T),drop=F,])
  var_importance[[idx]]$Code<-Code
  var_importance[[idx]]$colmod<-rownames(var_importance[[idx]])
  colnames(var_importance[[idx]])[1:3]<-c('Sky','Leaves','Stems')



  ##################################################
  ### second half: model application ###############
  ### Read in the image       ######################
  ##################################################
  indx<-which(gsub("_lab","",tools::file_path_sans_ext(basename(all_colourmodels)))==Code)
  srgb_indx<-which(gsub("_srgb","",tools::file_path_sans_ext(basename(all_colourmodels)))==Code)

  if(length(indx)>0){# check there is an image to classify....
    print('image found')
    }else{
    print(paste0('error for i = ',i,' - no file in rgb_tiffs named ',Code))
  }
    
  #read in LAB image
  LAB_tiff<-rast(all_colourmodels[[indx]])
  srgb_tiff<-rast(all_colourmodels[[srgb_indx]])
  
  names(LAB_tiff)<-c("L","A","B")
  names(srgb_tiff)<-c("red","green","blue")
  
  #add Chroma and Hue to image
  LAB_tiff$Chroma<-sqrt(LAB_tiff$A^2 + LAB_tiff$B^2)
  LAB_tiff$Hue <- (atan2(LAB_tiff$B, LAB_tiff$A) * 180 / pi + 360) %% 360
  
  srgb_tiff2<-srgb_tiff/65535*255
  RGB(srgb_tiff2)<-c(1,2,3)
  hsv_tiff<-colorize(srgb_tiff2,to='hsv')
  names(hsv_tiff)<-c('Hue_sRGB','Saturation_sRGB','Value_sRGB')
  
  combo_tif<-c(LAB_tiff,srgb_tiff,hsv_tiff)
  rm(LAB_tiff,srgb_tiff,hsv_tiff,srgb_tiff2)
  
  # mask the image
  combo_tif <- circ_mask(combo_tif)
  ###############
  
  print(paste0(Code,' starting_acc = ',acc_results[[acidx]]$accuracy,' starting best = ',best[[i]]))
  #####################################################################
  ################# apply the RF models
  #####################################################################
  #if(best[[i]]==1){
    # get RF classifier votes and maximum likelihood for each pixel
    whichmax_img <- terra::predict(combo_tif, classifier_RF_1, na.rm = T,type='response')
    plot(whichmax_img)
    writeRaster(whichmax_img,paste0(end_path,"classified_imgs_17Nov24/",Code,"_individ_seg_classified.tif"),overwrite=T)

    #best[which(names(SDF)%in%redos)]

    samps<-samps_all[j == k,]
    gref<-dummy_georef(whichmax_img)
    rm(whichmax_img)
    samps$predicted<-extract(gref,samps)[,2]
    samps$predicted<-as.character(class_to_factor(samps$predicted))
    sampsdf<-na.omit(as.data.frame(samps))
    colnames(sampsdf)[1]<-'observed'

    acc_results[[acidx]]<-data.frame(accuracy=(sum(sampsdf$predicted==sampsdf$observed)/length(sampsdf$observed)),Code=Code)
    tables[[acidx]]<-table(sampsdf$predicted,sampsdf$observed)
    names(tables)[acidx]<-Code
    sampsdf$Code<-Code
    results[[acidx]]<-sampsdf
    best[[i]]<-1

    if(acc_results[[acidx]]$accuracy<0.85){
      print(paste0(Code,' still sucks ',acc_results[[acidx]]$accuracy))
    }else{
      print(paste0(Code,' success ',acc_results[[acidx]]$accuracy))
    }
  #}

  #if(best[[i]]<=2){
    # train the classifier
    var.imp1 <- data.frame(importance(classifier_RF_1, type=2))
    var.imp1$Variables <- row.names(var.imp1)
    varimp1 <- var.imp1[order(var.imp1$MeanDecreaseGini,decreasing = T),]
    pixel_vars<-which(colnames(train)%in%varimp1$Variables[(varimp1$MeanDecreaseGini/sum(varimp1$MeanDecreaseGin))>.10])
    classifier_RF_2 = randomForest(x = train[,pixel_vars],
                                 y = train$Class,
                                 ntree = 200,mtry=2,importance = T)
    whichmax_img <- terra::predict(combo_tif, classifier_RF_2, na.rm = T,type='response')
    #plot(whichmax_img)
    writeRaster(whichmax_img,paste0(end_path,"classified_imgs_17Nov24/",Code,"_individ_seg_classified.tif"),overwrite=T)

    gref<-dummy_georef(whichmax_img)
    rm(whichmax_img)
    samps<-samps_all[j == k,]
    samps$predicted<-extract(gref,samps)[,2]
    samps$predicted<-as.character(class_to_factor(samps$predicted))
    sampsdf<-na.omit(as.data.frame(samps))
    colnames(sampsdf)[1]<-'observed'

    if((sum(sampsdf$predicted==sampsdf$observed)/length(sampsdf$observed))<acc_results[[acidx]]$accuracy|best[[i]]==2){
      best[[i]]<-1
      print('getting worse')
    }else{
      acc_results[[acidx]]<-data.frame(accuracy=(sum(sampsdf$predicted==sampsdf$observed)/length(sampsdf$observed)),Code=Code)
      tables[[acidx]]<-table(sampsdf$predicted,sampsdf$observed)
      names(tables)[acidx]<-Code
      sampsdf$Code<-Code
      results[[acidx]]<-sampsdf
      best[[i]]<-2
    }
    if(acc_results[[acidx]]$accuracy<0.85){
      print(paste0(Code,' still sucks ',acc_results[[acidx]]$accuracy))
    }else{
      print(paste0(Code,' success ',acc_results[[acidx]]$accuracy))
    }
  #}

  #if(best[[i]]==3){
    # train the classifier
    # read in shp
    idx_d<-which(gsub("_classes","",tools::file_path_sans_ext(basename(shapefiles)))==Code)
    Categorised_shp<-st_read(shapefiles[idx_d])
    Categorised_shp<-Categorised_shp[!is.na(Categorised_shp$Class),]
    Categorised_shp<-Categorised_shp[!is.null(Categorised_shp$Class),]
    Categorised_shp<-Categorised_shp[!Categorised_shp$Class=='',]
    Categorised_shp$Class<-as.factor(Categorised_shp$Class)
    Categorised_shp<-Categorised_shp[Categorised_shp$Class%in%c("Dark wood","Leaves","Light wood","Sky","Sun"),]
    Categorised_shp$Class<-droplevels(Categorised_shp$Class)


    # sample locations within each class
    grid = rast(Categorised_shp, nrow=nrow(combo_tif), ncol=ncol(combo_tif))
    ncr = terra::rasterize(Categorised_shp, grid, field="Class")

    samps_train<-samps_all[j!=k,]
    samps_train$min_class<-terra::extract(ncr,samps_train,ID=F)
    samps_train<-na.omit(as.data.frame(samps_train))
    samps_train$min_class<-as.factor(samps_train$min_class)
    samps_train$min_class<-droplevels(samps_train$min_class)

    #figure out most important variables
    var.imp1 <- data.frame(importance(classifier_RF_1, type=2))
    var.imp1$Variables <- row.names(var.imp1)
    varimp1 <- var.imp1[order(var.imp1$MeanDecreaseGini,decreasing = T),]
    pixel_vars<-which(colnames(samps_train)%in%varimp1$Variables[(varimp1$MeanDecreaseGini/sum(varimp1$MeanDecreaseGin))>.10])

    #retrain the classifier
    classifier_RF_3 = randomForest(x = samps_train[,pixel_vars],
                                 y = samps_train$min_class,
                                 ntree = 200,mtry=2,importance = T)

    # Predict raster levels using the classifier
    whichmax_img <- terra::predict(combo_tif, classifier_RF_3, na.rm = TRUE, type = 'response')
    plot(whichmax_img)

    # Define the broad categories and their numeric mappings
    classes <- data.frame(level = c(0, 1, 2), class = c('Sky', 'Leaves', 'Stems'))

    # Update levels of the raster to their broader categories
    mx_lev <- levels(whichmax_img)[[1]]
    mx_lev$class[mx_lev$class %in% c('Dark wood', 'Light wood')] <- 'Stems'
    mx_lev$class[mx_lev$class %in% c('Sky', 'Sun')] <- 'Sky'
    levels(whichmax_img)[[1]] <- mx_lev

    # Map the factor levels of the raster to their corresponding numeric indices
    factor_values <- values(whichmax_img)  # Extract factor values (as integers)
    factor_classes <- mx_lev$class[factor_values]  # Map to updated class names
    class_map <- setNames(classes$level, classes$class)  # Create mapping
    numeric_values <- class_map[factor_classes]  # Map class names to numeric indices
    numeric_values[!is.na(numeric_values)]


    # Assign numeric values back to a new raster
    values(whichmax_img) <- numeric_values
    writeRaster(whichmax_img,paste0(end_path,"classified_imgs_17Nov24/",Code,"_individ_seg_classified.tif"),overwrite=T)


    gref<-dummy_georef(whichmax_img)
    rm(whichmax_img)
    samps<-samps_all[j == k,]
    samps$predicted<-as.character(class_to_factor(extract(gref,samps)[,2]))

    # Convert the numeric observations to class names
    sampsdf<-na.omit(as.data.frame(samps))
    colnames(sampsdf)[1]<-'observed'

     if((sum(sampsdf$predicted==sampsdf$observed)/length(sampsdf$observed))<acc_results[[acidx]]$accuracy|best[[i]]==3){
       print('getting worse')
     }else{
      acc_results[[acidx]]<-data.frame(accuracy=(sum(sampsdf$predicted==sampsdf$observed)/length(sampsdf$observed)),Code=Code)
      tables[[acidx]]<-table(sampsdf$predicted,sampsdf$observed)
      names(tables)[acidx]<-Code
      sampsdf$Code<-Code
      results[[acidx]]<-sampsdf
      best[[i]]<-3
     }
    if(acc_results[[acidx]]$accuracy<0.85){
      print(paste0(Code,' still sucks ',acc_results[[acidx]]$accuracy))
    }else{
      print(paste0(Code,' success ',acc_results[[acidx]]$accuracy))
    }
  #}

  
  #if(best[[i]]==4){
    # train the classifier
    # read in shp
    idx_d<-which(gsub("_classes","",tools::file_path_sans_ext(basename(shapefiles)))==Code)
    Categorised_shp<-st_read(shapefiles[idx_d])
    Categorised_shp<-Categorised_shp[!is.na(Categorised_shp$Class),]
    Categorised_shp<-Categorised_shp[!is.null(Categorised_shp$Class),]
    Categorised_shp<-Categorised_shp[!Categorised_shp$Class=='',]
    Categorised_shp$Class<-as.factor(Categorised_shp$Class)
    Categorised_shp<-Categorised_shp[Categorised_shp$Class%in%c("Dark wood","Leaves","Light wood","Sky","Sun"),]
    Categorised_shp$Class<-droplevels(Categorised_shp$Class)
    
    
    # sample locations within each class
    grid = rast(Categorised_shp, nrow=nrow(combo_tif), ncol=ncol(combo_tif))
    ncr = terra::rasterize(Categorised_shp, grid, field="Class")
    samps_train<-samps_all[j!=k,]
    samps_train$min_class<-terra::extract(ncr,samps_train,ID=F)
    samps_train<-na.omit(as.data.frame(samps_train))
    samps_train$min_class<-as.factor(samps_train$min_class)
    samps_train$min_class<-droplevels(samps_train$min_class)
    pixel_vars<-which(colnames(samps_train)%in%c('L','A','B','Hue','Chroma',"red","green","blue", "Hue_sRGB","Saturation_sRGB","Value_sRGB"))

    classifier_RF_4 = randomForest(x = samps_train[,pixel_vars], 
                                 y = samps_train$min_class,
                                 ntree = 200,mtry=2,importance = T)
    
    # Predict raster levels using the classifier
    whichmax_img <- terra::predict(combo_tif, classifier_RF_4, na.rm = TRUE, type = 'response')
    
    # Define the broad categories and their numeric mappings
    classes <- data.frame(level = c(0, 1, 2), class = c('Sky', 'Leaves', 'Stems'))
    
    # Update levels of the raster to their broader categories
    mx_lev <- levels(whichmax_img)[[1]]
    mx_lev$class[mx_lev$class %in% c('Dark wood', 'Light wood')] <- 'Stems'
    mx_lev$class[mx_lev$class %in% c('Sky', 'Sun')] <- 'Sky'
    levels(whichmax_img)[[1]] <- mx_lev
    
    # Map the factor levels of the raster to their corresponding numeric indices
    factor_values <- values(whichmax_img)  # Extract factor values (as integers)
    factor_classes <- mx_lev$class[factor_values]  # Map to updated class names
    class_map <- setNames(classes$level, classes$class)  # Create mapping
    numeric_values <- class_map[factor_classes]  # Map class names to numeric indices
    numeric_values[!is.na(numeric_values)]
    
    
    # Assign numeric values back to a new raster
    values(whichmax_img) <- numeric_values
    writeRaster(whichmax_img,paste0(end_path,"classified_imgs_17Nov24/",Code,"_individ_seg_classified.tif"),overwrite=T)
    
    
    gref<-dummy_georef(whichmax_img)
    rm(whichmax_img)
    plot(gref)
    plot(samps_all[j!=k,],add=T)
    
    samps<-samps_all[j == k,]
    samps$predicted<-as.character(class_to_factor(extract(gref,samps)[,2]))
    
    
    # Convert the numeric observations to class names
    sampsdf<-as.data.frame(samps)
    sampsdf<-na.omit(sampsdf)
    #idxofNA<-which(is.na(sampsdf),arr.ind = T)
    #sampsdf[idxofNA[,1],]
    colnames(sampsdf)[1]<-'observed'

    if((sum(sampsdf$predicted==sampsdf$observed)/length(sampsdf$observed))<acc_results[[acidx]]$accuracy|best[[i]]==4){
      print('getting worse')
    }else{
      acc_results[[acidx]]<-data.frame(accuracy=(sum(sampsdf$predicted==sampsdf$observed)/length(sampsdf$observed)),Code=Code)
      tables[[acidx]]<-table(sampsdf$predicted,sampsdf$observed)
      names(tables)[acidx]<-Code
      sampsdf$Code<-Code
      results[[acidx]]<-sampsdf
      best[[i]]<-4
    }
    if(acc_results[[acidx]]$accuracy<0.85){
      print(paste0(Code,' still sucks ',acc_results[[acidx]]$accuracy))
    }else{
      print(paste0(Code,' success ',acc_results[[acidx]]$accuracy))
    }
  #}
  
  # 
  #if(best[[i]]<=5){
    # train the classifier
    # read in shp
    idx_d<-which(gsub("_classes","",tools::file_path_sans_ext(basename(shapefiles)))==Code)
    Categorised_shp<-st_read(shapefiles[idx_d])
    Categorised_shp<-Categorised_shp[!is.na(Categorised_shp$Class),]
    Categorised_shp<-Categorised_shp[!is.null(Categorised_shp$Class),]
    Categorised_shp<-Categorised_shp[!Categorised_shp$Class=='',]
    Categorised_shp$Class<-as.factor(Categorised_shp$Class)
    Categorised_shp<-Categorised_shp[Categorised_shp$Class%in%c("Dark wood","Leaves","Light wood","Sky","Sun"),]
    Categorised_shp$Class<-droplevels(Categorised_shp$Class)


    # sample locations within each class
    grid = rast(Categorised_shp, nrow=nrow(combo_tif), ncol=ncol(combo_tif))
    ncr = terra::rasterize(Categorised_shp, grid, field="Class")
    samps_train<-samps_all[j!=k,]
    samps_train$min_class<-terra::extract(ncr,samps_train,ID=F)
    samps_train<-na.omit(as.data.frame(samps_train))
    samps_train$min_class<-as.factor(samps_train$min_class)
    samps_train$min_class<-droplevels(samps_train$min_class)

    pixel_vars<-which(colnames(samps_train)%in%c('L','A','B','Hue','Chroma',"red","green","blue", "Hue_sRGB","Saturation_sRGB","Value_sRGB"))


    sky<-samps_train[samps_train$min_class=='Sky',]
    nosky<-samps_train[!samps_train$min_class=='Sky',]
    sky<-sky[order(sky$blue),]
    sky$min_class<-cut(sky$blue,breaks=3,labels=c('sky1','sky2','sky3'))
    samps_train<-rbind(sky,nosky)
    samps_train$min_class<-droplevels(samps_train$min_class)

    classifier_RF_5 = randomForest(x = samps_train[,pixel_vars],
                                 y = samps_train$min_class,
                                 ntree = 250,mtry=2,importance = T)

    # Predict raster levels using the classifier
    whichmax_img <- terra::predict(combo_tif, classifier_RF_5, na.rm = TRUE, type = 'response')

    # Define the broad categories and their numeric mappings
    classes <- data.frame(level = c(0, 1, 2), class = c('Sky', 'Leaves', 'Stems'))

    # Update levels of the raster to their broader categories
    mx_lev <- levels(whichmax_img)[[1]]
    mx_lev$class[mx_lev$class %in% c('Dark wood', 'Light wood')] <- 'Stems'
    mx_lev$class[mx_lev$class %in% c('Sky', 'Sun','sky1','sky2','sky3','sky4')] <- 'Sky'
    levels(whichmax_img)[[1]] <- mx_lev

    # Map the factor levels of the raster to their corresponding numeric indices
    factor_values <- values(whichmax_img)  # Extract factor values (as integers)
    factor_classes <- mx_lev$class[factor_values]  # Map to updated class names
    class_map <- setNames(classes$level, classes$class)  # Create mapping
    numeric_values <- class_map[factor_classes]  # Map class names to numeric indices
    numeric_values[!is.na(numeric_values)]


    # Assign numeric values back to a new raster
    values(whichmax_img) <- numeric_values
    plot(whichmax_img)
    writeRaster(whichmax_img,paste0(end_path,"classified_imgs_17Nov24/",Code,"_individ_seg_classified.tif"),overwrite=T)


    gref<-dummy_georef(whichmax_img)
    rm(whichmax_img)
    samps<-samps_all[j == k,]
    samps$predicted<-as.character(class_to_factor(extract(gref,samps)[,2]))

    # Convert the numeric observations to class names
    sampsdf<-na.omit(as.data.frame(samps))
    colnames(sampsdf)[1]<-'observed'

    if((sum(sampsdf$predicted==sampsdf$observed)/length(sampsdf$observed))<acc_results[[acidx]]$accuracy|best[[i]]==5){
      print('getting worse')
    }else{
      acc_results[[acidx]]<-data.frame(accuracy=(sum(sampsdf$predicted==sampsdf$observed)/length(sampsdf$observed)),Code=Code)
      tables[[acidx]]<-table(sampsdf$predicted,sampsdf$observed)
      names(tables)[acidx]<-Code
      sampsdf$Code<-Code
      results[[acidx]]<-sampsdf
      best[[i]]<-5
    }
    if(acc_results[[acidx]]$accuracy<0.85){
      print(paste0(Code,' still sucks ',acc_results[[acidx]]$accuracy))
    }else{
      print(paste0(Code,' success ',acc_results[[acidx]]$accuracy))
    }
  #}
    
    whichmax_img <- rast(paste0(end_path,"classified_imgs_17Nov24/",Code,"_individ_seg_classified.tif"))
    


  ############################################
  ## part 2 - run the best classifier
  ############################################
  if(best[[i]]==1){
    whichmax_img <- terra::predict(combo_tif, classifier_RF_1, na.rm = T,type='response')
    writeRaster(whichmax_img,paste0(end_path,"classified_imgs_17Nov24/",Code,"_individ_seg_classified.tif"),overwrite=T)
  }

  if(best[[i]]==2){
    whichmax_img <- terra::predict(combo_tif, classifier_RF_2, na.rm = T,type='response')
    writeRaster(whichmax_img,paste0(end_path,"classified_imgs_17Nov24/",Code,"_individ_seg_classified.tif"),overwrite=T)
  }

  if(best[[i]]==3){
    # Predict raster levels using the classifier
    whichmax_img <- terra::predict(combo_tif, classifier_RF_3, na.rm = TRUE, type = 'response')

    # Define the broad categories and their numeric mappings
    classes <- data.frame(level = c(0, 1, 2), class = c('Sky', 'Leaves', 'Stems'))

    # Update levels of the raster to their broader categories
    mx_lev <- levels(whichmax_img)[[1]]
    mx_lev$class[mx_lev$class %in% c('Dark wood', 'Light wood')] <- 'Stems'
    mx_lev$class[mx_lev$class %in% c('Sky', 'Sun')] <- 'Sky'
    levels(whichmax_img)[[1]] <- mx_lev

    # Map the factor levels of the raster to their corresponding numeric indices
    factor_values <- values(whichmax_img)  # Extract factor values (as integers)
    factor_classes <- mx_lev$class[factor_values]  # Map to updated class names
    class_map <- setNames(classes$level, classes$class)  # Create mapping
    numeric_values <- class_map[factor_classes]  # Map class names to numeric indices
    numeric_values[!is.na(numeric_values)]

    # Assign numeric values back to a new raster
    values(whichmax_img) <- numeric_values
    writeRaster(whichmax_img,paste0(end_path,"classified_imgs_17Nov24/",Code,"_individ_seg_classified.tif"),overwrite=T)
  }

  if(best[[i]]==4){
    # Predict raster levels using the classifier
    whichmax_img <- terra::predict(combo_tif, classifier_RF_4, na.rm = TRUE, type = 'response')

    # Define the broad categories and their numeric mappings
    classes <- data.frame(level = c(0, 1, 2), class = c('Sky', 'Leaves', 'Stems'))

    # Update levels of the raster to their broader categories
    mx_lev <- levels(whichmax_img)[[1]]
    mx_lev$class[mx_lev$class %in% c('Dark wood', 'Light wood')] <- 'Stems'
    mx_lev$class[mx_lev$class %in% c('Sky', 'Sun')] <- 'Sky'
    levels(whichmax_img)[[1]] <- mx_lev

    # Map the factor levels of the raster to their corresponding numeric indices
    factor_values <- values(whichmax_img)  # Extract factor values (as integers)
    factor_classes <- mx_lev$class[factor_values]  # Map to updated class names
    class_map <- setNames(classes$level, classes$class)  # Create mapping
    numeric_values <- class_map[factor_classes]  # Map class names to numeric indices
    numeric_values[!is.na(numeric_values)]

    # Assign numeric values back to a new raster
    values(whichmax_img) <- numeric_values
    writeRaster(whichmax_img,paste0(end_path,"classified_imgs_17Nov24/",Code,"_individ_seg_classified.tif"),overwrite=T)
  }

  #############################################
  ## read in the optimised image
  #############################################
  writeRaster(dummy_georef(whichmax_img),paste0(end_path,"georef_classified_imgs_17Nov24/",Code,"_individ_seg_classified.tif"),overwrite=T)
  rm(whichmax_img)
  whichmax_img <- rast(paste0(end_path,"classified_imgs_17Nov24/",Code,"_individ_seg_classified.tif"))
  
  ###########################################
  if(acc_results[[acidx]]$accuracy>0.85){
    print('starting seg')
    # segment image at 2 degree resolution
    Seg_img<-segment(whichmax_img)
      # quantify leaf to wood ratio in each segment
    Seg_img$azimuth_group<-as.factor(Seg_img$azimuth_group)
    Seg_img$ring<-as.factor(Seg_img$ring)
    Seg_img$class<-as.factor(Seg_img$class)
    rr<-list()
    for (r in 1:length(levels(Seg_img$ring))){
      tt<-list()
      for (t in 1:length(levels(Seg_img$azimuth_group))){
        ring<-levels(Seg_img$ring)[r]
        azimuth_group<-levels(Seg_img$azimuth_group)[t]
        class_prop<-as.data.frame(t(as.matrix(table(Seg_img[Seg_img$ring==ring&Seg_img$azimuth_group==azimuth_group,]$class))))
        colnames(class_prop)<-as.character(class_to_factor(colnames(class_prop)))
        a<-class_prop$Stems/(class_prop$Stems+class_prop$Leaves)
        Pgap<-class_prop$Sky/(class_prop$Sky+class_prop$Leaves)
        tt[[t]]<-cbind(ring,azimuth_group,class_prop,a,Pgap,Code)
      }
      rr[[r]]<-do.call(rbind,tt)
    }
    idx<-which(names(class_by_seg)==Code)
    if(length(idx)<1){idx<-length(class_by_seg)+1}
    class_by_seg[[idx]]<-do.call(rbind,rr)
    class_by_seg[[idx]]$ring<-as.factor(class_by_seg[[idx]]$ring)
    ringtot<-aggregate(class_by_seg[[idx]][,3:5],by=list(class_by_seg[[idx]]$ring),FUN=sum)
    total_a<-mean(ringtot$Stems/(ringtot$Stems+ringtot$Leaves),na.rm=T)
    sd_a<-sd(ringtot$Stems/(ringtot$Stems+ringtot$Leaves),na.rm=T)
    
    idx<-which(names(img_a)==Code)
    if(length(idx)<1){idx<-length(img_a)+1}
    img_a[[idx]]<-cbind(total_a,sd_a,Code)
    ring_a<-ringtot$Stems/(ringtot$Stems+ringtot$Leaves)
    
    idx<-which(names(class_by_ring)==Code)
    if(length(idx)<1){idx<-length(class_by_ring)+1}
    class_by_ring[[idx]]<-cbind(ringtot,ring_a,Code)
    class_by_ring[[idx]]$Code<-Code
  
    # segment image for LAI calculation  @ 15 degree azm
    seg_img_adj<-adj_seg(Seg_img,
                         nseg = (360/15),  # Number of azimuth segments
                         startVZA=0, # start viewing zenith angle
                         endVZA=70, # end viewing zenith angle
                         maxVZA=90, # maximum view angle of the camera
                         nrings=7) # for the number of zenith segments)
    seg_img_adj$azimuth_group<-as.factor(seg_img_adj$azimuth_group)
    seg_img_adj$ring<-as.factor(seg_img_adj$ring)
    rr<-list()
    for (r in 1:length(levels(seg_img_adj$ring))){
      tt<-list()
      for (t in 1:length(levels(seg_img_adj$azimuth_group))){
        ring<-levels(seg_img_adj$ring)[r]
        azimuth_group<-levels(seg_img_adj$azimuth_group)[t]
        class_prop<-as.data.frame(t(as.matrix(table(seg_img_adj[seg_img_adj$ring==ring&seg_img_adj$azimuth_group==azimuth_group,]$class))))
        colnames(class_prop)<-as.character(class_to_factor(colnames(class_prop)))
        a<-class_prop$Stems/(class_prop$Stems+class_prop$Leaves)
        Pgap<-class_prop$Sky/(class_prop$Sky+class_prop$Leaves)
        tt[[t]]<-cbind(ring,azimuth_group,class_prop,a,Pgap,Code)
      }
      rr[[r]]<-do.call(rbind,tt)
    }
    
    idx<-which(names(class_by_seg_4LAI)==Code)
    if(length(idx)<1){idx<-length(class_by_seg_4LAI)+1}
    class_by_seg_4LAI[[idx]]<-do.call(rbind,rr)
    class_by_seg_4LAI[[idx]]$ring<-as.factor(class_by_seg_4LAI[[idx]]$ring)
    
    idx<-which(names(img_a_4LAI)==Code)
    if(length(idx)<1){idx<-length(img_a_4LAI)+1}
    total_a<-mean(ringtot$Stems/(ringtot$Stems+ringtot$Leaves),na.rm=T)
    sd_a<-sd(ringtot$Stems/(ringtot$Stems+ringtot$Leaves),na.rm=T)
    img_a_4LAI[[idx]]<-cbind(total_a,sd_a,Code)
    
    
    # segment image for LAI calculation  @ 15 degree azm
    # seg_img_adj_60<-adj_seg(Seg_img,
    #                      nseg = (360/15),  # Number of azimuth segments
    #                      startVZA=0, # start viewing zenith angle
    #                      endVZA=60, # end viewing zenith angle
    #                      maxVZA=90, # maximum view angle of the camera
    #                      nrings=6,)# for the number of zenith segments
    # seg_img_adj_60$azimuth_group<-as.factor(seg_img_adj_60$azimuth_group)
    # seg_img_adj_60$ring<-as.factor(seg_img_adj_60$ring)
    # rr<-list()
    # for (r in 1:length(levels(seg_img_adj_60$ring))){
    #   tt<-list()
    #   for (t in 1:length(levels(seg_img_adj_60$azimuth_group))){
    #     ring<-levels(seg_img_adj_60$ring)[r]
    #     azimuth_group<-levels(seg_img_adj_60$azimuth_group)[t]
    #     class_prop<-as.data.frame(t(as.matrix(table(seg_img_adj_60[seg_img_adj_60$ring==ring&seg_img_adj_60$azimuth_group==azimuth_group,]$class))))
    #     colnames(class_prop)<-as.character(class_to_factor(colnames(class_prop)))
    #     a<-class_prop$Stems/(class_prop$Stems+class_prop$Leaves)
    #     Pgap<-class_prop$Sky/(class_prop$Sky+class_prop$Leaves)
    #     tt[[t]]<-cbind(ring,azimuth_group,class_prop,a,Pgap,Code)
    #   }
    #   rr[[r]]<-do.call(rbind,tt)
    # }
    # class_by_seg_4LAI60deg[[i]]<-do.call(rbind,rr)
    # class_by_seg_4LAI60deg[[i]]$ring<-as.factor(class_by_seg_4LAI60deg[[i]]$ring)
    # 
    
    # segment image for LAI2200 calculation
    seg_img_adj<-adj_seg(Seg_img,
                         nseg = (360/45),  # Number of azimuth segments
                         startVZA=0, # start viewing zenith angle
                         endVZA=75, # end viewing zenith angle
                         maxVZA=90, # maximum view angle of the camera
                         nrings=5)# for the number of zenith segments
    seg_img_adj$azimuth_group<-as.factor(seg_img_adj$azimuth_group)
    seg_img_adj$ring<-as.factor(seg_img_adj$ring)
    rr<-list()
    for (r in 1:length(levels(seg_img_adj$ring))){
      tt<-list()
      for (t in 1:length(levels(seg_img_adj$azimuth_group))){
        ring<-levels(seg_img_adj$ring)[r]
        azimuth_group<-levels(seg_img_adj$azimuth_group)[t]
        class_prop<-as.data.frame(t(as.matrix(table(seg_img_adj[seg_img_adj$ring==ring&seg_img_adj$azimuth_group==azimuth_group,]$class))))
        colnames(class_prop)<-as.character(class_to_factor(colnames(class_prop)))
        a<-class_prop$Stems/(class_prop$Stems+class_prop$Leaves)
        Pgap<-class_prop$Sky/(class_prop$Sky+class_prop$Leaves)
        tt[[t]]<-cbind(ring,azimuth_group,class_prop,a,Pgap,Code)
      }
      rr[[r]]<-do.call(rbind,tt)
    }
    idx<-which(names(class_by_seg_4LAI2200)==Code)
    if(length(idx)<1){idx<-length(class_by_seg_4LAI2200)+1}
    class_by_seg_4LAI2200[[idx]]<-do.call(rbind,rr)
    class_by_seg_4LAI2200[[idx]]$ring<-as.factor(class_by_seg_4LAI2200[[idx]]$ring)
    
    # segment image for 57deg calculation
    seg_img_adj<-adj_seg(Seg_img,
                         nseg = (360/45),  # Number of azimuth segments
                         startVZA=57.3-2.5, # start viewing zenith angle
                         endVZA=57.3+2.5, # end viewing zenith angle
                         maxVZA=90, # maximum view angle of the camera
                         nrings=1)# for the number of zenith segments
    seg_img_adj$azimuth_group<-as.factor(seg_img_adj$azimuth_group)
    seg_img_adj$ring<-as.factor(seg_img_adj$ring)
    rr<-list()
    for (r in 1:length(levels(seg_img_adj$ring))){
      tt<-list()
      for (t in 1:length(levels(seg_img_adj$azimuth_group))){
        ring<-levels(seg_img_adj$ring)[r]
        azimuth_group<-levels(seg_img_adj$azimuth_group)[t]
        class_prop<-as.data.frame(t(as.matrix(table(seg_img_adj[seg_img_adj$ring==ring&seg_img_adj$azimuth_group==azimuth_group,]$class))))
        colnames(class_prop)<-as.character(class_to_factor(colnames(class_prop)))
        a<-class_prop$Stems/(class_prop$Stems+class_prop$Leaves)
        Pgap<-class_prop$Sky/(class_prop$Sky+class_prop$Leaves)
        tt[[t]]<-cbind(ring,azimuth_group,class_prop,a,Pgap,Code)
      }
      rr[[r]]<-do.call(rbind,tt)
    }
    idx<-which(names(class_by_seg_57deg)==Code)
    if(length(idx)<1){idx<-length(class_by_seg_57deg)+1}
    class_by_seg_57deg[[idx]]<-do.call(rbind,rr)
    class_by_seg_57deg[[idx]]$ring<-as.factor(class_by_seg_57deg[[idx]]$ring)
    rm(Seg_img)
    rm(seg_img_adj)
  }

  print(i)
  end<-Sys.time()
  print(end-start)
}


#######
### save outputs
saveRDS(class_by_ring,paste0(end_path,'segmented_byzenith__45zen_noazm_maxlikelihood_v6.rds'))
saveRDS(img_a,paste0(end_path,'Total_a_by_img_zen0to90_v6.rds'))
saveRDS(img_a_4LAI,paste0(end_path,'Total_a_by_img_zen0to70_v6.rds'))
saveRDS(class_by_seg_4LAI,paste0(end_path,'segmented_7zen_24azm_maxlikelihood_v6.rds'))
saveRDS(class_by_seg_4LAI2200,paste0(end_path,'segmented_5zen_8azm_maxlikelihood_v6.rds'))
saveRDS(class_by_seg_57deg,paste0(end_path,'segmented_zen57deg_8azm_maxlikelihood_v6.rds'))
saveRDS(acc_results,paste0(end_path,'individual_RF_global_accuracy_results_v6.rds'))
saveRDS(tables,paste0(end_path,'missclass_tables_v6.rds'))
saveRDS(results,paste0(end_path,'individual_RF_sample_prediction_accuracy_results_v6.rds'))
saveRDS(var_importance,paste0(end_path,'colour_model_band_importance_v6.rds'))
saveRDS(best,paste0(end_path,'best_model_v6.rds'))
saveRDS(class_by_seg,paste0(end_path,'segmented_45zen_8azm_maxlikelihood_v6.rds'))







# 
# var_importance<-do.call(rbind,var_importance)
# var_importance<-merge(var_importance,img_details,by='Code',all.x=T)
# saveRDS(var_importance,paste0(end_path,'colour_model_band_importance_v5.rds'))
# #indexi<-which(lengths(class_by_seg)>0)
# class_by_seg_df<-do.call(rbind,class_by_seg)
# saveRDS(class_by_seg_df,paste0(end_path,'segmented_45zen_8azm_maxlikelihood_v5.rds'))
# indexi<-which(lengths(class_by_ring)>0)
# class_by_ring<-do.call(rbind,class_by_ring)
# saveRDS(class_by_ring,paste0(end_path,'segmented_byzenith__45zen_noazm_maxlikelihood_v5.rds'))
# saveRDS(tables,paste0(end_path,'missclass_tables_v3.rds'))
# acc_results<-do.call(rbind,acc_results)
# saveRDS(acc_results,paste0(end_path,'individual_RF_global_accuracy_results_v5.rds'))
# results<-do.call(rbind,results)
# saveRDS(results,paste0(end_path,'individual_RF_global_accuracy_results_v5.rds'))
# img_a<-do.call(rbind,img_a)
# saveRDS(img_a,paste0(end_path,'Total_a_by_img_zen0to90_v5.rds'))
# img_a_4LAI<-do.call(rbind,img_a_4LAI)
# saveRDS(img_a_4LAI,paste0(end_path,'Total_a_by_img_zen0to70_v5.rds'))
# class_by_seg_4LAI<-do.call(rbind,class_by_seg_4LAI)
# saveRDS(class_by_seg_4LAI,paste0(end_path,'segmented_7zen_24azm_maxlikelihood_v5.rds'))
# class_by_seg_4LAI2200<-do.call(rbind,class_by_seg_4LAI2200)
# saveRDS(class_by_seg_4LAI2200,paste0(end_path,'segmented_5zen_8azm_maxlikelihood_v5.rds'))
# # class_by_seg_4LAI60deg<-do.call(rbind,class_by_seg_4LAI60deg)
# # saveRDS(class_by_seg_4LAI60deg,paste0(end_path,'segmented_6zen_24azm_maxlikelihood_v5.rds'))
# class_by_seg_57deg<-do.call(rbind,class_by_seg_57deg)
# saveRDS(class_by_seg_57deg,paste0(end_path,'segmented_zen57deg_8azm_maxlikelihood_v5.rds'))
# saveRDS(best,paste0(end_path,'best_model_v5.rds'))

###################################################################################################

# class_by_ring_old<-readRDS(paste0(end_path,'segmented_byzenith__45zen_noazm_maxlikelihood_v3.rds'))
# class_by_ring_old<-class_by_ring_old[!which(class_by_ring_old$Code%in%unique(class_by_ring$Code)),]
# class_by_ring<-rbind(class_by_ring,class_by_ring_old)
# saveRDS(class_by_ring,paste0(end_path,'segmented_byzenith__45zen_noazm_maxlikelihood_v4.rds'))
# 
# img_a_old<-readRDS(paste0(end_path,'Total_a_by_img_zen0to90_v3.rds'))
# img_a_old<-img_a_old[!which(img_a_old$Code%in%unique(img_a$Code)),]
# img_a<-rbind(img_a,img_a_old)
# saveRDS(img_a,paste0(end_path,'Total_a_by_img_zen0to90_v4.rds'))
# 
# img_a_4LAI_old<-readRDS(paste0(end_path,'Total_a_by_img_zen0to70_v3.rds'))
# img_a_4LAI_old<-img_a_4LAI_old[!which(img_a_4LAI_old$Code%in%unique(img_a_4LAI$Code)),]
# img_a_4LAI<-rbind(img_a_4LAI,img_a_4LAI_old)
# saveRDS(img_a_4LAI,paste0(end_path,'Total_a_by_img_zen0to70_v4.rds'))
# 
# class_by_seg_4LAI_old<-readRDS(paste0(end_path,'segmented_7zen_24azm_maxlikelihood_v3.rds'))
# class_by_seg_4LAI_old<-class_by_seg_4LAI_old[!which(class_by_seg_4LAI_old$Code%in%unique(class_by_seg_4LAI$Code)),]
# class_by_seg_4LAI<-rbind(class_by_seg_4LAI,class_by_seg_4LAI_old)
# saveRDS(class_by_seg_4LAI,paste0(end_path,'segmented_7zen_24azm_maxlikelihood_v4.rds'))
# 
# class_by_seg_4LAI2200_old<-readRDS(paste0(end_path,'segmented_5zen_8azm_maxlikelihood_v3.rds'))
# class_by_seg_4LAI2200_old<-class_by_seg_4LAI2200_old[!which(class_by_seg_4LAI2200_old$Code%in%unique(class_by_seg_4LAI2200$Code)),]
# class_by_seg_4LAI2200<-rbind(class_by_seg_4LAI2200,class_by_seg_4LAI2200_old)
# saveRDS(class_by_seg_4LAI2200,paste0(end_path,'segmented_7zen_24azm_maxlikelihood_v4.rds'))
# 
# class_by_seg_57deg_old<-readRDS(paste0(end_path,'segmented_zen57deg_8azm_maxlikelihood_v3.rds'))
# class_by_seg_57deg_old<-class_by_seg_57deg_old[!which(class_by_seg_57deg_old$Code%in%unique(class_by_seg_57deg$Code)),]
# class_by_seg_57deg<-rbind(class_by_seg_57deg,class_by_seg_57deg_old)
# saveRDS(class_by_seg_57deg,paste0(end_path,'segmented_zen57deg_8azm_maxlikelihood_v4.rds'))
# 
# acc_results_old<-readRDS(paste0(end_path,'individual_RF_global_accuracy_results__v3.rds'))
# acc_results_old<-acc_results_old[!which(acc_results_old$Code%in%unique(acc_results$Code)),]
# acc_results<-rbind(acc_results,acc_results_old)
# saveRDS(acc_results,paste0(end_path,'individual_RF_global_accuracy_results__v4.rds'))
# 
# tables_old<-readRDS(paste0(end_path,'missclass_tables_v3.rds'))
# tables_old[which(lapply(tables,function(x)is.null(x)))]<-tables[which(lapply(tables,function(x)is.null(x)))]
# saveRDS(acc_results,paste0(end_path,'missclass_tables_v4.rds'))
# 
# results_old<-readRDS(paste0(end_path,'individual_RF_sample_prediction_accuracy_results_v3.rds'))
# results_old<-results_old[!which(results_old$Code%in%unique(results$Code)),]
# results<-rbind(results,results_old)
# saveRDS(results,paste0(end_path,'individual_RF_sample_prediction_accuracy_results_v4.rds'))
# 
# var_importance_old<-readRDS(paste0(end_path,'colour_model_band_importance_v3.rds'))
# var_importance_old<-var_importance_old[!which(var_importance_old$Code%in%unique(var_importance$Code)),]
# var_importance<-rbind(var_importance,var_importance_old)
# saveRDS(results,paste0(end_path,'colour_model_band_importance_v4.rds'))
# 
# class_by_ring_old<-readRDS(paste0(end_path,'segmented_byzenith__45zen_noazm_maxlikelihood_v3.rds'))
# class_by_ring_old<-class_by_ring_old[!which(class_by_ring_old$Code%in%unique(class_by_ring$Code)),]
# class_by_ring<-rbind(class_by_ring,class_by_ring_old)
# saveRDS(class_by_ring,paste0(end_path,'segmented_byzenith__45zen_noazm_maxlikelihood_v4.rds'))
# 













#test which ones were improved and overall improvement
acc_results_orig<-read.csv(paste0(end_path,'/individual_RF_global_accuracy_results.csv'))
acc_results_orig[which(acc_results_orig$accuracy<.9),]
acc_results_2<-do.call(rbind,readRDS(paste0(end_path,'individual_RF_global_accuracy_results_v6.rds')))
nrow(acc_results_2[acc_results_2$accuracy>0.85,])
nrow(acc_results_orig[acc_results_orig$accuracy>0.9,])

for (i in 1:nrow(acc_results_orig)){
  indexio<-which(acc_results$Code==acc_results_orig$Code[i])
  if(length(indexio)>0){
    if(acc_results_orig$accuracy[i]<.9&signif(acc_results$accuracy,2)[indexio]>=0.9){
      print(acc_results$Code[indexio])
    }
  }
}
acc_results_orig<-merge(acc_results_orig,img_details,'Code',all.x=T)
acc_results_2<-merge(acc_results_2,img_details,'Code',all.x=T)

t.test(acc_results_orig$accuracy,acc_results$accuracy)

t.test(acc_results_orig$accuracy[acc_results_orig$Sky=='clear'],acc_results_2$accuracy[acc_results_2$Sky=='clear'])
t.test(acc_results_orig$accuracy[acc_results_orig$Sky=='overcast'],acc_results_2$accuracy[acc_results_2$Sky=='overcast'])
t.test(acc_results_orig$accuracy[acc_results_orig$Sky=='partly cloudy'],acc_results_2$accuracy[acc_results_2$Sky=='partly cloudy'])
