#generalised classifier.
### Libraries + functions #####################################################################
library(randomForest)
library(sf)
library(terra)
source('LAI to PAI functions.r')
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
####################################################################################


### file paths #####################################################################

end_path<-"F:/Training tiffs/"
sample_dataframe<-readRDS(paste0(end_path,'Pixel_value_training_dataset_v3.rds'))
img_details<-read.csv(paste0(end_path,"img_details.csv"))
metadata<-read.csv(paste0(end_path,'img_ID_metadata.csv'))[,c(1,9,12)]
sample_dataframe<-merge(sample_dataframe,metadata,by='Code',all.x=T,all.y = F)
all_colourmodels<-c(Sys.glob(file.path("F:/PhotoproRGB_LAB_imgs",'*.tif')),Sys.glob(file.path("F:/PhotoproRGB_LAB_imgs",'*.tiff')))
updated_shapes="F:/Training tiffs/edited_shapes/"
shapefiles<-Sys.glob(file.path(updated_shapes, "*.shp"))
samp_list<-Sys.glob(file.path(paste0(end_path,"Sample locations"),'*.shp'))

#sample_dataframe<-merge(sample_dataframe,img_details[,c(1,5)])

# for (i in 1:nrow(sample_dataframe)){
#   if(sample_dataframe$Sky[i]==''){
#     sample_dataframe$Sky[i]<-img_details$Sky[which(img_details$Code==sample_dataframe$Code[i])]
#   }
# }



### editing sample_dataframe  #####################################################################
sample_dataframe$Class<-as.factor(sample_dataframe$Class)
sample_dataframe$Code<-as.factor(sample_dataframe$Code)
SDF <- split(sample_dataframe, sample_dataframe$Code)
rm(sample_dataframe)

# SDF_metadata<-list()
# for (i in 1:length(SDF)){
#   idx<-which(colnames(SDF[[i]])%in%c("Code","Class","L","A","B","Chroma","Hue"))
#   SDF_metadata[[i]]<- SDF[[i]][1,-idx]
#   SDF_metadata[[i]]<-cbind(data.frame(Code=names(SDF)[i]),SDF_metadata[[i]])
#   SDF[[i]]<- SDF[[i]][,idx]
# }
# SDF_metadata<-do.call(rbind,SDF_metadata)

# for(i in 1:length(SDF)){
#   # check there is a leaf sky and wood category in each image
#   if(!length(levels(SDF[[i]]$Class))==3){
#     print(paste0('Level missing in ',SDF[[i]][1,1]))
#     print(levels(SDF[[i]]$Class))
#     SDF<-SDF[[-i]]
#     #total_a[[i]]<-1
#     next
#   }
#   SDF[[i]]$Class<-droplevels(SDF[[i]]$Class)
# }


###################################################################################################
### General image RF classifier ################################################################
###################################################################################################

####################################################################
######    running test and train procedure for the RF model   ######
####################################################################

###### create groupings for test/train
set.seed(42)
# number of folds
k <- 2
# randomly assign test or train
j<-sample(rep(1:k,each=round(length(SDF))/k))
#splitting test/train dataset
train <- do.call(rbind,SDF[which(j!= k)])
test <- SDF[which(j== k)]
# specify variables to use in the model
pixel_vars<-which(colnames(train)%in%c('L','A','B','Hue','Chroma','red','green','blue','Hue_sRGB','Saturation_sRGB','Value_sRGB'))

# train the classifier
classifier_RF = randomForest(x = train[,pixel_vars], 
                           y = train$Class, 
                           ntree = 200,mtree=2)
## add something that captures 'importance'
gen_var_importance<-as.data.frame(importance(classifier_RF)[order(importance(classifier_RF)[,1],decreasing=T),drop=F,1])
gen_var_importance$colmod<-rownames(gen_var_importance)
#colnames(gen_var_importance)[1:3]<-c('Sky','Leaves','Stems')

gen_acc_results<-list()
gen_results<-list()
gen_tables<-list()
for (i in 2:length(test)){
  start<-Sys.time()
  print(i)
  Code<-names(test)[i]  
  samps_all<-terra::vect(samp_list[which(gsub('_SampleLocations','',tools::file_path_sans_ext(basename(samp_list)))==Code)])
  
  # check there is an image to classify....
  indx<-which(gsub("_lab","",tools::file_path_sans_ext(basename(all_colourmodels)))==Code)
  srgb_indx<-which(gsub("_srgb","",tools::file_path_sans_ext(basename(all_colourmodels)))==Code)
  if(length(indx)>0){
    print('image found')
  }else{
    print(paste0('error for i = ',i,' - no file in rgb_tiffs named ',Code))
  }
  
  #read in LAB image
  LAB_tiff<-rast(all_colourmodels[[indx]])
  srgb_tiff<-rast(all_colourmodels[[srgb_indx]])
  #pprgb_tiff<-rast(all_colourmodels[[pprgb_indx]])
  #xyz_tiff<-rast(all_colourmodels[[XYZ_indx]])
  
  names(LAB_tiff)<-c("L","A","B")
  names(srgb_tiff)<-c("red","green","blue")
  #names(pprgb_tiff)<-c("Red_PP","Green_PP","Blue_PP")
  #names(xyz_tiff)<-c("X","Y","Z")
  
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
  
  whichmax_img <- terra::predict(combo_tif, classifier_RF, na.rm = T,type='response')
  #plot(whichmax_img)
  #writeRaster(whichmax_img,paste0(end_path,"gen_classified_imgs_22Nov24/",Code,"_individ_seg_classified.tif"),overwrite=T)
  
  samps<-samps_all[j == k,]
  gref<-dummy_georef(whichmax_img)
  rm(whichmax_img)
  samps$predicted<-extract(gref,samps)[,2]
  samps$predicted<-as.character(class_to_factor(samps$predicted))
  sampsdf<-na.omit(as.data.frame(samps))
  colnames(sampsdf)[1]<-'observed'
  
  gen_acc_results[[i]]<-data.frame(accuracy=(sum(sampsdf$predicted==sampsdf$observed)/length(sampsdf$observed)),Code=Code)
  gen_tables[[i]]<-table(sampsdf$predicted,sampsdf$observed)
  names(gen_tables)[i]<-Code
  sampsdf$Code<-Code
  gen_results[[i]]<-sampsdf
  if(gen_acc_results[[i]]$accuracy<0.85){
    print(paste0(Code,' still sucks ',gen_acc_results[[i]]$accuracy))
  }else{
    print(paste0(Code,' success ',gen_acc_results[[i]]$accuracy))
  }
  end<-Sys.time()
  print(end-start)
}
saveRDS(do.call(rbind,gen_acc_results),paste0(end_path,'Generalised_classifier_accuracy.rds'))
saveRDS(do.call(rbind,gen_results),paste0(end_path,'Generalised_classifier_results.rds'))
saveRDS(gen_tables,paste0(end_path,'Generalised_classifier_tables.rds'))




SDF2 <- SDF[which(unlist(lapply(SDF,function(x)unique(x[,which(colnames(x)=='Sky')])=='clear')==T))]
#rm(sample_dataframe)

###################################################################################################
### General clear RF classifier ################################################################
###################################################################################################
###### create groupings for test/train
set.seed(42)
# number of folds
k <- 2
# randomly assign test or train
j<-sample(rep(1:k,each=round(length(SDF2))/k))
#splitting test/train dataset
train <- do.call(rbind,SDF2[which(j!= k)])
test <- SDF2[which(j== k)]
# specify variables to use in the model
pixel_vars<-which(colnames(train)%in%c('L','A','B','Hue','Chroma','red','green','blue','Hue_sRGB','Saturation_sRGB','Value_sRGB'))

# train the classifier
classifier_RF = randomForest(x = train[,pixel_vars], 
                             y = train$Class, 
                             ntree = 200,mtree=2)
## add something that captures 'importance'
gen_var_importance<-as.data.frame(importance(classifier_RF)[order(importance(classifier_RF)[,1],decreasing=T),drop=F,])
gen_var_importance$Code<-Code
gen_var_importance$colmod<-rownames(gen_var_importance)
colnames(gen_var_importance)[1:3]<-c('Sky','Leaves','Stems')

gen_acc_results<-list()
gen_results<-list()
gen_tables<-list()
for (i in 1:length(test)){
  start<-Sys.time()
  print(i)
  Code<-names(test)[i]  
  samps_all<-terra::vect(samp_list[which(gsub('_SampleLocations','',tools::file_path_sans_ext(basename(samp_list)))==Code)])
  
  # check there is an image to classify....
  indx<-which(gsub("_lab","",tools::file_path_sans_ext(basename(all_colourmodels)))==Code)
  srgb_indx<-which(gsub("_srgb","",tools::file_path_sans_ext(basename(all_colourmodels)))==Code)
  if(length(indx)>0){
    print('image found')
  }else{
    print(paste0('error for i = ',i,' - no file in rgb_tiffs named ',Code))
  }
  
  #read in LAB image
  LAB_tiff<-rast(all_colourmodels[[indx]])
  srgb_tiff<-rast(all_colourmodels[[srgb_indx]])
  #pprgb_tiff<-rast(all_colourmodels[[pprgb_indx]])
  #xyz_tiff<-rast(all_colourmodels[[XYZ_indx]])
  
  names(LAB_tiff)<-c("L","A","B")
  names(srgb_tiff)<-c("red","green","blue")
  #names(pprgb_tiff)<-c("Red_PP","Green_PP","Blue_PP")
  #names(xyz_tiff)<-c("X","Y","Z")
  
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
  
  whichmax_img <- terra::predict(combo_tif, classifier_RF, na.rm = T,type='response')
  #plot(whichmax_img)
  #writeRaster(whichmax_img,paste0(end_path,"gen_classified_imgs_22Nov24/",Code,"_individ_seg_classified.tif"),overwrite=T)
  
  samps<-samps_all[j == k,]
  gref<-dummy_georef(whichmax_img)
  rm(whichmax_img)
  samps$predicted<-extract(gref,samps)[,2]
  samps$predicted<-as.character(class_to_factor(samps$predicted))
  sampsdf<-na.omit(as.data.frame(samps))
  colnames(sampsdf)[1]<-'observed'
  
  gen_acc_results[[i]]<-data.frame(accuracy=(sum(sampsdf$predicted==sampsdf$observed)/length(sampsdf$observed)),Code=Code)
  gen_tables[[i]]<-table(sampsdf$predicted,sampsdf$observed)
  names(gen_tables)[i]<-Code
  sampsdf$Code<-Code
  gen_results[[i]]<-sampsdf
  print(end-start)
  print(i)
  if(gen_acc_results[[i]]$accuracy<0.85){
    print(paste0(Code,' still sucks ',gen_acc_results[[i]]$accuracy))
  }else{
    print(paste0(Code,' success ',gen_acc_results[[i]]$accuracy))
  }
  end<-Sys.time()
  print(end-start)
}
saveRDS(do.call(rbind,gen_acc_results),paste0(end_path,'Generalised_clear_classifier_accuracy.rds'))
saveRDS(do.call(rbind,gen_results),paste0(end_path,'Generalised_clear_classifier_results.rds'))
saveRDS(gen_tables,paste0(end_path,'Generalised_clear_classifier_tables.rds'))



SDF2 <- SDF[which(unlist(lapply(SDF,function(x)unique(x[,which(colnames(x)=='Sky')])=='overcast')==T))]
#rm(sample_dataframe)

###################################################################################################
### General overcast RF classifier ################################################################
###################################################################################################
###### create groupings for test/train
set.seed(42)
# number of folds
k <- 2
# randomly assign test or train
j<-sample(rep(1:k,each=round(length(SDF2))/k))
#splitting test/train dataset
train <- do.call(rbind,SDF2[which(j!= k)])
test <- SDF2[which(j== k)]
# specify variables to use in the model
pixel_vars<-which(colnames(train)%in%c('L','A','B','Hue','Chroma','red','green','blue','Hue_sRGB','Saturation_sRGB','Value_sRGB'))

# train the classifier
classifier_RF = randomForest(x = train[,pixel_vars], 
                             y = train$Class, 
                             ntree = 200,mtree=2)
## add something that captures 'importance'
gen_var_importance<-as.data.frame(importance(classifier_RF)[order(importance(classifier_RF)[,1],decreasing=T),drop=F,])
gen_var_importance$Code<-Code
gen_var_importance$colmod<-rownames(gen_var_importance)
colnames(gen_var_importance)[1:3]<-c('Sky','Leaves','Stems')

gen_acc_results<-list()
gen_results<-list()
gen_tables<-list()
for (i in 1:length(test)){
  start<-Sys.time()
  print(i)
  Code<-names(test)[i]  
  samps_all<-terra::vect(samp_list[which(gsub('_SampleLocations','',tools::file_path_sans_ext(basename(samp_list)))==Code)])
  
  # check there is an image to classify....
  indx<-which(gsub("_lab","",tools::file_path_sans_ext(basename(all_colourmodels)))==Code)
  srgb_indx<-which(gsub("_srgb","",tools::file_path_sans_ext(basename(all_colourmodels)))==Code)
  if(length(indx)>0){
    print('image found')
  }else{
    print(paste0('error for i = ',i,' - no file in rgb_tiffs named ',Code))
  }
  
  #read in LAB image
  LAB_tiff<-rast(all_colourmodels[[indx]])
  srgb_tiff<-rast(all_colourmodels[[srgb_indx]])
  #pprgb_tiff<-rast(all_colourmodels[[pprgb_indx]])
  #xyz_tiff<-rast(all_colourmodels[[XYZ_indx]])
  
  names(LAB_tiff)<-c("L","A","B")
  names(srgb_tiff)<-c("red","green","blue")
  #names(pprgb_tiff)<-c("Red_PP","Green_PP","Blue_PP")
  #names(xyz_tiff)<-c("X","Y","Z")
  
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
  
  whichmax_img <- terra::predict(combo_tif, classifier_RF, na.rm = T,type='response')
  #plot(whichmax_img)
  #writeRaster(whichmax_img,paste0(end_path,"gen_classified_imgs_22Nov24/",Code,"_individ_seg_classified.tif"),overwrite=T)
  
  samps<-samps_all[j == k,]
  gref<-dummy_georef(whichmax_img)
  rm(whichmax_img)
  samps$predicted<-extract(gref,samps)[,2]
  samps$predicted<-as.character(class_to_factor(samps$predicted))
  sampsdf<-na.omit(as.data.frame(samps))
  colnames(sampsdf)[1]<-'observed'
  
  gen_acc_results[[i]]<-data.frame(accuracy=(sum(sampsdf$predicted==sampsdf$observed)/length(sampsdf$observed)),Code=Code)
  gen_tables[[i]]<-table(sampsdf$predicted,sampsdf$observed)
  names(gen_tables)[i]<-Code
  sampsdf$Code<-Code
  gen_results[[i]]<-sampsdf
  print(end-start)
  print(i)
  if(gen_acc_results[[i]]$accuracy<0.85){
    print(paste0(Code,' still sucks ',gen_acc_results[[i]]$accuracy))
  }else{
    print(paste0(Code,' success ',gen_acc_results[[i]]$accuracy))
  }
  end<-Sys.time()
  print(end-start)
}
saveRDS(do.call(rbind,gen_acc_results),paste0(end_path,'Generalised_overcast_classifier_accuracy.rds'))
saveRDS(do.call(rbind,gen_results),paste0(end_path,'Generalised_overcast_classifier_results.rds'))
saveRDS(gen_tables,paste0(end_path,'Generalised_overcast_classifier_tables.rds'))












#################################################
### Figures
#################################################
end_path<-"F:/Training tiffs/"
titles<-c('All','Clear','Overcast')
subs<-c('_','_clear_','_overcast_')
acc_res_all<-list()
for(mmm in 1:length(titles)){
  acc_res_all[[mmm]]<-readRDS(paste0(end_path,'Generalised',subs[mmm],'classifier_accuracy.rds'))
  acc_res_all[[mmm]]$Sky<-titles[mmm]
}
acc_res_all<-do.call(rbind,acc_res_all)


plot_list<-list()

plot_list[[4]]<-ggplot(data=acc_res_all,aes(x=Sky,y=accuracy*100,fill=Sky))+
  geom_boxplot(outlier.size=2,outlier.fill = 'black',alpha=0.9,linewidth=.7,varwidth = TRUE)+theme_bw()+ # varwidth should change the width of each box proportional to the number of obs...
  scale_fill_manual(values=c("white","#56B4E9","#999999"))+scale_y_continuous(limits=c(30,103),position='right')+
  ggtitle(paste0('Total'),)+
  theme(axis.text.x=element_text(size=txt.size,angle = 90, hjust = 1,vjust = 0.5,colour='black'),
        axis.text.y=element_text(colour='black'),
        text = element_text(size=txt.size+2),
        legend.position = "none",
        plot.title = element_text(size = txt.size+2,hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title=element_text(size=15))+
  labs(x='',y='Accuracy (%)')+  scale_x_discrete(labels = label_wrap(10))



for(mmm in 1:length(titles)){
  acc_results<-readRDS(paste0(end_path,'Generalised',subs[mmm],'classifier_accuracy.rds'))
  results<-readRDS(paste0(end_path,'Generalised',subs[mmm],'classifier_results.rds'))
  tables<-readRDS(paste0(end_path,'Generalised',subs[mmm],'classifier_tables.rds'))
  
  print(titles[mmm])
  print(mean(acc_results$accuracy))
  print(standard_error(acc_results$accuracy))
  

  #######################################################################################
  # proportion tables
  txt.size=12
  prop_sum_table<-list()
  for (i in 1:length(tables)){
    if(!is.null(tables[[i]])){
        x<-tables[[i]]
        
        for (j in 1:nrow(x)){
          x[j,]<-round(x[j,]/sum(x[j,]),digits=2)
        }
        prop_sum_table[[i]]<-x
    }
  }
  
  # get an average from all tables
  combos<-matrix(nrow=0,ncol=2)
  for (r in 1:length(levels(as.factor(colnames(tables[[1]]))))){
    for(c in 1:length(levels(as.factor(colnames(tables[[1]]))))){
      x<-cbind(levels(as.factor(colnames(tables[[1]])))[r],levels(as.factor(colnames(tables[[1]])))[c])
      combos<-rbind(combos,x)
    }
  }
  end_df<-as.data.frame(matrix(ncol=4,nrow=nrow(combos)))
  colnames(end_df)<-c('observed','predicted','proportion','count')
  end_df$count<-0
  end_df$proportion<-0
  for(k in 1:nrow(combos)){
    end_df[k,1]<-combos[k,1]
    end_df[k,2]<-combos[k,2]
  }
  for(i in 1:length(prop_sum_table)){
    if(!is.null(prop_sum_table[[i]])){
      for (k in 1:nrow(combos)){
        r<-which(rownames(prop_sum_table[[i]])==combos[k,1])
        c<-which(colnames(prop_sum_table[[i]])==combos[k,2])
        if(!identical(r, integer(0))&!identical(c, integer(0))){
          end_df[k,3]<-prop_sum_table[[i]][r,c]+end_df[k,3]
          end_df[k,4]<-end_df[k,4]+1
        }
      }
    }
  }
  
  
  end_df$mean.proportion<-end_df$proportion/end_df$count
  ######turn the above into a table
  tab<-as.data.frame(matrix(nrow = length(levels(as.factor(end_df$observed))),ncol=length(levels(as.factor(end_df$observed)))))
  colnames(tab)<-levels(as.factor(end_df$observed))
  rownames(tab)<-levels(as.factor(end_df$observed))
  for (r in 1:length(rownames(tab))){
    for(c in 1:length(colnames(tab))){
      for(i in 1:nrow(end_df)){
        if(end_df[i,1]==rownames(tab)[r]&end_df[i,2]==colnames(tab)[c]){
          tab[r,c]<-end_df$mean.proportion[i]
        }
      }
    }
  }
  
  
  title<-paste0(titles[mmm],' (n = ',nrow(acc_results),')')
  hm<-end_df[,c(1,2,5)]
  hm$mean.proportion<-round(hm$mean.proportion,2)
  txt.size<-15
  
  plot_list[[mmm]]<-ggplot(hm, aes(x=predicted, y=observed)) + 
    ggtitle(paste0(title),)+
    geom_tile(aes(fill=mean.proportion)) +
    theme_bw()+geom_text(aes(label=mean.proportion))+
    scale_fill_gradientn(colours = terrain.colors(10,rev = TRUE),name = "Accuracy")+
    theme(text = element_text(size=txt.size),
          axis.text = element_text(size=txt.size),
          axis.text.x=element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = 'none',
          plot.title = element_text(size = txt.size),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.title=element_text(size=txt.size))
  
} 
((plot_list[[1]] / plot_list[[2]] /plot_list[[3]] + plot_layout(guides = 'collect')) | wrap_elements(plot_list[[4]]) + 
    plot_layout(guides = 'collect',widths = c(2, 1)))

ggsave(paste0(end_path,'Plots/Alt Review/',"Sky_accuracy_plots_gen_class.tiff"), width = 28.5, height = 22.5, units = "cm")
