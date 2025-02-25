# Install and load required packages if not already installed
library(terra)
library(sf)
library(rpart)
library(class)
library(fishmethods)
library(randomForest) 
library(reshape2)
library(ggplot2)
library(viridis)
library(cowplot)
library(scales)

segment<-function(img,
                  xc=2049, # x coordinate of mask centre
                  yc=2970, # y coordinate of mask centre
                  rc=1920, # radius of mask
                  nseg = 8,  # Number of azimuth segments
                  startVZA=0, # start viewing zenith angle
                  endVZA=90, # end viewing zenith angle
                  maxVZA=90, # maximum view angle of the camera
                  nrings=45, # n-1 for the number of zenith segments
                  add2img=F, # should the output be a dataframe or added to the image? default is false 
                  reduce_outputs=T,  # just keep the necessary outputs (discard intermediate steps)
                  x_col='X',  # x column name
                  y_col='Y'  # y column name
){
  # determine whether it's a classified image or a dataframe
  if(class(img)[1]=="SpatRaster"){
    imgdf<-terra::as.data.frame(img, xy=TRUE)
    base::names(imgdf)<-c('x','y',names(img))
  }
  if(class(img)[1]=='data.frame'){
    imgdf<-img
    colnames(imgdf)[which(colnames(imgdf)==x_col)]<-'x'
    colnames(imgdf)[which(colnames(imgdf)==y_col)]<-'y'
  }
  ###################################################
  
  VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
  maxVZAmm<-2*sin(maxVZA*pi/180)
  rset<-round(rc*2*sin(VZAbins*pi/180)/maxVZAmm) 
  
  
  
  #create radial distance from centre
  imgdf$dx<-imgdf$x-xc
  imgdf$dy<-imgdf$y-yc
  imgdf$r<-round(sqrt(imgdf$dx^2+imgdf$dy^2)) # radial distance of each pixel to the centre
  
  
  # Initialize a new column 'azimuth' in the data frame 'imgdf' and set all its values to NA
  imgdf$azimuth = NA
  
  # Calculate angles for points in the data frame based on differences in x and y coordinates
  # For points in the first quadrant (dx > 0 and dy >= 0)
  imgdf$azimuth[imgdf$dx > 0 & imgdf$dy >= 0] <- atan(imgdf$dy[imgdf$dx > 0 & imgdf$dy >= 0] / imgdf$dx[imgdf$dx > 0 & imgdf$dy >= 0])
  
  # For points in the second quadrant (dx > 0 and dy < 0)
  imgdf$azimuth[imgdf$dx > 0 & imgdf$dy < 0] <- atan(imgdf$dy[imgdf$dx > 0 & imgdf$dy < 0] / imgdf$dx[imgdf$dx > 0 & imgdf$dy < 0]) + 2 * pi
  
  # For points in the third quadrant (dx < 0)
  imgdf$azimuth[imgdf$dx < 0] <- atan(imgdf$dy[imgdf$dx < 0] / imgdf$dx[imgdf$dx < 0]) + pi
  
  # For points on the positive y-axis (dx == 0 and dy > 0)
  imgdf$azimuth[imgdf$dx == 0 & imgdf$dy > 0] <- pi / 2
  
  # For points on the negative y-axis (dx == 0 and dy < 0)
  imgdf$azimuth[imgdf$dx == 0 & imgdf$dy < 0] <- 3 * pi / 2
  
  # Convert angles from radians to degrees
  imgdf$azimuth <- imgdf$azimuth * 180 / pi
  
  # Create a new column 'ring' by cutting the 'r' column into segments
  imgdf$ring <- cut(imgdf$r, rset, include.lowest = TRUE,
                    labels = seq(startVZA + ((endVZA - startVZA) / 2 / nrings),
                                 endVZA - ((endVZA - startVZA) / 2 / nrings),
                                 (endVZA - startVZA) / nrings))
  
  if(add2img==F){
    # Drop rows with missing values in the 'ring' column
    imgdf <- imgdf[!is.na(imgdf$ring), ]
    
    # Convert 'ring' to numeric
    imgdf$ring <- as.numeric(as.character(imgdf$ring))
    
    # Create a new column 'azimuth_group' by cutting the 'azimuth' column into segments
    imgdf$azimuth_group <- cut(imgdf$azimuth, seq(0, 2 * pi, 2 * pi / nseg) * 180 / pi,
                             include.lowest = TRUE,
                             labels = rep(seq(2 * pi / nseg, 2 * pi, 2 * pi / nseg)) * 180 / pi)
    
    # Convert 'alpha.to' to numeric
    imgdf$azimuth_group <- as.numeric(as.character(imgdf$azimuth_group))
    
    #reduce outputs
    if(reduce_outputs==T){
      idx<-which(colnames(imgdf)%in%c('Leaves','Sky','Stems','class','r','azimuth','ring','azimuth_group'))
      imgdf<-imgdf[,idx]
    }
    
    return(imgdf)
  }
  if(add2img==T){
    # Convert 'ring' to numeric
    imgdf$ring <- as.numeric(as.character(imgdf$ring))
    
    # Create a new column 'azimuth_group' by cutting the 'azimuth' column into segments
    imgdf$azimuth_group <- cut(imgdf$azimuth, seq(0, 2 * pi, 2 * pi / nseg) * 180 / pi,
                             include.lowest = TRUE,
                             labels = rep(seq(2 * pi / nseg, 2 * pi, 2 * pi / nseg)) * 180 / pi)
    
    # Convert 'alpha.to' to numeric
    imgdf$azimuth_group <- as.numeric(as.character(imgdf$azimuth_group))
    
    # Convert 'azimuth_group' to a raster layer
    idx_1<-which(colnames(imgdf)%in%c('azimuth_group','ring','x','y'))
    img$r <- imgdf$r # this works do for all check where the rgb vals have gone wtf...
    img$ring<-imgdf$ring
    img$azimuth<-imgdf$azimuth
    img$azimuth_group<-imgdf$azimuth_group
    return(img)
  }
}


toNumerics <- function(Date) {                  ######function to seperate out numbers in date strings
  stopifnot(inherits(Date, c("Date", "POSIXt")))
  day <- as.numeric(strftime(Date, format = "%d"))
  month <- as.numeric(strftime(Date, format = "%m"))
  year <- as.numeric(strftime(Date, format = "%Y"))
  hour <- as.numeric(strftime(Date, format = "%H.%M"))
  list(year = year, month = month, day = day,hour = hour)
}

nor <-function(x) { (x -min(x))/(max(x)-min(x))   }
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

gamma_cor<-function(rgb,gamma=2.2){
  for(i in 1:3){
  minmB=min(na.omit(terra::values(rgb[[i]])))
  maxmB=max(na.omit(terra::values(rgb[[i]])))
  terra::values(rgb[[i]])=(maxmB-minmB)*((terra::values(rgb[[i]])/(maxmB-minmB))^gamma)
}
  return(rgb)
}

add_details<-function(img_details,image_list){
  img_details<-merge(data.frame(Code=gsub("_georef","",tools::file_path_sans_ext(basename(georef_rgb_path)))),img_details,by='Code',all.x=T)
  for (i in 1:nrow(img_details)){
    if(is.na(img_details$Site[i])){
      idx<-which(image_list$Code==img_details$Code[i])
      img_details[i,3:4]<-image_list[idx,2:3]
      img_details$file_type[i]<-'ARW'
    }
    if(is.na(img_details$Sky[i])){
      idx<-which(image_list$Code==img_details$Code[i])
      img_details[i,]$Sky<-image_list[idx,]$Sky
    }
  }
  return(img_details)
  write.csv(img_details,paste0(end_path,"img_details.csv"),row.names = F,na='')
}

dummy_georef<-function(raster_img){
  # basically making a spatial object
  geom<-matrix(data=c(116.286427,-29.440805)) # centred on perenjori hotel
  middle = st_sf( index = 1, geom = st_sfc(st_point(geom[,1])))
  st_crs(middle)<-'EPSG:4326' #set crs = coordinate reference system
  
  # change the crs to something that uses metres
  trans_middle<-st_transform(middle,crs='EPSG:3857')
  metres<-matrix(data=st_coordinates(trans_middle),nrow=1)
  
  #list the extent (multiplied by 100) of the plot
  ########### note #################
  # extent is 50 m by 50 m not 50 x 50 cm
  
  
  xmin <- metres[1,1]-ncol(raster_img)/2  # Adjust if needed
  xmax <- metres[1,1]+ncol(raster_img)/2  # Adjust
  ymin <- metres[1,2]-nrow(raster_img)/2  # Adjust
  ymax <- metres[1,2]+nrow(raster_img)/2   # Adjust
  
  
  # create an extent object
  extent<-terra::ext(xmin,xmax,ymin,ymax)
  
  
  #### set the extent of the images
  terra::ext(raster_img)<-extent
  terra::crs(raster_img)<-'EPSG:3857'
  return(raster_img)
}
#sampdata,'X_coord','Y_coord'
undo_georef<-function(df=sampdata,x_col='X_coord',y_col='Y_coord',xcols=4024,ycols=6024){
# basically making a spatial object
geom<-matrix(data=c(116.286427,-29.440805)) # centred on perenjori hotel
middle = st_sf( index = 1, geom = st_sfc(st_point(geom[,1])))
st_crs(middle)<-'EPSG:4326' #set crs = coordinate reference system

# change the crs to something that uses metres
trans_middle<-st_transform(middle,crs='EPSG:3857')
metres<-matrix(data=st_coordinates(trans_middle),nrow=1)

xmin <- metres[1,1]-xcols/2  # Adjust if needed
ymin <- metres[1,2]-ycols/2  # Adjust

df[,which(colnames(df)==x_col)]<-round(df[,which(colnames(df)==x_col)]-xmin+0.5)
df[,which(colnames(df)==y_col)]<-round(df[,which(colnames(df)==y_col)]-ymin+0.5)

return(df)
}


# Function to convert LAB to RGB for plotting
lab_to_rgb <- function(lab_raster) {
  lab_raster$L <- terra::clamp(lab_raster$L, 0, 100)  # Clip to valid LAB ranges
  lab_raster$A <- terra::clamp(lab_raster$A, -128, 127)  # Clip to valid LAB ranges
  lab_raster$B <- terra::clamp(lab_raster$B, -128, 127)  # Clip to valid LAB ranges
  values(lab_raster)[,2:3] <- values(lab_raster)[,2:3]+128  # Shift to positive values
  values(lab_raster)[,2:3] <- values(lab_raster)[,2:3] / 255  # Scale to [0, 1]
  values(lab_raster)[,1] <- values(lab_raster)[,1] / 100  # Scale to [0, 1]
  values(lab_raster)<-values(lab_raster)*255
  terra::setMinMax(lab_raster)
  RGB(lab_raster)<-c(1,2,3)
  return(lab_raster)
}

mask<-function(rgb_raster,  
               xc=2049,
               yc=2970,
               rc=1920){
  # mask
  xy <- terra::xyFromCell(RGB_HSV_rast,1:terra::ncell(rgb_raster)) #create xy coordinates
  circular.mask = (xy[,1] - xc)^2 + (xy[,2] - yc)^2 <= rc^2 #determine which pixels are within the circular mask
  terra::values(rgb_raster)[!circular.mask] <- NA #everything outside of that is NA
}

#################################################################################
#################################################################################
#################################################################################
### creating the trained model for all images
#################################################################################
#################################################################################
#################################################################################

#path for all images
path1="F:/Oct 23 Hemi tiffs LAI paper/"
path2="F:/Field sites and data/Hemi cam/June 2023/"
rgb_tiffs<-Sys.glob(file.path(c(path1,path2), "*.tif"))

#training images
end_path="F:/Training tiffs/"
updated_shapes="F:/Training tiffs/edited_shapes/"

shapefiles<-Sys.glob(file.path(updated_shapes, "*.shp"))
georef_rgb_path<-Sys.glob(file.path(end_path, "*_georef.tiff"))
image_list<-read.csv(paste0(end_path,"img_ID_metadata.csv"))
img_details<-read.csv(paste0(end_path,"img_details.csv"))
img_details<-add_details(img_details,image_list)
#write.csv(img_details,paste0(end_path,"img_details.csv"),row.names = F,na='')
all_colourmodels<-c(Sys.glob(file.path("F:/PhotoproRGB_LAB_imgs",'*.tif')),Sys.glob(file.path("F:/PhotoproRGB_LAB_imgs",'*.tiff')))


###########
########### K fold test and train procedure
########### treatments: 
###########   sky mask/unmasked -> 
###########   Recursive Partitioning and Regression Trees/
########### 
pixval_list<-list()
sampdata_list<-list()

# redo_idx<-c(123,105,102,101,93,92,90,89,88,84,85,86,87,83,81,79,75,72,70,69,68,63,59,57,44,46,47,48,49,50,52,53,54,56,43,32,28,24,39,31,27,23,38,30,26,36,29,25,15)
# tools::file_path_sans_ext(basename(georef_rgb_path[redo_idx]))
# i=1
#problem in i = 38 DSC04847_classes

#4736 needs to be redone - 21.11.24
for (i in 1:length(georef_rgb_path)){
  idx<-which(img_details$Code==gsub("_georef","",tools::file_path_sans_ext(basename(georef_rgb_path[i]))))
  if(img_details$file_type[idx]=='ARW'){ 
    georef_rgb <- rast(georef_rgb_path[i])

    idx_d<-which(gsub("_classes","",tools::file_path_sans_ext(basename(shapefiles)))==gsub("_georef","",tools::file_path_sans_ext(basename(georef_rgb_path[i]))))
    Categorised_shp<-st_read(shapefiles[idx_d])
    #Categorised_shp$Class<-as.factor(Categorised_shp$Class)

    if(img_details$Sky[idx]=='overcast'&"Sun"%in%levels(Categorised_shp$Class)){ 
      print(paste0('Check for sun in ',shapefiles[idx_d]))
      break
      }

    # create a hsv image
    if(FALSE %in% as.logical(names(georef_rgb)==c('red','green','blue'))){
      names(georef_rgb)<-c('red','green','blue')
    }
    
    georef_rgb2<-georef_rgb/65535*255
    RGB(georef_rgb2)<-c(1,2,3)
    georef_hsv<-colorize(georef_rgb2,to='hsv')
    names(georef_hsv)<-c('Hue_sRGB','Saturation_sRGB','Value_sRGB')

    
    ### alternative colour models:
    indx<-which(gsub("_lab","",tools::file_path_sans_ext(basename(all_colourmodels)))==img_details$Code[idx])
    pprgb_indx<-which(gsub("_PPrgb","",tools::file_path_sans_ext(basename(all_colourmodels)))==img_details$Code[idx])
    XYZ_indx<-which(gsub("_xyz","",tools::file_path_sans_ext(basename(all_colourmodels)))==img_details$Code[idx])
    
    #read in LAB image
    LAB_tiff<-dummy_georef(rast(all_colourmodels[[indx]]))
    pprgb_tiff<-dummy_georef(rast(all_colourmodels[[pprgb_indx]]))
    xyz_tiff<-dummy_georef(rast(all_colourmodels[[XYZ_indx]]))
    
    names(LAB_tiff)<-c("L","A","B")
    names(pprgb_tiff)<-c("Red_PP","Green_PP","Blue_PP")
    names(xyz_tiff)<-c("X","Y","Z")
    
    #add Chroma and Hue to image
    LAB_tiff$Chroma<-sqrt(LAB_tiff$A^2 + LAB_tiff$B^2)
    LAB_tiff$Hue <- (atan2(LAB_tiff$B, LAB_tiff$A) * 180 / pi + 360) %% 360
    

    #Fix regular naming issues.
    if ('Insect' %in% levels(Categorised_shp$Class)){
      levels(Categorised_shp$Class)<-c(levels(Categorised_shp$Class),"Insects")
      Categorised_shp$Class[Categorised_shp$Class=='Insect']<-"Insects"
      Categorised_shp$Class<-droplevels(Categorised_shp$Class)
    }
    if ('Senesced leaves' %in% levels(Categorised_shp$Class)){
      levels(Categorised_shp$Class)<-c(levels(Categorised_shp$Class),"Brown leaves")
      Categorised_shp$Class[Categorised_shp$Class=='Senesced leaves']<-"Brown leaves"
      Categorised_shp$Class<-droplevels(Categorised_shp$Class)
    }
    
    #Delete anything with a null geometry 
    if(TRUE %in% st_is_empty(Categorised_shp$geometry)){
      Categorised_shp<-Categorised_shp[!st_is_empty(Categorised_shp$geometry),]
      #Categorised_shp$Class<-droplevels(Categorised_shp$Class)
      st_write(Categorised_shp,shapefiles[idx_d],append=F)
      print("deleted empty vector layer")
    }
    
    # Check there are no additional odd categories.
    if(FALSE %in% as.logical(levels(Categorised_shp$Class) %in% c("Dark wood","Leaves","Light wood","Sky","Sun","Flowers","Insects","Brown leaves"))){
      print(paste0('ERROR: odd category in ',basename(shapefiles[idx_d])))
      break
    }
    if(TRUE %in% as.logical(is.na(Categorised_shp$Class)&is.na(Categorised_shp$solar_disk))){
      print(paste0('Null category in ',basename(shapefiles[idx_d])))
      break
    }
    #Categorised_shp$Class<-as.factor(as.character(Categorised_shp$Class))
    #Categorised_shp$Class[Categorised_shp$Class=="Dark wood\\"]<-"Dark wood" 
    #Categorised_shp$Class<-droplevels(Categorised_shp$Class)
    #st_write(Categorised_shp,shapefiles[idx_d],append=F,SHAPE_RESTORE_SHX=T)
    
    Categorised_shp<-Categorised_shp[!is.na(Categorised_shp$Class),]
    Categorised_shp$Class<-as.factor(Categorised_shp$Class)
    
    # sample locations within each class
    grid = rast(Categorised_shp, nrow=nrow(georef_rgb), ncol=ncol(georef_rgb))
    ncr = terra::rasterize(Categorised_shp, grid, field="Class")
    #plot(ncr)
    
    # check for levels with too few values
    if(T %in% as.logical(as.data.frame(table(terra::as.data.frame(ncr)$Class))$Freq<1000)){
      print(paste0('ERROR: not enough values in ',basename(shapefiles[idx_d])," ",
                   levels(terra::as.data.frame(ncr)$Class)[which(as.logical(as.data.frame(table(terra::as.data.frame(ncr)$Class))$Freq<500)==T)]))
    }
    ########################################################################################
    samples <- spatSample(ncr, size = 100, method="stratified",as.points=T,xy=T)
    #help("spatSample")
    #help("writeRaster")
    # use names rather than levels for the sample group
    samples$Class<-as.factor(samples$Class)
    levels(samples$Class)<-levels(Categorised_shp$Class)
    #extract spectral values from each image
    df <- terra::extract(c(georef_rgb,georef_hsv,LAB_tiff,pprgb_tiff,xyz_tiff), samples,ID=F)
    sampdata <- data.frame(Class = samples$Class, df,x_coord=samples$x,y_coord=samples$y)
    sampdata$Code <- gsub("_georef",tools::file_path_sans_ext(basename(georef_rgb_path[i])),replacement = "")
    pixval_list[[i]]<-sampdata
    #########################################################################################
    Categorised_shp$Class<-as.character(Categorised_shp$Class)
    Categorised_shp$Min_class<-as.factor(Categorised_shp$Class)
    Categorised_shp$Maj_class<-NA
    for (j in 1:nrow(Categorised_shp)){
      if(Categorised_shp$Class[j]=='Leaves'){
        Categorised_shp$Maj_class[j]<-'Leaves'
      }
      if(Categorised_shp$Class[j] %in% c("Dark wood","Light wood")){
        Categorised_shp$Maj_class[j]<-'Stems'
      }
      if(Categorised_shp$Class[j] %in% c("Sky","Sun")){
        Categorised_shp$Maj_class[j]<-'Sky'
      }
    }
    Categorised_shp<-Categorised_shp[!as.logical(is.na(Categorised_shp$Maj_class)),]
    Categorised_shp$Maj_class<-as.factor(Categorised_shp$Maj_class)
    #levels(Categorised_shp$Maj_class)
    
    ncr = terra::rasterize(Categorised_shp, grid, field="Maj_class")
    set.seed(42)
    samples <- spatSample(ncr$Maj_class, size = 1000, method="stratified",values=T,as.points=T,xy=T,replace=F) # removed xy=T
    
    # add the minor class
    #ncr_2 = as.data.frame(terra::rasterize(Categorised_shp, grid, field="Min_class"),xy=T)
    #samples<-merge(samples,ncr_2,all.x=T)
    
    # use names rather than levels for the sample group
    samples$Class<-as.factor(samples$Maj_class)
    levels(samples$Class)<-levels(ncr$Maj_class)[[1]][,2]
    
    

    #extract spectral and zenith values from each image
    df <- terra::extract(c(georef_rgb,georef_hsv,LAB_tiff,pprgb_tiff,xyz_tiff), samples,ID=T)
    sampdata <- data.frame(Class = samples$Class, df,X_coord = samples$x,Y_coord = samples$y)
    samp_vec<-terra::vect(sampdata,geom=c("X_coord", "Y_coord"),crs=crs(samples))
    
    # plot(georef_rgb2)
    # plot(samp_vec,col=as.numeric(samp_vec$Class)+1,add=T,pch=1)
    # plot(ncr,add=T,col=alpha(unique(as.numeric(samp_vec$Class)+1),0.7))
    #terra::writeVector(samp_vec,paste0(end_path,'Sample locations/',gsub("_georef",tools::file_path_sans_ext(basename(georef_rgb_path[i])),replacement = ""),'_SampleLocations.shp'),overwrite=T)
    sampdata<-undo_georef(sampdata)
    
    sampdata<-segment(sampdata,
                      xc=2049, # x coordinate of mask centre
                      yc=2970, # y coordinate of mask centre
                      rc=1920, # radius of mask
                      nseg = 8,  # Number of azimuth segments
                      startVZA=0, # start viewing zenith angle
                      endVZA=90, # end viewing zenith angle
                      maxVZA=90, # maximum view angle of the camera
                      nrings=45, # n-1 for the number of zenith segments
                      add2img=T, # should the output be a dataframe or added to the image? default is false 
                      reduce_outputs=T,  # just keep the necessary outputs (discard intermediate steps)
                      x_col='X_coord',  # x column name
                      y_col='Y_coord'  # y column name
    )
    
    sampdata$Code <- gsub("_georef",tools::file_path_sans_ext(basename(georef_rgb_path[i])),replacement = "")
    sampdata_list[[i]]<-sampdata
    
    samp_vec<-cbind(samp_vec,sampdata[,which(colnames(sampdata)%in%c("X_coord","Y_coord","r","ring","azimuth","azimuth_group","Code"))])
    terra::writeVector(samp_vec,paste0(end_path,'Sample locations/',gsub("_georef",tools::file_path_sans_ext(basename(georef_rgb_path[i])),replacement = ""),'_SampleLocations.shp'),overwrite=T)

    #########################################################################################
    print(i)
  
  }# ARW if function....

}         ################### main pixel sampling loop
sample_dataframe<-do.call(rbind,sampdata_list)
pixval_dataframe<-do.call(rbind,pixval_list)

###############################################################
# add astro data
###############################################################
astro_data<-read.csv("F:/Training tiffs/Astro_details.csv")
astro_data<-astro_data[,which(colnames(astro_data)%in%c('year','month','day','hour','Code','file_type','Site','Trap','Long','Lat','Solar_azimuth','Solar_zenith','PAR'))]
# samples
sample_dataframe<-merge(sample_dataframe,astro_data,by='Code')
sample_dataframe$file_type<-as.factor(sample_dataframe$file_type)
sample_dataframe$Class<-as.factor(sample_dataframe$Class)

# pixels
pixval_dataframe<-merge(pixval_dataframe,astro_data,by='Code')
pixval_dataframe$file_type<-as.factor(pixval_dataframe$file_type)
pixval_dataframe$Class<-as.factor(pixval_dataframe$Class)

# save so we dont have to repeat above.
#write.csv(sample_dataframe,paste0("H:/My Documents/PhD 2020/Field sites and data/LAI validation chapter compiled/",'Pixel_value_training_dataset.csv'),na='',row.names = F)

saveRDS(sample_dataframe,paste0(end_path,'Pixel_value_training_dataset_v4.rds'))
saveRDS(pixval_dataframe, paste0(end_path,'Pixel_value_for_visualisation_v4.rds'))

# end
#sampdata_list<-readRDS(sampdata_listpaste0(end_path,'Pixel_value_training_dataset_v4.rds'))
#pixval_list<-readRDS(paste0(end_path,'Pixel_value_for_visualisation_v4.rds'))




























































































#sample_dataframe<-read.csv(paste0(end_path,'Pixel_value_training_dataset'))
levels(sample_dataframe$Class)<-c(levels(sample_dataframe$Class),'Stems')
sample_dataframe$Class[sample_dataframe$Class=='Sun']<-'Sky'
sample_dataframe$Class[sample_dataframe$Class=='Light wood']<-'Stems'
sample_dataframe$Class[sample_dataframe$Class=='Dark wood']<-'Stems'
sample_dataframe<-sample_dataframe[!sample_dataframe$Class%in%c('Flowers','Insects','Brown leaves'),]
sample_dataframe$Class<-droplevels(sample_dataframe$Class)

sample_dataframe<-sample_dataframe[sample_dataframe$file_type=='ARW',]
sample_dataframe<-sample_dataframe[sample_dataframe$Sky=='clear',]

####################################################################
#run test and train procedure for the random forest model ##########
####################################################################
# create groupings for test/train
# new procedure groups images rather than pixels
set.seed(42)
# number of folds
k <- 5
j <- sample(rep(1:k, each = round(length(unique(sample_dataframe$Code)))/k))
p<-cbind(unique(sample_dataframe$Code),j) #will give a warning message as it's repeating j at the end
colnames(p)[1]<-"Code"
p<-merge(p,sample_dataframe,by='Code')
j<-as.numeric(p$j)

####################################################################
#run test and train procedure for the k nearest neighbor model #####
####################################################################


#train[which(TRUE == is.na(train[,pixel_vars])),]
x_rf <- list()
pixel_vars<-which(colnames(sample_dataframe)%in%c('hue','saturation','value'))
#for (k_ in 1:10) {
  train <- sample_dataframe[j!= k_, ]
  train$Class<-droplevels(train$Class)
  test <- sample_dataframe[j == k_, ]
  classifier_RF = randomForest(x = train[,pixel_vars], 
                               y = train$Class, 
                               ntree = 500) 
  predicted <-predict(classifier_RF, newdata = test[,pixel_vars])
  x_rf <- cbind(predicted,test)
#  x_rf[[k_]] <- cbind(predicted,test)
#}

#add each run to a master df
#x_rf <- do.call(rbind, x_rf)
x_rf$algorithm<-'Random_forest'
colnames(x_rf)[3]<-'observed'
agreement(x_rf)


agrmt<-x_rf$observed==x_rf$predicted

# agreement by image - summary
agrmt_sum<-aggregate(agrmt,by=list(x_rf$Code),sum)
agrmt_length<-aggregate(agrmt,by=list(x_rf$Code),length)
agrmt_sum$prop<-agrmt_sum$x/agrmt_length$x
colnames(agrmt_sum)[c(1,2)]<-c("Code",'n correct')

# merge with image details
agrmt_sum<-merge(agrmt_sum,astro_data,by='Code',all.x=T,all.y=F)
agrmt_sum<-merge(agrmt_sum,img_details,by='Code',all.x=T,all.y=F)

# plot agreement....
with(agrmt_sum,plot(prop~zenith))
with(agrmt_sum,text(zenith,prop,agrmt_sum$Site.x))
with(agrmt_sum,text(zenith+2.5,prop,agrmt_sum$Trap.y))

summary(with(agrmt_sum,lm(prop~zenith)))
agrmt_sum$prop


  x<-table(x_rf$predicted,x_rf$observed)
  for (j in 1:nrow(x)){
    x[j,]<-round(x[j,]/sum(x[j,]),digits=2)
  }

save(classifier_RF,file = "F:/Training tiffs/classifier_RF.RData")
print(classifier_RF)
    
    
###########################################################
###########################################################
### proportion
###########################################################
###########################################################
  txt.size<-16
  hm<-melt(x)
  colnames(hm)<-c('observed','predicted','value')
  table_plot<-ggplot(hm, aes(x=predicted, y=observed)) + 
    ggtitle(paste0('Predicting new images'))+
    geom_tile(aes(fill=value)) +
    theme_bw()+geom_text(aes(label=value))+
    scale_fill_gradientn(colours = terrain.colors(10,rev = TRUE))+
    theme(text = element_text(size=txt.size-2),
          legend.position = 'none',
          plot.title = element_text(size = txt.size))
################################################
######## predicting for a given image ##########
################################################
agrmt_sum
for (i in 1:length(agrmt_sum)){
Code<-agrmt_sum$Code[i]
img<-rast(rgb_tiffs[[which(tools::file_path_sans_ext(basename(rgb_tiffs))==Code)]])
#add hsv to image
RGB(img)<-c(1,2,3)
hsv_img<-colorize(img, "hsv")
RGB_HSV_rast<-c(img,hsv_img)

#coordinates for mask (centre and radius)
xc=2049
yc=2970
rc=1920
xy <- terra::xyFromCell(RGB_HSV_rast,1:terra::ncell(RGB_HSV_rast)) #create xy coordinates
circular.mask = (xy[,1] - xc)^2 + (xy[,2] - yc)^2 <= rc^2 #determine which pixels are within the circular mask
terra::values(RGB_HSV_rast)[!circular.mask] <- NA #everything outside of that is NA
plot(RGB_HSV_rast)

# add zenith angle and file type
 #which(gsub("_masked","",tools::file_path_sans_ext(basename(masked_rgb_path)))
idx<-which(Code==astro_data$Code)
astro_data$file_type<-as.factor(astro_data$file_type)
RGB_HSV_rast[["zenith"]]<-astro_data[idx,]$zenith
RGB_HSV_rast[['file_type']]<-astro_data[idx,]$file_type

classified_RF <- terra::predict(RGB_HSV_rast, classifier_RF, na.rm = TRUE,type="prob")
writeRaster(classified_RF,paste0(end_path,Code,"_classified.tif"),overwrite=T)
}
#classified_RF[classified_RF$Stems>.8,]


#writeRaster(classified_RF,"F:/Training tiffs/classified.tiff")
#head(na.omit(terra::values(classified[["class"]])))
#df<-as.data.frame(levels(classified[["class"]]))
#df$count<-NA
#for (i in 1:nrow(df)){
#  df$count[i]<-sum(na.omit(values(classified[["class"]]))==i)
#}
#df

#KNN_test<-as.data.frame(as.matrix(nor(na.omit(terra::values(RGB_HSV_rast)[,1:6]))))
#KNN_test$zenith<-astro_data[idx,]$zenith
#KNN_test$file_type<-astro_data[idx,]$file_type
#classified_KNN <- knn(normalised_dataframe[,pixel_vars],KNN_test,cl=normalised_dataframe$Class,k=11)





