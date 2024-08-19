# Step 3 (After manual labelling):
# extract pixel values from each file for training and testing RF algorithm
# currently extracts gamma corrected and brightened RGB values, linear RGB, HSV, LABCH and XYZ spectral values
# read carefully: extracts 2 datasets based on categories and sub-categories of canopy components. this may need to be edited.


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

segment<-function(img,
                  xc=2049, # x coordinate of mask centre
                  yc=2970, # y coordinate of mask centre
                  rc=1920, # radius of mask
                  nseg = 8,  # Number of azimuth segments
                  startVZA=0, # start viewing zenith angle
                  endVZA=90, # end viewing zenith angle
                  maxVZA=90, # maximum view angle of the camera
                  nrings=45, # n-1 for the number of zenith segments
                  add2img=F # should the output be a dataframe or added to the image? default is false 
){
  # determine whether it's a classified image or a dataframe
  if(class(img)[1]=="SpatRaster"){
    imgdf<-terra::as.data.frame(img, xy=TRUE)
    base::names(imgdf)<-c('x','y',names(img))
  }
  if(class(img)[1]=='data.frame'){
    imgdf<-img
    colnames(imgdf)[which(colnames(imgdf)=='X')]<-'x'
    colnames(imgdf)[which(colnames(imgdf)=='Y')]<-'y'
  }
  ###################################################
  
  VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
  maxVZAmm<-2*sin(maxVZA*pi/180)
  rset<-round(rc*2*sin(VZAbins*pi/180)/maxVZAmm) 
  
  
  
  #create radial distance from centre
  imgdf$dx<-imgdf$x-xc
  imgdf$dy<-imgdf$y-yc
  imgdf$r<-round(sqrt(imgdf$dx^2+imgdf$dy^2)) # radial distance of each pixel to the centre
  
  
  # Initialize a new column 'theta' in the data frame 'imgdf' and set all its values to NA
  imgdf$theta = NA
  
  # Calculate angles for points in the data frame based on differences in x and y coordinates
  # For points in the first quadrant (dx > 0 and dy >= 0)
  imgdf$theta[imgdf$dx > 0 & imgdf$dy >= 0] <- atan(imgdf$dy[imgdf$dx > 0 & imgdf$dy >= 0] / imgdf$dx[imgdf$dx > 0 & imgdf$dy >= 0])
  
  # For points in the second quadrant (dx > 0 and dy < 0)
  imgdf$theta[imgdf$dx > 0 & imgdf$dy < 0] <- atan(imgdf$dy[imgdf$dx > 0 & imgdf$dy < 0] / imgdf$dx[imgdf$dx > 0 & imgdf$dy < 0]) + 2 * pi
  
  # For points in the third quadrant (dx < 0)
  imgdf$theta[imgdf$dx < 0] <- atan(imgdf$dy[imgdf$dx < 0] / imgdf$dx[imgdf$dx < 0]) + pi
  
  # For points on the positive y-axis (dx == 0 and dy > 0)
  imgdf$theta[imgdf$dx == 0 & imgdf$dy > 0] <- pi / 2
  
  # For points on the negative y-axis (dx == 0 and dy < 0)
  imgdf$theta[imgdf$dx == 0 & imgdf$dy < 0] <- 3 * pi / 2
  
  # Convert angles from radians to degrees
  imgdf$theta <- imgdf$theta * 180 / pi
  
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
    
    # Create a new column 'theta_group' by cutting the 'theta' column into segments
    imgdf$theta_group <- cut(imgdf$theta, seq(0, 2 * pi, 2 * pi / nseg) * 180 / pi,
                             include.lowest = TRUE,
                             labels = rep(seq(2 * pi / nseg, 2 * pi, 2 * pi / nseg)) * 180 / pi)
    
    # Convert 'alpha.to' to numeric
    imgdf$theta_group <- as.numeric(as.character(imgdf$theta_group))
    
    
    
    return(imgdf)
  }
  if(add2img==T){  #######################this is bronken you dumb shit
    # Convert 'ring' to numeric
    imgdf$ring <- as.numeric(as.character(imgdf$ring))
    
    # Create a new column 'theta_group' by cutting the 'theta' column into segments
    imgdf$theta_group <- cut(imgdf$theta, seq(0, 2 * pi, 2 * pi / nseg) * 180 / pi,
                             include.lowest = TRUE,
                             labels = rep(seq(2 * pi / nseg, 2 * pi, 2 * pi / nseg)) * 180 / pi)
    
    # Convert 'alpha.to' to numeric
    imgdf$theta_group <- as.numeric(as.character(imgdf$theta_group))
    
    # Convert 'theta_group' to a raster layer
    idx_1<-which(colnames(imgdf)%in%c('theta_group','ring','x','y'))
    img$r <- imgdf$r # this works do for all check where the rgb vals have gone wtf...
    img$ring<-imgdf$ring
    img$theta<-imgdf$theta
    img$theta_group<-imgdf$theta_group
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

#path for processed images
root_dir = "C:/Users/22064705/R directory/Fortescue hemi img processing/Submission scripts/"
start_path = paste0(root_dir,"training_labelled/")
processed_tiffs<-paste0(root_dir,'processed_tiffs/') # specify whether you want the autobright file HERE by replacing "*.tif" with "*bright.tif"


#training images
shapefiles<-Sys.glob(file.path(start_path, "*.shp"))
georef_rgb_path<-Sys.glob(file.path(start_path, "*_georef.tiff"))

# read in metadata
img_details<-read.csv(paste0(root_dir,"img_ID_metadata.csv"))


# set number of pixels to sample
npixels = 500

# set the names of the sub-categories used. light and dark wood are combined later in this script.
unique_cats = c("Dark wood","Leaves","Light wood","Sky","Sun","Flowers","Insects","Brown leaves")
# note: this creates 2 outputs: one that has the subcategory


# automated stuff
xyz_colourmodel<-Sys.glob(file.path(processed_tiffs,'*xyz.tif'))
lab_colourmodel<-Sys.glob(file.path(processed_tiffs,'*lab.tif'))
PPrgb_colourmodel<-Sys.glob(file.path(processed_tiffs,'*PPrgb.tif')) # uses linear not autobrightened tiff






pixval_list<-list()
sampdata_list<-list()
for (i in 1:length(georef_rgb_path)){
  
  # find the details for the image
  idx<-which(img_details$Code==gsub("_georef","",tools::file_path_sans_ext(basename(georef_rgb_path[i]))))
  
  # read in image
  georef_rgb <- rast(georef_rgb_path[i])
  
  # find the corresponding shapefile
  idx_d<-which(gsub("_classes","",tools::file_path_sans_ext(basename(shapefiles)))==gsub("_georef","",tools::file_path_sans_ext(basename(georef_rgb_path[i]))))
  Categorised_shp<-st_read(shapefiles[idx_d])
  Categorised_shp$Class<-as.factor(Categorised_shp$Class)
  
  
  
  ################################################### read in alternative colour models
  # create a hsv image
  if(FALSE %in% as.logical(names(georef_rgb)==c('red','green','blue'))){
    names(georef_rgb)<-c('red','green','blue')
  }
  RGB(georef_rgb)<-c(1,2,3)
  georef_hsv<-colorize(georef_rgb,to='hsv')
  names(georef_hsv)<-c('Hue_sRGB','Saturation_sRGB','Value_sRGB')
  
  ### alternative colour models:
  idx_a<-which(gsub("_PPrgb","",tools::file_path_sans_ext(basename(PPrgb_colourmodel)))==img_details$Code[idx])
  PPrgb_tif <- rast(PPrgb_colourmodel[idx_a])
  PPrgb_tif<-dummy_georef(PPrgb_tif)
  names(PPrgb_tif)<-c('Red_PP','Green_PP','Blue_PP')
  RGB(PPrgb_tif)<-NA
  
  idx_a<-which(gsub("_lab","",tools::file_path_sans_ext(basename(lab_colourmodel)))==img_details$Code[idx])
  lab_tif <- rast(lab_colourmodel[idx_a])
  lab_tif<-dummy_georef(lab_tif)
  names(lab_tif)<-c('L','A','B')
  RGB(lab_tif)<-NA
  lab_tif$Chroma<-sqrt(lab_tif$A^2 + lab_tif$B^2)
  lab_tif$Hue <- (atan2(lab_tif$B, lab_tif$A) * 180 / pi + 360) %% 360
  
  idx_a<-which(gsub("_xyz","",tools::file_path_sans_ext(basename(xyz_colourmodel)))==img_details$Code[idx])
  xyz_tif <- rast(xyz_colourmodel[idx_a])
  xyz_tif<-dummy_georef(xyz_tif)
  names(xyz_tif)<-c('X','Y','Z')
  RGB(xyz_tif)<-NA
  
  
  
  # Code to fix regular naming issues.
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
    Categorised_shp$Class<-droplevels(Categorised_shp$Class)
    st_write(Categorised_shp,shapefiles[idx_d],append=F)
    print("deleted empty vector layer")
  }
  
  # Check there are no additional odd categories.
  if(FALSE %in% as.logical(levels(Categorised_shp$Class) %in% unique_cats)){
    print(paste0('ERROR: odd category in ',basename(shapefiles[idx_d])))
    break
  }
  
  # check there's nothing with no category label
  if(TRUE %in% as.logical(is.na(Categorised_shp$Class))){
    print(paste0('Null category in ',basename(shapefiles[idx_d])))
    break
  }

  # sample locations within each class
  grid = rast(Categorised_shp, nrow=nrow(georef_rgb), ncol=ncol(georef_rgb))
  ncr = terra::rasterize(Categorised_shp, grid, field="Class")
  
  
  # check for levels with too few values
  if(T %in% as.logical(as.data.frame(table(terra::as.data.frame(ncr)$Class))$Freq<npixels)){
    print(paste0('ERROR: not enough values in ',basename(shapefiles[idx_d])," ",
                 levels(terra::as.data.frame(ncr)$Class)[which(as.logical(as.data.frame(table(terra::as.data.frame(ncr)$Class))$Freq<npixels)==T)]))
    break
  }
  
  
  # sample pixels in all categories listed
  ########################################################################################
  samples <- spatSample(ncr, size = npixels, method="stratified",as.points=T,xy=T)


  # use names rather than levels for the sample group (align categories)
  samples$Class<-as.factor(samples$Class)
  levels(samples$Class)<-levels(Categorised_shp$Class)
  
  
  #extract spectral values from each image
  df <- terra::extract(c(georef_rgb,georef_hsv,lab_tif,PPrgb_tif,xyz_tif), samples,ID=F)
  sampdata <- data.frame(Class = samples$Class, df,X=samples$x,Y=samples$y)
  sampdata$Code <- gsub("_georef",tools::file_path_sans_ext(basename(georef_rgb_path[i])),replacement = "")
  pixval_list[[i]]<-sampdata
  #########################################################################################
  
  # create a second set of descriptors for image classification
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
  Categorised_shp$Maj_class<-as.factor(Categorised_shp$Maj_class)
  Categorised_shp<-Categorised_shp[!as.logical(is.na(Categorised_shp$Maj_class)),]
  
  
  ncr = terra::rasterize(Categorised_shp, grid, field="Maj_class")
  samples <- spatSample(ncr, size = npixels, method="stratified",as.points=T) # removed xy=T
  # use names rather than levels for the sample group
  samples$Class<-as.factor(samples$Maj_class)
  levels(samples$Class)<-levels(Categorised_shp$Maj_class)
  
  #extract spectral values from each image
  df <- terra::extract(c(georef_rgb,georef_hsv,lab_tif,PPrgb_tif,xyz_tif), samples, ID=F)
  sampdata <- data.frame(Class = samples$Class, df)
  sampdata$Code <- gsub("_georef",tools::file_path_sans_ext(basename(georef_rgb_path[i])),replacement = "")
  sampdata_list[[i]]<-sampdata
  
  
  #########################################################################################
  print(i)
  
  
  
  
}         ################### main pixel sampling loop
sample_dataframe<-do.call(rbind,sampdata_list)
pixval_list<-do.call(rbind,pixval_list)

# fix class.
sample_dataframe$Class<-as.factor(sample_dataframe$Class)
pixval_list$Class<-as.factor(pixval_list$Class)

# save so we dont have to repeat above.
saveRDS(sample_dataframe,paste0(root_dir,'Pixel_value_training_dataset.rds'))

# save so we dont have to repeat above.
saveRDS(pixval_list, paste0(root_dir,'Pixel_value_for_visualisation.rds'))
