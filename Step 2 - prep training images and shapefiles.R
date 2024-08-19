# step 2 - create a database of files prepped for training

# Training database prep script
# fix integer issue from using python to read in .arw files
# Create georef images and blank shapefiles
# After creating this, open the paired shapefile and brightened RGB in qgis and label each category
# we have used the labels: c('Leaves','Dark wood','Light wood','Sun','Sky') although you should determine what is appropriate (e.g., using a single bark or wood category and determining if a flowering/new shoot/seed/infrastructure category is needed)
# labelling notes: - blank shapefiles have blank numeric 'ID' and a text "Class" columns
#                  - there is 1 row with Class='txt' that should be deleted when labelling





# determine where you want the blank shapefiles to end up
# the best approach is to have two locations
# the first is where you save the blank shapes and georef images
# then the second is where you'll move them after they're labelled
blank_shapes="C:/Users/22064705/R directory/Fortescue hemi img processing/Submission scripts/training_needs_labelling/"
end_path="C:/Users/22064705/R directory/Fortescue hemi img processing/Submission scripts/training_labelled/"



# Load required packages if not already installed
library(terra)
library(sf)

# additional functions
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
fix_python_ints<-function(raster_img,bits=16){
  # 16 bit rasters from python are problematic, as terra::writeRaster uses the 0 for NAs where they are present
  # we therefore need to subtract 1 from saturated values anything that has saturated values.
  if(TRUE %in% as.logical(terra::values(raster_img)==65535)&bits==16){
    terra::values(raster_img)<-round(terra::values(raster_img)*(65534/65535))
  }
  if(TRUE %in% as.logical(terra::values(raster_img)==256)&bits==8){
    terra::values(raster_img)<-round(terra::values(raster_img)*(255/256))
  }
  return(raster_img)
}



##########################################################
###### read in files
##########################################################
path=c("F:/Oct 23 Hemi tiffs LAI paper")
root_dir = "C:/Users/22064705/R directory/Fortescue hemi img processing/Submission scripts/"
rgb_tiffs<-Sys.glob(file.path(paste0(root_dir,'processed_tiffs/'), "*.tif")) # specify whether you want the autobright file HERE by replacing "*.tif" with "*bright.tif"
image_list<-read.csv(paste0(root_dir,"ID.csv")) # this is a list of hemispheric images w/metadata for each file

# to filter images. it's likely you'll only use a subset so choose a couple worth looking at!
# image_list<-image_list[image_list$Position=='C',] # filter based on position


################################################################
################# cleaned loop chunk ###########################
################################################################
#idx<-which(tools::file_path_sans_ext(basename(unstretched))==image_list$Code[which(!image_list$Code%in%gsub("_georef","",tools::file_path_sans_ext(basename(georef_path))))[101]])

for (im_path in 1:nrow(image_list)){    # START
  # select an image. can remove index step if image_list is just raw files
  idx<-which(tools::file_path_sans_ext(basename(rgb_tiffs))==image_list$Code[im_path])
  code<-tools::file_path_sans_ext(basename(rgb_tiffs))[idx]
  raster_img <- rast(rgb_tiffs[[idx]])
  
  # fix python number coding issue
  raster_img<-fix_python_ints(raster_img,bits=16)
  
  # name the layers
  names(raster_img)<-c('red','green','blue')
  
  # check for NAs add break after print if it needs to be removed
  if(TRUE %in% as.logical(is.na(terra::values(raster_img)))){
    print('NAs in raster ')
  }
  
  # add the dummy georeference so you can read it into QGIS
  raster_img<-dummy_georef(raster_img)
  
  # check the file isn't already there and save it.
  # important if you're running this twice for some reason and don't want to delete your training data!
  if(!file.exists(paste0(end_path,code,"_georef.tiff"))){
    writeRaster(raster_img,paste0(blank_shapes,code,"_georef.tiff"),overwrite=T,datatype="INT2U")
  }else{
    print('File already exists!')
  }
  
  # save an empty shapefile with the same dimensions and CRS
  nrows=1
  df <- st_sf(id = 1:nrows,Class="txt", geometry = st_sfc(lapply(1:nrows, function(x) st_multipolygon())))
  st_crs(df)<-'EPSG:3857'
  
  if(!file.exists(paste0(end_path,code,"_classes.shp"))){
    write_sf(df,paste0(blank_shapes,code,"_classes.shp"),append=F)
  }
  print(code)
}

