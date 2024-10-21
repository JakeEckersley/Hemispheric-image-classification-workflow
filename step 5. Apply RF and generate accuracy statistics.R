# imagewise RF model training and application
### Libraries + functions #####################################################################
library(randomForest)
library(terra)
library(scales)
library(dismo)
library(gtools)
library(sf)
library("segmented")

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
brighten<-function(img.rgb,brightness=0,gamma=1.8){ #brightness is on a scale of 0 to 100%
  img.values <- terra::values(img.rgb)
  if (is.numeric(gamma) & gamma !=1){
    minm<-min(img.values,na.rm=TRUE)
    maxm<-max(img.values,na.rm=TRUE)
    img.values=(maxm-minm)*((img.values/(maxm-minm))^gamma)
  }
  terra::values(img.rgb)<-img.values
  RGB(img.rgb)<-c(1,2,3)
  if(is.numeric(brightness)&brightness>0){
    img.hsv<-colorize(img.rgb,to='hsv')
    img.values <- terra::values(img.hsv)
    minm<-min(img.values[,3],na.rm=TRUE)
    maxm<-max(img.values[,3],na.rm=TRUE)
    adjustment<-maxm/100*brightness
    img.values[,3]<-img.values[,3]+adjustment
    img.values[,3][img.values[,3]>maxm]<-maxm
    terra::values(img.hsv)<-img.values
    RGB(img.hsv,type='hsv')<-c(1,2,3)
    img.rgb<-colorize(img.hsv,to='rgb')
  }
  return(img.rgb)
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



# calculates effective L P and W area index and clumping corrections
# annoyingly i've used 'Theta' for the fkn azimuth here and for the zenith in the paper. RIP.
getLAI <- function(data_frame,
                   sat_LAI=5,sat_WAI=5,sat_PAI=8, # user defined inputs for saturated cells
                   group_column='Code'){   # grouping. can be either an image or site code..
  # split the DF by code
  class_by_seg_lst<-split(data_frame,f=data_frame[,group_column])
  LAI_by_img<-list()
  for (i in 1:length(class_by_seg_lst)){
    fcs<-class_by_seg_lst[[i]]
    fcs$Code<-droplevels(as.factor(fcs$Code))
    Code<-unique(fcs$Code)
    LAI_by_ring<-list()
    # calculate weights and add a gap to saturated pixels based on saturated PAI vals
    azimuth<-(as.numeric(unique(fcs$ring))* pi) / 180
    azimuth<-azimuth[order(azimuth)]
    wi<-list()
    max_Lgap<-list()
    max_Wgap<-list()
    max_Pgap<-list()
    for(rngs in 1:length(azimuth)){
      wi[[rngs]]<-sin(azimuth[rngs])/sum(sin(azimuth))
      max_Lgap[[rngs]] <- exp(-sat_LAI / (2 * cos(azimuth[rngs])))
      max_Wgap[[rngs]] <- exp(-sat_WAI / (2 * cos(azimuth[rngs])))
      max_Pgap[[rngs]] <- exp(-sat_PAI / (2 * cos(azimuth[rngs])))
    }
    wts<-do.call(rbind,wi)
    # get gap fractions
    zen_group<-as.numeric(unique(as.character(fcs$ring)))
    zen_group<-zen_group[order(zen_group)]
    for (rngs in 1:length(unique(fcs$ring))){
      idx_a<-which(fcs$ring==zen_group[rngs])
      gap_by_azm<-list()
      fcs_2<-fcs[idx_a,]
      
      
      
      # calculate Pgap; 
      # also calculate Pgap for LX clumping that gives numbers to saturated azum segs
      satL_azm<-list()
      satW_azm<-list()
      satP_azm<-list()
      for (azm in 1:length(fcs_2$azimuth_group)){
        azumith_seg<-fcs_2[azm,]
        
        
        #Leaf gap
        PgapL<-azumith_seg$Sky/(azumith_seg$Sky+azumith_seg$Leaves)
        PgapL_LX<-PgapL
        if(azumith_seg$Sky==0&azumith_seg$Leaves>0) {
          PgapL_LX<-max_Lgap[[rngs]]
          satL_azm[[azm]]<-1
        }else{satL_azm[[azm]]<-0}
        if(!is.finite(PgapL)&(azumith_seg$Sky+azumith_seg$Leaves)==0){ 
          PgapL<-NA
        }
        
        #Wood gap
        PgapW<-azumith_seg$Sky/(azumith_seg$Sky+azumith_seg$Stems)
        PgapW_LX<-PgapW
        if(azumith_seg$Sky==0&azumith_seg$Stems>0){
          PgapW_LX<-max_Wgap[[rngs]]  
          satW_azm[[azm]]<-1
        }else{satW_azm[[azm]]<-0}
        if(!is.finite(PgapW)&(azumith_seg$Sky+azumith_seg$Stems)==0){ 
          PgapW<-NA
        }
        
        PgapT<-azumith_seg$Sky/(azumith_seg$Sky+azumith_seg$Leaves+azumith_seg$Stems)
        PgapT_LX<-PgapT
        if(azumith_seg$Sky==0&azumith_seg$Stems+azumith_seg$Leaves>0){ 
          PgapT_LX<-max_Pgap[[rngs]]
          satP_azm[[azm]]<-1
        }else{satP_azm[[azm]]<-0}
        
        gap_by_azm[[azm]]<-as.data.frame(cbind(PgapL,PgapW,PgapT,PgapL_LX,PgapW_LX,PgapT_LX))
      }
      gap_by_azm<-do.call(rbind,gap_by_azm)
      satL_azm<-sum(do.call(rbind,satL_azm))
      satW_azm<-sum(do.call(rbind,satW_azm))
      satP_azm<-sum(do.call(rbind,satP_azm))
      
      if(satL_azm>0){gap_by_azm$PgapL[gap_by_azm$PgapL==0]<-max_Lgap[[rngs]]}
      if(satW_azm>0){gap_by_azm$PgapW[gap_by_azm$PgapW==0]<-max_Wgap[[rngs]]}
      if(satP_azm>0){gap_by_azm$PgapT[gap_by_azm$PgapT==0]<-max_Pgap[[rngs]]}
      
      
      #effective LAI, PAI and WAI calcs
      LAIe_rng<-(-log(mean(na.omit(gap_by_azm$PgapL))))*cos(azimuth[rngs])*wts[rngs]
      PAIe_rng<-(-log(mean(na.omit(gap_by_azm$PgapT))))*cos(azimuth[rngs])*wts[rngs]
      WAIe_rng<-(-log(mean(na.omit(gap_by_azm$PgapW))))*cos(azimuth[rngs])*wts[rngs]
      
      #if any of them are zeros the above code will shit itself so work this out values
      if(!is.finite(LAIe_rng)&mean(na.omit(gap_by_azm$PgapL))==0){LAIe_rng<-0}
      if(!is.finite(PAIe_rng)&mean(na.omit(gap_by_azm$PgapT))==0){PAIe_rng<-0}
      if(!is.finite(WAIe_rng)&mean(na.omit(gap_by_azm$PgapW))==0){WAIe_rng<-0}
      
      #LX clumping adjusted LAI calcs
      LAI_rng<-(-mean(log(na.omit(gap_by_azm$PgapL_LX))))*cos(azimuth[rngs])*wts[rngs]
      WAI_rng<-(-mean(log(na.omit(gap_by_azm$PgapW_LX))))*cos(azimuth[rngs])*wts[rngs]
      PAI_rng<-(-mean(log(na.omit(gap_by_azm$PgapT_LX))))*cos(azimuth[rngs])*wts[rngs]  
      
      
      
      # LXG clumping measure
      focus_list<-list(gap_by_azm$PgapL,gap_by_azm$PgapW,gap_by_azm$PgapT)
      LXG1_clmp<-list()
      LXG2_clmp<-list()
      for (gapF in 1:length(focus_list)){
        #calculate weights (removing NA segments)
        n_seg<-length(na.omit(focus_list[[gapF]]))
        ORwt1<-(2*(n_seg:1)/(n_seg*(n_seg+1)))
        rcumsum <- function(x) rev(cumsum(rev(x)))
        ORwt2<-rcumsum((1/1:n_seg)/n_seg)
        
        
        asc<-sort(focus_list[[gapF]],decreasing = F)
        desc<-sort(focus_list[[gapF]],decreasing = T)
        LXG1_clmp[[gapF]]<-((-log(sum(desc*ORwt1)))*cos(azimuth[rngs])*wts[rngs])/((-log(sum(asc*ORwt1)))*cos(azimuth[rngs])*wts[rngs])
        LXG2_clmp[[gapF]]<-((-log(sum(desc*ORwt2)))*cos(azimuth[rngs])*wts[rngs])/((-log(sum(asc*ORwt2)))*cos(azimuth[rngs])*wts[rngs])
      }
      
      
      #LXG clumping adjusted LAI  calcs
      LAI_LXG1_rng<-LAIe_rng*LXG1_clmp[[1]]
      WAI_LXG1_rng<-WAIe_rng*LXG1_clmp[[2]]
      PAI_LXG1_rng<-PAIe_rng*LXG1_clmp[[3]]
      
      LAI_LXG2_rng<-LAIe_rng*LXG2_clmp[[1]]
      WAI_LXG2_rng<-WAIe_rng*LXG2_clmp[[2]]
      PAI_LXG2_rng<-PAIe_rng*LXG2_clmp[[3]]
      
      # edit saturated rings
      if(mean(na.omit(gap_by_azm$PgapL))==1){
        LAI_LXG1_rng<-(-log(max_Lgap[[rngs]]))*cos(azimuth[rngs])*wts[rngs]
        LAI_LXG2_rng<-(-log(max_Lgap[[rngs]]))*cos(azimuth[rngs])*wts[rngs]
        Lsat_ring<-1
      }else{Lsat_ring<-0}
      
      if(mean(na.omit(gap_by_azm$PgapW))==1){
        WAI_LXG1_rng<-(-log(max_Wgap[[rngs]]))*cos(azimuth[rngs])*wts[rngs]
        WAI_LXG2_rng<-(-log(max_Wgap[[rngs]]))*cos(azimuth[rngs])*wts[rngs]
        Wsat_ring<-1
      }else{Wsat_ring<-0}
      
      if(mean(na.omit(gap_by_azm$PgapT))==1){
        PAI_LXG1_rng<-(-log(max_Pgap[[rngs]]))*cos(azimuth[rngs])*wts[rngs]
        PAI_LXG2_rng<-(-log(max_Pgap[[rngs]]))*cos(azimuth[rngs])*wts[rngs]
        Psat_ring<-1
      }else{Psat_ring<-0}
      
      
      LAI_by_ring[[rngs]]<-as.data.frame(cbind(LAIe_rng,WAIe_rng,PAIe_rng,
                                               LAI_rng,WAI_rng,PAI_rng,
                                               LAI_LXG1_rng,WAI_LXG1_rng,PAI_LXG1_rng,
                                               LAI_LXG2_rng,WAI_LXG2_rng,PAI_LXG2_rng,
                                               Lsat_ring,Wsat_ring,Psat_ring,
                                               satL_azm,satW_azm,satP_azm
      ))
    }
    LAI_by_ring<-do.call(rbind,LAI_by_ring)
    LAI_by_img_temp<-as.data.frame(t(as.matrix(2*colSums(LAI_by_ring))))
    LAI_by_img_temp$Code<-as.vector(Code)
    LAI_by_img[[i]]<-LAI_by_img_temp
  }
  LAI_by_img<-do.call(rbind,LAI_by_img)
  colnames(LAI_by_img)[1:12]<-c('LAIe','WAIe','PAIe',
                                'LAI_LX','WAI_LX','PAI_LX',
                                'LAI_LXG1','WAI_LXG1','PAI_LXG1',
                                'LAI_LXG2','WAI_LXG2','PAI_LXG2')
  return(LAI_by_img)
}

deWitt_G<-function(){
  deWit<-data.frame(
    mu=c(2.77,1.172,3.326,.433,1,1.101),
    nu=c(1.172,2.77,3.326,.433,1,1.93),
    LAD=c('Planophile','Erectophile','Plagiophile','Extremophile','Uniform','Spherical')
  )
  
  out<-list()
  for(idx in 1:nrow(deWit)){
    mu <- deWit$mu[idx]
    nu <- deWit$nu[idx]
    
    rad=pi/180
    
    realG<-NULL
    zen_angles<-seq(0.5,89.5,by=1)
    
    
    ft=matrix(NA,ncol=length(zen_angles))
    for (i in 1:length(zen_angles)){
      z<-zen_angles[i]
      t=2*z*rad/pi
      a=(1-t)^(mu-1)
      b=t^(nu-1)
      ft[i]=2/pi*(1/beta(mu,nu))*a*b# probability density function
    }
    frel=ft/sum(ft)
    
    Gi <- NULL
    for(i in 1:length(zen_angles)){
      z<-zen_angles[i]
      
      for(j in 1:length(zen_angles)){
        z2<-zen_angles[j]
        (cot_i <- 1/tan(z*rad))
        (cot_fl <- 1/tan(z2*rad))
        
        if(abs(cot_i*cot_fl)>1) {
          A = cos(z*rad)*cos(z2*rad)
        } else {
          A = cos(z*rad)*cos(z2*rad)*(1+(2/pi)*((tan(acos(cot_i*cot_fl)))-acos(cot_i*cot_fl)))
        }
        phi=A*frel[i]
        
        Gi <- c(Gi, phi)
      }
    }
    
    Gmat <- matrix(Gi, ncol=90)
    Gfun<-apply(Gmat, 1, sum)# G-function
    realG<-rbind(realG,round(Gfun,3))
    
    out[[idx]]<-data.frame(zenith=seq(0.5,89.5,by=1),pdf=round(as.numeric(frel),4),G=as.numeric(realG),type=deWit$LAD[idx])
    names(out)[idx]<-deWit$LAD[idx]
  }
  return(out)
}
getLAINils <- function(data_frame,
                       sat_LAI=5,sat_WAI=5,sat_PAI=8, # user defined inputs for saturated cells
                       group_column='Code',G_funct='Extremophile',Z_max=70){   # grouping. can be either an image or site code..
  # split the DF by code
  class_by_seg_lst<-split(data_frame,f=data_frame[,group_column])
  LAI_by_img<-list()
  
  G_list<-deWitt_G()
  G_focus<-G_list[[which(names(G_list)==G_funct)]]
  
  
  for (i in 1:length(class_by_seg_lst)){
    fcs<-class_by_seg_lst[[i]]
    fcs$Code<-droplevels(as.factor(fcs$Code))
    Code<-unique(fcs$Code)
    
    # get the G function
    ring_bins<-unique(fcs$ring)
    G_focus$ring_group<-NA
    for(k in 1:length(G_focus$zenith)){
      G_focus$ring_group[k]<- ring_bins[which((ring_bins - G_focus$zenith[k])^2 == min((ring_bins - G_focus$zenith[k])^2))]
    }
    G_0<-G_focus[G_focus$zenith<Z_max,]
    G_0$zenith_rad<-(G_0$zenith* pi) / 180
    
    LAI_by_ring<-list()
    # calculate weights and add a gap to saturated pixels based on saturated PAI vals
    azimuth<-(as.numeric(unique(ring_bins))* pi) / 180
    azimuth<-azimuth[order(azimuth)]
    
    max_Lgap<-list()
    max_Wgap<-list()
    max_Pgap<-list()
    for(rngs in 1:length(azimuth)){
      max_Lgap[[rngs]] <- exp(sum(-sat_LAI*G_0$G[G_0$ring_group==ring_bins[rngs]])/ sum(cos(G_0$zenith_rad[G_0$ring_group==ring_bins[rngs]])))
      max_Wgap[[rngs]] <- exp(sum(-sat_WAI*0.5) / (cos(azimuth[rngs])))
      max_Pgap[[rngs]] <- exp(sum(-sat_PAI*0.5) / (cos(azimuth[rngs])))
    }
    
    # get gap fractions
    zen_group<-as.numeric(unique(as.character(fcs$ring)))
    zen_group<-zen_group[order(zen_group)]
    for (rngs in 1:length(unique(fcs$ring))){
      idx_a<-which(fcs$ring==zen_group[rngs])
      gap_by_azm<-list()
      fcs_2<-fcs[idx_a,]
      
      
      
      # calculate Pgap; 
      # also calculate Pgap for LX clumping that gives numbers to saturated azum segs
      satL_azm<-list()
      satW_azm<-list()
      satP_azm<-list()
      for (azm in 1:length(fcs_2$azimuth_group)){
        azumith_seg<-fcs_2[azm,]
        
        
        #Leaf gap
        PgapL<-azumith_seg$Sky/(azumith_seg$Sky+azumith_seg$Leaves)
        PgapL_LX<-PgapL
        if(azumith_seg$Sky==0&azumith_seg$Leaves>0) {
          PgapL_LX<-max_Lgap[[rngs]]
          satL_azm[[azm]]<-1
        }else{satL_azm[[azm]]<-0}
        if(!is.finite(PgapL)&(azumith_seg$Sky+azumith_seg$Leaves)==0){ 
          PgapL<-NA
        }
        
        #Wood gap
        PgapW<-azumith_seg$Sky/(azumith_seg$Sky+azumith_seg$Stems)
        PgapW_LX<-PgapW
        if(azumith_seg$Sky==0&azumith_seg$Stems>0){
          PgapW_LX<-max_Wgap[[rngs]]  
          satW_azm[[azm]]<-1
        }else{satW_azm[[azm]]<-0}
        if(!is.finite(PgapW)&(azumith_seg$Sky+azumith_seg$Stems)==0){ 
          PgapW<-NA
        }
        
        PgapT<-azumith_seg$Sky/(azumith_seg$Sky+azumith_seg$Leaves+azumith_seg$Stems)
        PgapT_LX<-PgapT
        if(azumith_seg$Sky==0&azumith_seg$Stems+azumith_seg$Leaves>0){ 
          PgapT_LX<-max_Pgap[[rngs]]
          satP_azm[[azm]]<-1
        }else{satP_azm[[azm]]<-0}
        
        gap_by_azm[[azm]]<-as.data.frame(cbind(PgapL,PgapW,PgapT,PgapL_LX,PgapW_LX,PgapT_LX))
      }
      gap_by_azm<-do.call(rbind,gap_by_azm)
      satL_azm<-sum(do.call(rbind,satL_azm))
      satW_azm<-sum(do.call(rbind,satW_azm))
      satP_azm<-sum(do.call(rbind,satP_azm))
      
      
      if(satL_azm>0){gap_by_azm$PgapL[gap_by_azm$PgapL==0]<-max_Lgap[[rngs]]}
      if(satW_azm>0){gap_by_azm$PgapW[gap_by_azm$PgapW==0]<-max_Wgap[[rngs]]}
      if(satP_azm>0){gap_by_azm$PgapT[gap_by_azm$PgapT==0]<-max_Pgap[[rngs]]}
      
      
      #effective LAI, PAI and WAI calcs
      LAIe_rng<-sum(-log(mean(na.omit(gap_by_azm$PgapL)))*cos(G_0$zenith_rad[G_0$ring_group==ring_bins[rngs]]))/sum(G_0$G[G_0$ring_group==ring_bins[rngs]])
      PAIe_rng<-sum(-log(mean(na.omit(gap_by_azm$PgapT)))*cos(G_0$zenith_rad[G_0$ring_group==ring_bins[rngs]]))/(0.5*nrow(G_0[G_0$ring_group==ring_bins[rngs],]))
      WAIe_rng<-sum(-log(mean(na.omit(gap_by_azm$PgapW)))*cos(G_0$zenith_rad[G_0$ring_group==ring_bins[rngs]]))/(0.5*nrow(G_0[G_0$ring_group==ring_bins[rngs],]))
      
      
      
      #if any of them are 0s the above code will shit itself so this gives it a 0 value for LAI - meaning the LAI is very low...
      if(!is.finite(LAIe_rng)&mean(na.omit(gap_by_azm$PgapL))==0){LAIe_rng<-0}
      if(!is.finite(PAIe_rng)&mean(na.omit(gap_by_azm$PgapT))==0){PAIe_rng<-0}
      if(!is.finite(WAIe_rng)&mean(na.omit(gap_by_azm$PgapW))==0){WAIe_rng<-0}
      
      #LX clumping adjusted LAI calcs
      LAI_rng<-sum(-mean(log(na.omit(gap_by_azm$PgapL_LX)))*cos(G_0$zenith_rad[G_0$ring_group==ring_bins[rngs]]))/sum(G_0$G[G_0$ring_group==ring_bins[rngs]])
      WAI_rng<-sum(mean(-log(na.omit(gap_by_azm$PgapW_LX)))*cos(G_0$zenith_rad[G_0$ring_group==ring_bins[rngs]]))/(0.5*nrow(G_0[G_0$ring_group==ring_bins[rngs],]))
      PAI_rng<-sum(mean(-log(na.omit(gap_by_azm$PgapT_LX)))*cos(G_0$zenith_rad[G_0$ring_group==ring_bins[rngs]]))/(0.5*nrow(G_0[G_0$ring_group==ring_bins[rngs],]))
      
      # edit saturated rings
      if(mean(na.omit(gap_by_azm$PgapL))==1){
        Lsat_ring<-1
      }else{Lsat_ring<-0}
      
      if(mean(na.omit(gap_by_azm$PgapW))==1){
        Wsat_ring<-1
      }else{Wsat_ring<-0}
      
      if(mean(na.omit(gap_by_azm$PgapT))==1){
        Psat_ring<-1
      }else{Psat_ring<-0}
      
      
      
      LAI_LXG1_rng=NA
      WAI_LXG1_rng=NA
      PAI_LXG1_rng=NA
      LAI_LXG2_rng=NA
      WAI_LXG2_rng=NA
      PAI_LXG2_rng=NA
      
      
      LAI_by_ring[[rngs]]<-as.data.frame(cbind(LAIe_rng,WAIe_rng,PAIe_rng,
                                               LAI_rng,WAI_rng,PAI_rng,
                                               LAI_LXG1_rng,WAI_LXG1_rng,PAI_LXG1_rng,
                                               LAI_LXG2_rng,WAI_LXG2_rng,PAI_LXG2_rng,
                                               Lsat_ring,Wsat_ring,Psat_ring,
                                               satL_azm,satW_azm,satP_azm
      ))
    }
    LAI_by_ring<-do.call(rbind,LAI_by_ring)
    LAI_by_img_temp<-as.data.frame(t(as.matrix(colMeans(LAI_by_ring))))
    LAI_by_img_temp$Code<-as.vector(Code)
    LAI_by_img[[i]]<-LAI_by_img_temp
  }
  LAI_by_img<-do.call(rbind,LAI_by_img)
  colnames(LAI_by_img)[1:12]<-c('LAIe','WAIe','PAIe',
                                'LAI_LX','WAI_LX','PAI_LX',
                                'LAI_LXG1','WAI_LXG1','PAI_LXG1',
                                'LAI_LXG2','WAI_LXG2','PAI_LXG2')
  return(LAI_by_img)
}




segment<-function(img=whichmax_img,
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

adj_seg<-function(Seg_img,
                  xc=2049, # x coordinate of mask centre
                  yc=2970, # y coordinate of mask centre
                  rc=1920, # radius of mask
                  nseg = 8,  # Number of azimuth segments
                  startVZA=0, # start viewing zenith angle
                  endVZA=90, # end viewing zenith angle
                  maxVZA=90, # maximum view angle of the camera
                  nrings=45,
                  reduce_outputs=T){ # n-1 for the number of zenith segments){
  imgdf<-Seg_img
  VZAbins<-seq(startVZA,endVZA,(endVZA-startVZA)/nrings)
  maxVZAmm<-2*sin(maxVZA*pi/180)
  rset<-round(rc*2*sin(VZAbins*pi/180)/maxVZAmm)
  
  # Create a new column 'ring' by cutting the 'r' column into segments
  imgdf$ring <- cut(imgdf$r, rset, include.lowest = TRUE,
                    labels = seq(startVZA + ((endVZA - startVZA) / 2 / nrings),
                                 endVZA - ((endVZA - startVZA) / 2 / nrings),
                                 (endVZA - startVZA) / nrings))
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
    idx<-which(colnames(imgdf)%in%c('class','ring','azimuth_group'))
    imgdf<-imgdf[,idx]
  }
  
  return(imgdf)
}



circ_mask<-function(tiff_img,
                    #coordinates for mask (X & Y centre and radius)
                    xc=2049,
                    yc=2970,
                    rc=1920){
  # applymask
  xy <- terra::xyFromCell(tiff_img,1:terra::ncell(tiff_img)) #create xy coordinates
  circular.mask = (xy[,1] - xc)^2 + (xy[,2] - yc)^2 <= rc^2 #determine which pixels are within the circular mask
  terra::values(tiff_img)[!circular.mask] <- NA #everything outside of that is NA
  return(tiff_img)
}
####################################################################################



### file paths #####################################################################

end_path<-"F:/Training tiffs/"
sample_dataframe<-readRDS(paste0(end_path,'Pixel_value_training_dataset.rds'))
img_details<-read.csv(paste0(end_path,"img_details.csv"))
metadata<-read.csv(paste0(end_path,'img_ID_metadata.csv'))[,c(1,9,12)]
all_colourmodels<-c(Sys.glob(file.path("F:/PhotoproRGB_LAB_imgs",'*.tif')),Sys.glob(file.path("F:/PhotoproRGB_LAB_imgs",'*.tiff')))


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
acc_results<-list()
tables<-list()
results<-list()
class_by_seg<-list()
img_a<-list()
class_by_ring<-list()
class_by_seg_4LAI<-list()
img_a_4LAI<-list()
class_by_seg_4LAI2200<-list()
class_by_seg_4LAI45deg<-list()
class_by_seg_57deg<-list()
var_importance<-list()

###################################################################################################
### first half: model training     ################################################################
###################################################################################################
#classify_img<-function(SDF_i,end_path,all_colourmodels,Code,classifier_RF){}
for(i in 1:length(SDF)){
  if(nrow(SDF[[i]])>0){
    #split off each image dataset
    sdf_i<-SDF[[i]]
    sdf_i$Class<-droplevels(sdf_i$Class)
    Code<-names(SDF)[i]
  }else{
    next
  }
  
  # check there is a leaf sky and wood category in each image
  if(!length(levels(sdf_i$Class))==3){
    print(paste0('Level missing in ',sdf_i[1,1]))
    print(levels(sdf_i$Class))
    #total_a[[i]]<-1
    next
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
  
  # specify variables to use in the model
  pixel_vars<-which(colnames(sdf_i)%in%c('L','A','B','Hue','Chroma','red','green','blue','Hue_sRGB','Saturation_sRGB','Value_sRGB','Red_PP','Green_PP','Blue_PP','X','Y','Z'))
  
  #splitting test/train dataset
  train <- sdf_i[j!= k, ]
  test <- sdf_i[j == k, ]
  
  # train the classifier
  classifier_RF = randomForest(x = train[,pixel_vars], 
                               y = train$Class, 
                               ntree = 700,mtry=4,importance = T) 
  
  
  # predict categories of test values
  predicted <-predict(classifier_RF, newdata = test[,pixel_vars])
  
  # combine test/train in a data frame
  x_rf <- cbind(predicted,test)
  colnames(x_rf)[3]<-'observed'
  x_rf$Code<-Code
  results[[i]]<-x_rf
  
  # add overall classification accuracy to a list
  acc_results[[i]]<-cbind(agreement(x_rf),Code)
  
  # add directional misclassification table to a list
  tables[[i]]<-table(x_rf$predicted,x_rf$observed)
  names(tables)[i]<-Code
  
  ## add something that captures 'importance'
  var_importance[[i]]<-as.data.frame(classifier_RF$importance[order(classifier_RF$importance[,4],decreasing=T),drop=F,])
  var_importance[[i]]$Code<-Code
  var_importance[[i]]$colmod<-rownames(var_importance)
  
  # use this if you just want 'tables'
  
  #    print(i)
  #}
  
  ###################################################################################################
  ### second half: model application ################################################################
  ###################################################################################################
  
  ##################################################
  ### Classify the remainder of the image       ####
  ##################################################
  indx<-which(gsub("_lab","",tools::file_path_sans_ext(basename(all_colourmodels)))==Code)
  srgb_indx<-which(gsub("_srgb","",tools::file_path_sans_ext(basename(all_colourmodels)))==Code)
  pprgb_indx<-which(gsub("_PPrgb","",tools::file_path_sans_ext(basename(all_colourmodels)))==Code)
  XYZ_indx<-which(gsub("_xyz","",tools::file_path_sans_ext(basename(all_colourmodels)))==Code)
  
  if(length(indx)>0){# check there is an image to classify....
    print('image found')
    }else{
    print(paste0('error for i = ',i,' - no file in rgb_tiffs named ',Code))
  }
    
  #read in LAB image
  LAB_tiff<-rast(all_colourmodels[[indx]])
  srgb_tiff<-rast(all_colourmodels[[srgb_indx]])
  pprgb_tiff<-rast(all_colourmodels[[pprgb_indx]])
  xyz_tiff<-rast(all_colourmodels[[XYZ_indx]])
  
  RGB(srgb_tiff)<-c(1,2,3)
  
  names(LAB_tiff)<-c("L","A","B")
  names(srgb_tiff)<-c("red","green","blue")
  names(pprgb_tiff)<-c("Red_PP","Green_PP","Blue_PP")
  names(xyz_tiff)<-c("X","Y","Z")
  
  #add Chroma and Hue to image
  LAB_tiff$Chroma<-sqrt(LAB_tiff$A^2 + LAB_tiff$B^2)
  LAB_tiff$Hue <- (atan2(LAB_tiff$B, LAB_tiff$A) * 180 / pi + 360) %% 360
  
  RGB(srgb_tiff)<-c(1,2,3)
  hsv_tiff<-colorize(srgb_tiff,to='hsv')
  names(hsv_tiff)<-c('Hue_sRGB','Saturation_sRGB','Value_sRGB')
  
  combo_tif<-c(LAB_tiff,srgb_tiff,pprgb_tiff,xyz_tiff,hsv_tiff)
  rm(LAB_tiff,pprgb_tiff,srgb_tiff,hsv_tiff,xyz_tiff)
  
  # mask the image
  combo_tif <- circ_mask(combo_tif)

  # get RF classifier votes and maximum likelihood for each pixel 
  #     classified_img <- terra::predict(LAB_tiff, classifier_RF, na.rm = T,type='prob')
  whichmax_img <- terra::predict(combo_tif, classifier_RF, na.rm = T,type='response')
  #      classified_img<-c(classified_img,whichmax_img)
  rm(combo_tif)
  
  # write the maximum likelihood image to file for display
  writeRaster(whichmax_img,paste0(end_path,"classified_imgs/",Code,"_individ_seg_classified.tif"),overwrite=T)
 
  # segment image at 2 degree resolution
  Seg_img<-segment(whichmax_img)
  # quantify leaf to wood ratio in each segment
  Seg_img$azimuth_group<-as.factor(Seg_img$azimuth_group)
  Seg_img$ring<-as.factor(Seg_img$ring)
  rr<-list()
  for (r in 1:length(levels(Seg_img$ring))){
    tt<-list()
    for (t in 1:length(levels(Seg_img$azimuth_group))){
      ring<-levels(Seg_img$ring)[r]
      azimuth_group<-levels(Seg_img$azimuth_group)[t]
      class_prop<-as.data.frame(t(as.matrix(table(Seg_img[Seg_img$ring==ring&Seg_img$azimuth_group==azimuth_group,]$class))))
      a<-class_prop$Stems/(class_prop$Stems+class_prop$Leaves)
      Pgap<-class_prop$Sky/(class_prop$Sky+class_prop$Leaves)
      tt[[t]]<-cbind(ring,azimuth_group,class_prop,a,Pgap,Code)
    }
    rr[[r]]<-do.call(rbind,tt)
  }
  class_by_seg[[i]]<-do.call(rbind,rr)
  class_by_seg[[i]]$ring<-as.factor(class_by_seg[[i]]$ring)
  ringtot<-aggregate(class_by_seg[[i]][,3:5],by=list(class_by_seg[[i]]$ring),FUN=sum)
  total_a<-mean(ringtot$Stems/(ringtot$Stems+ringtot$Leaves),na.rm=T)
  sd_a<-sd(ringtot$Stems/(ringtot$Stems+ringtot$Leaves),na.rm=T)
  img_a[[i]]<-cbind(total_a,sd_a,Code)
  ring_a<-ringtot$Stems/(ringtot$Stems+ringtot$Leaves)
  class_by_ring[[i]]<-cbind(ringtot,ring_a,Code)
  class_by_ring[[i]]$Code<-Code
  
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
      a<-class_prop$Stems/(class_prop$Stems+class_prop$Leaves)
      Pgap<-class_prop$Sky/(class_prop$Sky+class_prop$Leaves)
      tt[[t]]<-cbind(ring,azimuth_group,class_prop,a,Pgap,Code)
    }
    rr[[r]]<-do.call(rbind,tt)
  }
  class_by_seg_4LAI[[i]]<-do.call(rbind,rr)
  class_by_seg_4LAI[[i]]$ring<-as.factor(class_by_seg_4LAI[[i]]$ring)
  total_a<-mean(ringtot$Stems/(ringtot$Stems+ringtot$Leaves),na.rm=T)
  sd_a<-sd(ringtot$Stems/(ringtot$Stems+ringtot$Leaves),na.rm=T)
  img_a_4LAI[[i]]<-cbind(total_a,sd_a,Code)
  
  
  # segment image for LAI calculation  @ 45 degree azm
  seg_img_adj<-adj_seg(Seg_img,
                       nseg = (360/45),  # Number of azimuth segments
                       startVZA=0, # start viewing zenith angle
                       endVZA=70, # end viewing zenith angle
                       maxVZA=90, # maximum view angle of the camera
                       nrings=7,)# for the number of zenith segments
  seg_img_adj$azimuth_group<-as.factor(seg_img_adj$azimuth_group)
  seg_img_adj$ring<-as.factor(seg_img_adj$ring)
  rr<-list()
  for (r in 1:length(levels(seg_img_adj$ring))){
    tt<-list()
    for (t in 1:length(levels(seg_img_adj$azimuth_group))){
      ring<-levels(seg_img_adj$ring)[r]
      azimuth_group<-levels(seg_img_adj$azimuth_group)[t]
      class_prop<-as.data.frame(t(as.matrix(table(seg_img_adj[seg_img_adj$ring==ring&seg_img_adj$azimuth_group==azimuth_group,]$class))))
      a<-class_prop$Stems/(class_prop$Stems+class_prop$Leaves)
      Pgap<-class_prop$Sky/(class_prop$Sky+class_prop$Leaves)
      tt[[t]]<-cbind(ring,azimuth_group,class_prop,a,Pgap,Code)
    }
    rr[[r]]<-do.call(rbind,tt)
  }
  class_by_seg_4LAI45deg[[i]]<-do.call(rbind,rr)
  class_by_seg_4LAI45deg[[i]]$ring<-as.factor(class_by_seg_4LAI45deg[[i]]$ring)
  
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
      a<-class_prop$Stems/(class_prop$Stems+class_prop$Leaves)
      Pgap<-class_prop$Sky/(class_prop$Sky+class_prop$Leaves)
      tt[[t]]<-cbind(ring,azimuth_group,class_prop,a,Pgap,Code)
    }
    rr[[r]]<-do.call(rbind,tt)
  }
  class_by_seg_4LAI2200[[i]]<-do.call(rbind,rr)
  class_by_seg_4LAI2200[[i]]$ring<-as.factor(class_by_seg_4LAI2200[[i]]$ring)
  
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
      a<-class_prop$Stems/(class_prop$Stems+class_prop$Leaves)
      Pgap<-class_prop$Sky/(class_prop$Sky+class_prop$Leaves)
      tt[[t]]<-cbind(ring,azimuth_group,class_prop,a,Pgap,Code)
    }
    rr[[r]]<-do.call(rbind,tt)
  }
  class_by_seg_57deg[[i]]<-do.call(rbind,rr)
  class_by_seg_57deg[[i]]$ring<-as.factor(class_by_seg_57deg[[i]]$ring)
  
  # could turn this into something that operates over multiple cores but will likely run into memory issues.    
  # return(
  #   list(acc_results,
  #         tables,
  #         results,
  #         class_by_seg,
  #         img_a,
  #         class_by_ring,
  #         img_a_4LAI,
  #         class_by_seg_4LAI,
  #         class_by_seg_4LAI2200,
  #         class_by_seg_4LAI45deg,
  #         class_by_seg_57deg))

  print(i)
  rm(Seg_img)
  rm(seg_img_adj)
  
  # for some reason this is georeferencing whichmax_img outside the function? weird.
  writeRaster(dummy_georef(whichmax_img),paste0(end_path,"georef_classified_imgs/",Code,"_individ_seg_classified_georef.tif"),overwrite=T)
  
  rm(whichmax_img)
}








#######
### save outputs
var_importance<-do.call(rbind,var_importance)
var_importance<-merge(var_importance,img_details,by='Code',all.x=T)
aggregate(var_importance$MeanDecreaseAccuracy,by=list(var_importance))
saveRDS(var_importance,paste0(end_path,'colour__model_band_importance.rds'))



indexi<-which(lengths(class_by_seg)>0)
class_by_seg_df<-do.call(smartbind,class_by_seg[indexi])
saveRDS(class_by_seg_df,paste0(end_path,'segmented_45zen_8azm_maxlikelihood_v2.rds'))
indexi<-which(lengths(class_by_ring)>0)
class_by_ring<-do.call(smartbind,class_by_ring[indexi])
saveRDS(class_by_ring,paste0(end_path,'segmented_byzenith__45zen_noazm_maxlikelihood_v2.rds'))
saveRDS(tables,paste0(end_path,'missclass_tables_v2.rds'))
acc_results<-do.call(rbind,acc_results)
saveRDS(acc_results,paste0(end_path,'individual_RF_global_accuracy_results_v2.rds'))
results<-do.call(rbind,results)
saveRDS(results,paste0(end_path,'individual_RF_sample_prediction_accuracy_results_v2.rds'))
img_a<-do.call(rbind,img_a)
saveRDS(img_a,paste0(end_path,'Total_a_by_img_zen0to90_v2.rds'))
img_a_4LAI<-do.call(rbind,img_a_4LAI)
saveRDS(img_a_4LAI,paste0(end_path,'Total_a_by_img_zen0to70_v2.rds'))
class_by_seg_4LAI<-do.call(smartbind,class_by_seg_4LAI[indexi])
saveRDS(class_by_seg_4LAI,paste0(end_path,'segmented_7zen_24azm_maxlikelihood_v2.rds'))
class_by_seg_4LAI2200<-do.call(smartbind,class_by_seg_4LAI2200[indexi])
saveRDS(class_by_seg_4LAI2200,paste0(end_path,'segmented_5zen_8azm_maxlikelihood_v2.rds'))
class_by_seg_4LAI45deg<-do.call(smartbind,class_by_seg_4LAI45deg[indexi])
saveRDS(class_by_seg_4LAI45deg,paste0(end_path,'segmented_7zen_8azm_maxlikelihood_v2.rds'))
class_by_seg_57deg<-do.call(smartbind,class_by_seg_57deg[indexi])
saveRDS(class_by_seg_57deg,paste0(end_path,'segmented_zen57deg_8azm_maxlikelihood_v2.rds'))
###################################################################################################

#test which ones were improved and overall improvement
acc_results_orig<-read.csv(paste0(end_path,'/individual_RF_global_accuracy_results.csv'))
acc_results_orig[which(acc_results_orig$accuracy<.9),]
acc_results_2<-do.call(rbind,acc_results)

for (i in 1:nrow(acc_results_orig)){
  indexio<-which(acc_results_2$Code==acc_results_orig$Code[i])
  if(length(indexio)>0){
    if(acc_results_orig$accuracy[i]<.9&signif(acc_results_2$accuracy,2)[indexio]>=0.9){
      print(acc_results_2$Code[indexio])
    }
  }
}
acc_results_orig<-merge(acc_results_orig,img_details,'Code',all.x=T)
acc_results_2<-merge(acc_results_2,img_details,'Code',all.x=T)

t.test(acc_results_orig$accuracy,acc_results_2$accuracy)

t.test(acc_results_orig$accuracy[acc_results_orig$Sky=='clear'],acc_results_2$accuracy[acc_results_2$Sky=='clear'])
t.test(acc_results_orig$accuracy[acc_results_orig$Sky=='overcast'],acc_results_2$accuracy[acc_results_2$Sky=='overcast'])
t.test(acc_results_orig$accuracy[acc_results_orig$Sky=='partly cloudy'],acc_results_2$accuracy[acc_results_2$Sky=='partly cloudy'])
