
library(TMB)
source('/Users/jiecao/Desktop/UW_work/model_devp/simulator/cal_GM.r')
source('/Users/jiecao/Desktop/UW_work/model_devp/SSST/Codes/Data_Fn2.r')
source('/Users/jiecao/Desktop/UW_work/model_devp/SSST/Codes/Build_TMB_Fn2.r')
Version = c("spatial_size_8")[1]

n_x = c(30, 50, 75, 100, 150, 200, 300)[2] # Number of stations
Kmeans_Config = list( "randomseed"=1, "nstart"=100, "iter.max"=1e3 )

Nfactors_est = 5        # Number of dynamic factors in process error
Use_REML = FALSE
Estimate_Phi = TRUE   # Phi is the offset of initial abundance
Estimate_select = TRUE
#Estimate_recruitment = TRUE
Kappa = c("constant", "spatial_vs_spatiotemporal", "different")[1]  
ObsModel = c("Poisson","Lognormal","ZILN","LNP","Normal")[2]
#sampling_error = FALSE #(sampling error for catch data)
#logcatch_p = 0.25 # generating data

#M_pkt = array(0.23, dim=c(5,66,10))

low_s = 9  # lower boundary of size range
up_s  = 34 # upper boundary of size range
bin_width = 5 
n_p = (up_s-low_s)/bin_width # number of bins
s_mid_points = seq(low_s+bin_width/2, up_s-bin_width/2, bin_width)
size_midpoints = s_mid_points

n_r = 1
R_size = c(1,0,rep(0,n_p-2))
mature = rep(0.5,n_p)

#growth_par = c(100,0.1,0.5,0.01,0.9)
growthmale_tran = cal_GM(s_mid_points[1],s_mid_points[n_p],bin_width,36,1.0,0.3,0.2,0.9)
#growthmale_tran = cal_GM(s_mid_points[1],s_mid_points[n_p],bin_width,36,1,0.3,0.2,0.9)

strata.limits <- data.frame('STRATA'="All_areas")
Region = "user"

if(determ){
  DateFile = paste0('/Users/jiecao/Desktop/UW_work/model_devp/SSST','/SSST_output/',format(Sys.time(), "%d-%b-%Y %H.%M"))
}else{
  DateFile = paste0('/Users/jiecao/Desktop/UW_work/model_devp/SSST','/SSST_output/','iterations',n_loc_sample)
}

#DateFile = paste0('/Users/jiecao/Desktop/UW_work/model_devp/SSST','/SSST_output/',format(Sys.time(), "%d-%b-%Y %H.%M"))
dir.create(DateFile)

DF = read.csv('/Users/jiecao/Desktop/UW_work/model_devp/SSST/Simulation_Data/shrimp/survey_data.csv')
Data_Geostat = cbind( "Size_class"=DF[,"size_class"], "Year"=DF[,"year"], "Catch_N"=DF[,"count"], "AreaSwept_km2"=DF[,"area_swept"], "Vessel"=0, "Lat"=DF[,"lat"], "Lon"=DF[,"lon"] )

pander::pandoc.table( head(Data_Geostat), digits=3 )
pander::pandoc.table( tail(Data_Geostat), digits=3 )

load('/Users/jiecao/Google Drive/UW_work/snow_crab/simulator/shrimp_example/shrimp_grid.rda')
Extrapolation_List = SpatialDeltaGLMM::Prepare_Extrapolation_Data_Fn( Region=Region, strata.limits=strata.limits, input_grid = input_grid )

Spatial_List = SpatialDeltaGLMM::Spatial_Information_Fn( grid_size_km=50, n_x=n_x, Method="Mesh", Lon=Data_Geostat[,'Lon'], Lat=Data_Geostat[,'Lat'], Extrapolation_List=Extrapolation_List, randomseed=Kmeans_Config[["randomseed"]], nstart=Kmeans_Config[["nstart"]], iter.max=Kmeans_Config[["iter.max"]], DirPath=DateFile )
# Add knots to Data_Geostat
Data_Geostat = cbind( Data_Geostat, "knot_i"=Spatial_List$knot_i )

# aggregate catch
catch_raw = read.csv('/Users/jiecao/Desktop/UW_work/model_devp/SSST/Simulation_Data/shrimp/catch_data.csv')
#pander::pandoc.table( head(catch_raw), digits=3 )
#pander::pandoc.table( tail(catch_raw), digits=3 )

n_year = length(unique(Data_Geostat[,'Year']))

loc_catch = SpatialDeltaGLMM::Convert_LL_to_UTM_Fn(Lon = catch_raw$lon, Lat = catch_raw$lat, zone = Extrapolation_List$zone, flip_around_dateline = Extrapolation_List$flip_around_dateline)
NN_catch  = RANN::nn2(data = Spatial_List$Kmeans$centers[, c("E_km", "N_km")], query = loc_catch[, c("X", "Y")], k = 1)
catch_raw = cbind(catch_raw,NN_catch$nn.idx)
c_pkt = array(0, dim=c(n_p,Spatial_List$MeshList$anisotropic_mesh$n,n_year) )

n_boundary = Spatial_List$MeshList$anisotropic_mesh$n-nrow(Spatial_List$Kmeans$centers)

for (i in 1:n_year){
  temp = subset(catch_raw,year ==i)
    for (k in 1:n_p){
      temp1 = tapply(temp[, (4+k)], INDEX = factor(temp[,'NN_catch$nn.idx'], levels = 1:nrow(Spatial_List$Kmeans$centers)), FUN = sum)
      if (sampling_error) {
        for (l in 1:(nrow(Spatial_List$Kmeans$centers))){
          temp1[l] = rlnorm(1,log(temp1[l]),logcatch_p)
        }
      }
      c_pkt[k,,i] = c(temp1,rep(0,n_boundary))
    }
}

c_pkt[is.na(c_pkt)] <- 0

#catchlogsd = 0.01 # used in EM

#c_pkt = NULL

# Build and run model

# Compile TMB software

#compile( '/Users/jiecao/Desktop/SSST/SSST_output/spatial_size_1.cpp' )
#compile(paste('/Users/jiecao/Desktop/SSST/SSST_output/',Version,'.cpp',sep=''))
TmbDir = c('/Users/jiecao/Desktop/UW_work/model_devp/SSST/Codes')
file.copy( from=paste0(TmbDir,"/",Version,".cpp"), to=paste0(DateFile,"/",Version,".cpp"), overwrite=FALSE)
## Build model

ObsModel_p = rep( switch(ObsModel,"Poisson"=0, "Lognormal"=1, "ZILN"=2, "LNP"=3), length(unique(Data_Geostat[,'Size_class'])) )
TmbData = Data_Fn( "Version"=Version, "obsmodel_p"=ObsModel_p, "b_i"=Data_Geostat[,'Catch_N'], "c_pkt"=c_pkt, "s_i"=Data_Geostat[,'knot_i'], "t_i"=Data_Geostat[,'Year'], "p_i"=Data_Geostat[,'Size_class'], "a_x"=Spatial_List$a_xl[,1], "n_r"=n_r, "R_sex"=R_sex, "R_size"=R_size, "mature"=mature, "size_midpoints"=size_midpoints,"M_pkt"=NULL, "growthmale_tran"=growthmale_tran, "MeshList"=Spatial_List$MeshList, "n_factors"=Nfactors_est, catchlogsd=catchlogsd )
TmbList = Build_TMB_Fn( "TmbData"=TmbData, "Version"=Version, "use_REML"=ifelse(is.na(Use_REML),TRUE,Use_REML), "loc_x"=Spatial_List$MeshList$loc_x, "estimate_phi"=Estimate_Phi, "estimate_select"=Estimate_select,"Kappa"=Kappa, "eigenbounds"=EigenBounds, "RunDir"=DateFile )
obj = TmbList$Obj
#obj_normalized = TMB::normalize(obj, flag="include_data", value=FALSE)
runSymbolicAnalysis(obj)
Opt = TMBhelper::Optimize( obj=obj, lower=TmbList$Lower, upper=TmbList$Upper, getsd=TRUE, savedir=DateFile, newtonsteps=1 )

Report = obj$report()
ParHat = obj$env$parList()

# Save everything in object "Save" for it to be reloaded
Save = list("Opt"=Opt, "Report"=Report, "ParHat"=ParHat, "TmbData"=TmbData, "Data_Geostat"=Data_Geostat)
save(Save, file=paste0(DateFile,'/',"Save",iter,".RData"))

require(plotrix)
png(paste('Total_abunance_N',n_loc_sample,'K',n_x,'.png',sep=''), height = 3, width = 12, units = 'in', res=600)


par(mfrow=c(1,5))
for (k in 1:n_p){
  min = min(index[,k]); max = max(index[,k])
  plot(1:n_t,index[,k],type='l',ylim = c(0.8*min,1.2*max),ylab='Total Abundance',xlab='Year',main=paste('Size class','',k),col='red')
  U = Report$Index_tp[,k] + 1.96*Opt$SD$sd[(n_t*(k-1)+1):(n_t*(k-1)+n_t)]
  L = Report$Index_tp[,k] - 1.96*Opt$SD$sd[(n_t*(k-1)+1):(n_t*(k-1)+n_t)]
  plotCI(1:n_t, Report$Index_tp[,k], ui=U, li=L,add=TRUE,pch=16)
  if (k==n_p){
    legend("topright", legend = c("True", "Estimated"), col = c('red','black'), lty = c(1,NA), pch=c(NA,16), bty = "n", text.col = "black", horiz = F , inset = c(0.1, 0.1))
  }
}
dev.off()

# catch plots
catch_est = Save$Report$chat_pkt
catch_ob = Save$TmbData$c_pkt

n_t = Save$TmbData$n_t; n_p = Save$TmbData$n_p
catch_pt_ob = matrix(NA,nrow=Save$TmbData$n_t,ncol=Save$TmbData$n_p)
for (t in 1:n_t){
  for (p in 1:n_p){
    catch_pt_ob[t,p] = sum(catch_ob[p,,t])
  }
}
catch_pt_est = Save$Report$Catch_tp

require(plotrix)
png(paste('Total_Catch_fit','-K',n_x,'.png',sep=''), height = 3, width = 12, units = 'in', res=600)
par(mfrow=c(1,5))
for (k in 1:n_p){
  ind = 2*n_t*n_p
  min = min(catch_pt_ob[,k]); max = max(catch_pt_ob[,k])
  plot(1:n_t,catch_pt_ob[,k],type='l',ylim = c(0.8*min,1.2*max),ylab='Total Catch',xlab='Year',main=paste('Size class','',k),col='red')
  U = Save$Report$Catch_tp[,k] + 1.96*Save$Opt$SD$sd[(n_t*(k-1)+ind+1):(n_t*(k-1)+ind+n_t)]
  L = Save$Report$Catch_tp[,k] - 1.96*Save$Opt$SD$sd[(n_t*(k-1)+ind+1):(n_t*(k-1)+ind+n_t)]
  plotCI(1:n_t, Save$Report$Catch_tp[,k], ui=U, li=L,add=TRUE,pch=16)
  if (k==n_p){
    legend("topright", legend = c("Observed", "Estimated"), col = c('red','black'), lty = c(1,NA), pch=c(NA,16), bty = "n", text.col = "black", horiz = F , inset = c(0.1, 0.1))
  }
}
dev.off()


# plot data

if (plot_figures == TRUE)
{
SpatialDeltaGLMM::Plot_data_and_knots(Extrapolation_List=Extrapolation_List, Spatial_List=Spatial_List, Data_Geostat=Data_Geostat, PlotDir=paste0(DateFile,'/') )
MapDetails_List = SpatialDeltaGLMM::MapDetails_Fn( "Region"=Region, "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap, "Extrapolation_List"=Extrapolation_List )
# 
# Mat_xt = Report$d_pkt[5, ,]
# category_names = 1:dim(Report$d_pkt)[1]
# 
# Year_Set = 1:dim(Report$d_pkt)[3]
# Years2Include = 1:dim(Report$d_pkt)[3]
# Ncategories = dim(Report$d_pkt)[1]
# FileName = paste0(getwd(), "/")
# mfrow = c(ceiling(sqrt(length(Years2Include))), ceiling(length(Years2Include)/ceiling(sqrt(length(Years2Include)))))
# 
# PlotMap_Fn(MappingDetails = MapDetails_List[["MappingDetails"]], 
#            Mat = Mat_xt[, Years2Include, drop = FALSE], 
#            PlotDF = MapDetails_List[["PlotDF"]], MapSizeRatio = MapDetails_List[["MapSizeRatio"]], 
#            Xlim = MapDetails_List[["Xlim"]], Ylim = MapDetails_List[["Ylim"]], FileName = paste0(FileName, 
#                                              ifelse(Ncategories > 1, paste0("size_class_", category_names[5]), "")), 
#            Year_Set = Year_Set[Years2Include], Rescale = FALSE, 
#            Rotate = MapDetails_List[["Rotate"]], Format = "png", Res = 200, 
#            zone = MapDetails_List[["Zone"]], Cex = MapDetails_List[["Cex"]], textmargin = NULL, 
#            add = FALSE, pch = NULL, Legend = MapDetails_List[["Legend"]], mfrow = mfrow, 
#            plot_legend_fig = FALSE)

# comparison plots

# catch on the map
library(raster); library(maptools)
ocean<-readShapeSpatial("/Users/jiecao/Desktop/UW_work/model_devp/simulator/50m_land.shp")
PlotDF1 = subset(MapDetails_List$PlotDF,Include==TRUE)
coords2 = cbind(PlotDF1$Lon,PlotDF1$Lat)
coords2 = coords2[coords2[,1]<0,]
P2 = SpatialPoints(coords2)
rast <- raster(ncol=80,nrow=80)
extent(rast) <- extent(coords2)
xlim = c(-72,-67) ; ylim = c(40,46)
col=colorRampPalette(c("blue", "cyan", "yellow", "red"), bias=1)
#breakpoints <- c(0,250,500,750,1000,1250,1400,1500,1600,1750,2000,3000,5000,10000)

for (t in 1:n_t){
  png(paste('catch','_year-',t,'.png',sep=''), height = 8, width = 6, units = 'in', res=600)
  par(mfrow=c(n_p,2))
  par(mar = c(0, 0, 0, 0), oma = c(4, 4, 4, 0.5))
  par(tcl = -0.25)
  par(mgp = c(2, 0.6, 0))
  for (p in 1:n_p){
    P2$data = catch_ob[p,,t][PlotDF1$x2i]
    rast.temp <- rasterize(P2, rast, P2$data, fun = mean)
    plot(ocean,col="dark gray",axes=F,xlim=xlim,ylim=ylim)
    image(rast.temp,col=col(13),axes=TRUE,add=T,xlim=c(-180,-150),ylim=c(55,65),zlim = c(min,max))
    mtext(paste('size class',p,sep=''),1,-2,adj=0.95,cex=0.8)
    mtext(paste('Observed',sep=''),1,-2,adj=0.1,cex=0.8)
    
    P2$data = catch_est[p,,t][PlotDF1$x2i]
    rast.temp <- rasterize(P2, rast, P2$data, fun = mean)
    plot(ocean,col="dark gray",axes=F,xlim=xlim,ylim=ylim)
    image(rast.temp,col=col(13),axes=TRUE,add=T,xlim=c(-180,-150),ylim=c(55,65),zlim = c(min,max))
    mtext(paste('size class',p,sep=''),1,-2,adj=0.95,cex=0.8)
    mtext(paste('Estimated',sep=''),1,-2,adj=0.1,cex=0.8)
  }
  dev.off()
}

# density on the map
name = 'den_est' 
den.est = Save$Report$d_pkt

PlotDF1 = subset(MapDetails_List$PlotDF,Include==TRUE)
coords2 = cbind(PlotDF1$Lon,PlotDF1$Lat)
coords2 = coords2[coords2[,1]<0,]
P2 = SpatialPoints(coords2)
rast <- raster(ncol=80,nrow=80)
extent(rast) <- extent(coords2)
col=colorRampPalette(c("blue", "cyan", "yellow", "red"), bias=1)
#breakpoints <- c(0,0.1,0.5,1,1.5,2,3,5,10,50,150)
breakpoints <- c(0,0.3,0.5,0.7,0.9,1.1,1.3,1.5,2,2.5,3,3.5,4,5)

png(paste(name,'-N',n_loc_sample,'K',n_x,'.png',sep=''), height = 9, width = 18, units = 'in', res=600)
par(mfrow=c(n_p,n_t))
par(mar = c(0, 0, 0, 0), oma = c(4, 4, 4, 0.5))
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))
for(p in 1:n_p){
  for(t in 1:n_t){
    P2$data = exp(den.est[p,,t][PlotDF1$x2i])
    rast.temp <- rasterize(P2, rast, P2$data, fun = mean)
    plot(ocean,col="dark gray",axes=F,xlim=xlim,ylim=ylim)
    image(rast.temp,col=col(13),axes=TRUE,add=T,xlim=xlim,ylim=ylim,zlim = c(min,max))
    mtext(paste('Size class',p,',','Year',t,sep=''),1,-2,adj=1,cex=0.5)
  }
}
dev.off()

for(p in 1:n_p){
  png(paste(name,'_size-',p,'.png',sep=''), height = 8, width = 6, units = 'in', res=600)
  par(mfrow=c(4,3))
  par(mar = c(0, 0, 0, 0), oma = c(4, 4, 4, 0.5))
  par(tcl = -0.25)
  par(mgp = c(2, 0.6, 0))
  
  for(t in 1:n_t){
    P2$data = exp(den.est[p,,t][PlotDF1$x2i]) 
    rast.temp <- rasterize(P2, rast, P2$data, fun = mean)
    plot(ocean,col="dark gray",axes=F,xlim=xlim,ylim=ylim)
    image(rast.temp,col=col(13),axes=TRUE,add=T,xlim=xlim,ylim=ylim,zlim = c(min,max))
    mtext(paste('Year',t,sep=''),1,-2,adj=0.95,cex=0.8)
  }
  dev.off()
}

# den.est.exp = array(data = 0, dim = c(n_s,n_p,n_t))
# for(p in 1:n_p){
#   for(t in 1:n_t){
#     den.est.exp[,p,t]=exp(den.est[p,,t][PlotDF1$x2i])
#   }
# }
# 
# break_tot <- c(1,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,15)
# Plot_annimation(P,den.est.exp,n_t,Dir,name,breaks=break_tot)

# fishing mortality

F_est = exp(Save$Report$logF_male_kt)
breakpoints=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.5,2)

png(paste('F_est','.png',sep=''), height = 4, width = 8, units = 'in', res=600)
par(mfrow=c(2,5))
par(mar = c(0, 0, 0, 0), oma = c(1, 1, 1, 0.5))
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))

for (t in 1:n_t){
    P2$data = F_est[,t][PlotDF1$x2i]
    rast.temp <- rasterize(P2, rast, P2$data, fun = mean)
    plot(ocean,col="dark gray",axes=F,xlim=xlim,ylim=ylim)
    image(rast.temp,col=col(14),axes=TRUE,add=T,xlim=c(-180,-150),ylim=c(55,65))
    mtext(paste('year',t,sep=''),1,-2,adj=0.95,cex=0.8)
  }
  dev.off()

  # animation
  for (i in 1:n_t){
    P2$data = F_est[,i][PlotDF1$x2i]
    col=colorRampPalette(c("blue", "cyan", "yellow", "red"), bias=1)
    png(paste('F_est','_year',i,'.png',sep=''), height = 6, width = 6, units = 'in', res=900)
    par(mar = c(0, 0, 0, 0), oma = c(1, 1, 1, 0.5))
    par(tcl = -0.25)
    par(mgp = c(2, 0.6, 0))
    rast.temp <- rasterize(P2, rast, P2$data, fun = mean)
    plot(ocean,col="dark gray",axes=F,xlim=xlim,ylim=ylim)
    image(rast.temp,col=col(14),axes=TRUE,add=T,xlim=xlim,ylim=ylim)
    mtext(paste('Year',i,sep=''),1,-2,adj=0.95,cex=0.8)
    dev.off()
    img_temp = image_read(paste('F_est','_year',i,'.png',sep=''))
    if (i==1) img = img_temp
    img = c(img,img_temp)
  }
  img <- image_scale(img, "900")
  animation = image_animate(img, fps = 0.5, dispose = "previous")
  image_write(animation, "F_est.gif")

# covariance among size classess
  pander::pandoc.table( Save$Report$Cov_pp, digits=3 )
 
  png("corplot.png", width=5, height=5, units="in", res=600)
  op <- par(mar=c(6,6,1,1), ps=10)
  COR <- Save$Report$Cov_pp
  name_size = c('size 1','size 2','size 3','size 4','size 5')
  row.names(COR) = name_size
  colnames(COR) = name_size
  
  #image(COR,  col = cols2)
  image(x=seq(nrow(COR)), y=seq(ncol(COR)), z=COR, axes=F, xlab="", ylab="")
  text(expand.grid(x=seq(dim(COR)[1]), y=seq(dim(COR)[2])), labels=round(c(COR),3))
  box()
  axis(1, at=seq(nrow(COR)), labels = rownames(COR), las=2)
  axis(2, at=seq(ncol(COR)), labels = colnames(COR), las=1)
  par(op)
  dev.off()

}





