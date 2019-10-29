
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

strata.limits <- data.frame('STRATA'="All_areas")
Region = "user"
DateFile = paste0('/Users/jiecao/Desktop/UW_work/model_devp/comparison')

#DF = survey_data
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
    #if (sampling_error) {
    #  for (l in 1:(nrow(Spatial_List$Kmeans$centers))){
    #    temp1[l] = rlnorm(1,log(temp1[l]),logcatch_p)
    #  }
    #}
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
save(Save, file=paste0(DateFile,'/',"Save",dorun,".RData"))


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


