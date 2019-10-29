
source('/Users/jiecao/Desktop/UW_work/model_devp/simulator/cal_GM.r')
source('/Users/jiecao/Desktop/UW_work/model_devp/SSST/Codes/Data_Fn2.r')
source('/Users/jiecao/Desktop/UW_work/model_devp/SSST/Codes/Build_TMB_Fn2.r')
source('/Users/jiecao/Desktop/UW_work/model_devp/SSST/Codes/dist_matrix_nearest_neighbor.r')
Version = c("spatial_size_8","spatial_size_9")[1]

n_x = c(30, 50, 75, 100, 150, 200, 300)[2] # Number of stations
Kmeans_Config = list( "randomseed"=1, "nstart"=100, "iter.max"=1e3 )

Nfactors_est = 5        # Number of dynamic factors in process error
Use_REML = FALSE
Estimate_Phi = TRUE   # Phi is the offset of initial abundance
Estimate_select = TRUE
#Estimate_recruitment = TRUE
Kappa = c("constant", "spatial_vs_spatiotemporal", "different")[1]  
ObsModel = c("Poisson", "LNP", "ZILN", "Lognormal")[3]

low_s = 9  # lower boundary of size range
up_s  = 34 # upper boundary of size range
bin_width = 5 
n_p = (up_s-low_s)/bin_width # number of bins
s_mid_points = seq(low_s+bin_width/2, up_s-bin_width/2, bin_width)
size_midpoints = s_mid_points

n_r = 1
R_size = c(1,0,rep(0,n_p-2))
mature = rep(0.5,n_p)

growthmale_tran = cal_GM(s_mid_points[1],s_mid_points[n_p],bin_width,35,1.0,0.2,0.1,0.9)

strata.limits <- data.frame('STRATA'="All_areas")
Region = "user"
DateFile = paste0('/Users/jiecao/Desktop/UW_work/model_devp/comparison')
#dir.create(DateFile)

DF = survey_data
Data_Geostat = cbind( "Size_class"=DF[,"size_class"], "Year"=DF[,"year"], "Catch_N"=DF[,"count"], "AreaSwept_km2"=DF[,"area_swept"], "Vessel"=0, "Lat"=DF[,"lat"], "Lon"=DF[,"lon"] )

pander::pandoc.table( head(Data_Geostat), digits=3 )
pander::pandoc.table( tail(Data_Geostat), digits=3 )

load('/Users/jiecao/Google Drive/UW_work/snow_crab/simulator/shrimp_example/shrimp_grid.rda')
Extrapolation_List = SpatialDeltaGLMM::Prepare_Extrapolation_Data_Fn( Region=Region, strata.limits=strata.limits, input_grid = input_grid )

Spatial_List = SpatialDeltaGLMM::Spatial_Information_Fn( grid_size_km=50, n_x=n_x, Method="Mesh", Lon=Data_Geostat[,'Lon'], Lat=Data_Geostat[,'Lat'], Extrapolation_List=Extrapolation_List, randomseed=Kmeans_Config[["randomseed"]], nstart=Kmeans_Config[["nstart"]], iter.max=Kmeans_Config[["iter.max"]], DirPath=DateFile )
# Add knots to Data_Geostat
Data_Geostat = cbind( Data_Geostat, "knot_i"=Spatial_List$knot_i )

# aggregate catch
# if(sampling_error==TRUE){
#   catch_raw = catch_data_error
# }else{
#   catch_raw = read.csv('/Users/jiecao/Desktop/UW_work/model_devp/SSST/Simulation_Data/catch_data.csv')
# }

catch_raw = catch_data
n_year = length(unique(Data_Geostat[,'Year']))

loc_catch = SpatialDeltaGLMM::Convert_LL_to_UTM_Fn(Lon = catch_raw$lon, Lat = catch_raw$lat, zone = Extrapolation_List$zone, flip_around_dateline = Extrapolation_List$flip_around_dateline)
NN_catch  = RANN::nn2(data = Spatial_List$Kmeans$centers[, c("E_km", "N_km")], query = loc_catch[, c("X", "Y")], k = 1)
catch_raw = cbind(catch_raw,NN_catch$nn.idx)
c_pkt = array(0, dim=c(n_p,Spatial_List$MeshList$anisotropic_mesh$n,n_year) )

n_boundary = Spatial_List$MeshList$anisotropic_mesh$n-nrow(Spatial_List$Kmeans$centers)

comp_temp = matrix(NA, nrow = n_x, ncol = n_p)
tot_c_temp = c()

for (i in 1:n_year){
  temp = subset(catch_raw,year ==i)
  
  #if(!sampling_error)
  #  {
      for (k in 1:n_p)
        {
          temp1 = tapply(temp[, (4+k)], INDEX = factor(temp[,'NN_catch$nn.idx'], levels = 1:nrow(Spatial_List$Kmeans$centers)), FUN = sum)
          c_pkt[k,,i] = c(temp1,rep(0,n_boundary))
        }
  #  }
  
  # if (sampling_error) 
  #   {
  #     tot_c_temp = tapply(temp[, (6+n_p)], INDEX = factor(temp[,'NN_catch$nn.idx'], levels = 1:nrow(Spatial_List$Kmeans$centers)), FUN = sum)
  #     sample.comp[,,i][,n_p+2] = temp[,'NN_catch$nn.idx'][sample.comp[,,i][,n_p+1]]
  #     temp_1 = sample.comp[,,i]
  #     for (k in 1:n_p)
  #       {
  #         comp_temp[,k] = tapply(temp_1[, k], INDEX = factor(temp_1[,n_p+2], levels = 1:nrow(Spatial_List$Kmeans$centers)), FUN = sum)
  #       
  #             if (any(is.na(comp_temp[,k]))){
  #                   na.id = which(is.na(comp_temp[,k])) 
  #                   distNN = as.matrix(cal_dist(Spatial_List$Kmeans$centers,k=4)$nearest_n[,na.id])
  #                   na.n = length(comp_temp[,k][is.na(comp_temp[,k])])
  #                         for (j in 1:na.n){
  #                                comp_temp[,k][na.id[j]] = mean(comp_temp[,k][distNN[,j]],na.rm = T)
  #                             }
  #                 }
  #        }
  #     comp_norl = comp_temp/apply(comp_temp,1,sum)
  #     c_pkt[,,i] = cbind(t(comp_norl*c(tot_c_temp)),matrix(0,nrow=n_p,ncol=n_boundary))
  #   }
  }

  c_pkt[is.na(c_pkt)] <- 0

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















