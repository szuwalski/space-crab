
# example spatial points
library(sp);library(INLA);library(RandomFields)
library(magick); library(maptools)
ocean<-readShapeSpatial("/Users/jiecao/Desktop/UW_work/model_devp/simulator/50m_land.shp")
setwd('/Users/jiecao/Desktop/UW_work/model_devp/simulator')
SimFile = ('/Users/jiecao/Desktop/UW_work/model_devp/simulator/results')
SimFile = paste0(SimFile,'/',format(Sys.time(), "%d-%b-%Y %H.%M"))
dir.create(SimFile)
#source('Calc_Kmeans.r'); source('Calc_Anisotropic_Mesh.r'); 
source('/Users/jiecao/Desktop/UW_work/model_devp/simulator/simulator_rcodes/grow_tran.r')

#movement = FALSE
set.seed(10)
# scale pars controling decorrelation distance
scale_r = 1.5
scale_f = 1.2
scale_m = 1.5
scale_g = 1.5

library(VAST)
strata.limits <- data.frame('STRATA'="All_areas")
Region = "Eastern_Bering_Sea"
Extrapolation_List = SpatialDeltaGLMM::Prepare_Extrapolation_Data_Fn( Region=Region, strata.limits=strata.limits )
Data_Extrap = Extrapolation_List$Data_Extrap[Extrapolation_List$Data_Extrap$Area_in_survey_km2>0,]
coords = cbind(Data_Extrap$Lon,Data_Extrap$Lat)
coords = coords[coords[,1]<0,]

x=coords[,1]
y=coords[,2]
grid_dataframe = data.frame(x=x,y=y)

P1 = SpatialPoints(coords)

# data input 

sp.points = P1
n_s = nrow(sp.points@coords)
loc_x = P1@coords

if(movement==TRUE) source('/Users/jiecao/Desktop/UW_work/model_devp/simulator/simulator_rcodes/movement.r')

# define population structure
n_sex = 1 # number of sexes: male
n_g = 2  # 1-immature; 2-mature

low_s = 25  # lower boundary of size range
up_s  = 125 # upper boundary of size range
bin_width = 20 

n_p = (up_s-low_s)/bin_width # number of bins
s_mid_points = seq(low_s+bin_width/2, up_s-bin_width/2, bin_width)

# define time priod
n_t = 10
mat_at_size_male = c(0,0.05,0.1,0.4,1)

# R
R_mean_v = rlnorm(n_t,log(2e7/n_s),0.3); # per grid 

SD_R_v = rep(0.3,n_t) # spatial variation for each year 
R_size_pert_v = c(1,rep(0,n_p-1)) # percentages for each size bin
R_at_size_s = array (data = NA, dim = c(n_s,n_p,n_t))

for (i in 1:n_t){
  spde_model_R <- RMgauss(var=SD_R_v[i]^2, scale=scale_r)
  R_kmeans = R_mean_v[i] * exp(RFsimulate(model = spde_model_R, x=loc_x[,1], y=loc_x[,2])@data[,1])
  R_at_size_s[,,i] = as.matrix(R_kmeans)%*%R_size_pert_v
}

R_male_at_size_s = R_at_size_s

# specify life history parameters
# M
M_spatial = FALSE # spatial or non-spatial
M_mean_v = rep(0.23,n_t); 
SD_M_v = rep(0.1,n_t) # spatial variation of M for each year
M_size_v = rep(1,n_p) # size-specific M

M_at_size_s = array (data = NA, dim = c(n_s,n_p,n_t))
if (!M_spatial) {
  M_at_size_s = array (data = 0.23, dim = c(n_s,n_p,n_t))
}else{
  for(i in 1:n_t){
    spde_model_M <- RMgauss(var=SD_M_v[i]^2, scale=scale_m)
    M_kmeans = M_mean_v[i] * exp(RFsimulate(model = spde_model_M, x=loc_x[,1], y=loc_x[,2])@data[,1])
    M_at_size_s[,,i] = as.matrix(M_kmeans)%*%M_size_v
  }
}
  
# F

# selectivity logistic
############# selectivities ##############

n_f = 1 # total catch; retained catch; bycatch
#pars_male = c(105,0.2) 
pars_male = c(75,0.05);
sels_male <- 1.0/(1.0+exp((pars_male[1]- s_mid_points)* pars_male[2]));
#sels_male = c(0.05,0.1,0.5,0.8,1)

F_mean_v =  rnorm(n_t,mean=0.5,sd=0.2); 
F_spatial = TRUE # spatial or non-spatial
SD_F_v = rep(0.2,n_t) # spatial variation by year
if (!F_spatial)  SD_F_v = rep(0.0,n_t)

F_at_size_s = array (data = NA, dim = c(n_s,n_p,n_t,n_sex))
F_kmeans = matrix(NA, nrow=n_s, ncol=n_t)

for(i in 1:n_t){
  spde_model_F <- RMgauss(var=SD_F_v[i]^2, scale=scale_f)
  F_kmeans[,i] = F_mean_v[i] * exp(RFsimulate(model = spde_model_F, x=loc_x[,1], y=loc_x[,2])@data[,1])
  F_at_size_s[,,i,1] = as.matrix(F_kmeans[,i])%*%sels_male
}

# Growth
G_spatial = FALSE # spatial or non-spatial

n_grpar = 3
int.male = rep(15,n_t); int.female = rep(1,n_t)
slope.male = rep(1,n_t); slope.female = rep(1.5,n_t)
beta.male = rep(0.1,n_t); beta.female = rep(0.5,n_t)

sd.int.male = rep(0.2,n_t); sd.int.female = rep(0.2,n_t)
sd.slope.male = rep(0.2,n_t); sd.slope.female = rep(0.2,n_t)
sd.beta.male = rep(0.2,n_t); sd.beta.female = rep(0.2,n_t)

if (!G_spatial){
  GM = growth_trans(int.male[1],slope.male[1],beta.male[1],low_s,up_s,bin_width)
}else{
  growth_par = array (data = NA, dim = c(n_s,n_grpar,n_t,n_sex)) # 1-male; 2-female
  for(i in 1:n_t){
    spde_model_G_int_male        <- RMgauss(var=sd.int.male[i]^2, scale=scale_g)
    spde_model_G_slope_male      <- RMgauss(var=sd.slope.male[i]^2, scale=scale_g)
    spde_model_G_beta_male       <- RMgauss(var=sd.beta.male[i]^2, scale=scale_g)
    int_kmeans_male = int.male[i] + RFsimulate(model = spde_model_G_int_male, x=loc_x[,1], y=loc_x[,2])@data[,1]
    slope_kmeans_male = slope.male[i] + RFsimulate(model = spde_model_G_slope_male, x=loc_x[,1], y=loc_x[,2])@data[,1]
    beta_kmeans_male = beta.male[i] + RFsimulate(model = spde_model_G_beta_male, x=loc_x[,1], y=loc_x[,2])@data[,1]
    growth_par[,1,i,1] = int_kmeans_male
    growth_par[,2,i,1] = slope_kmeans_male
    growth_par[,3,i,1] = beta_kmeans_male
    
    spde_model_G_int_female        <- RMgauss(var=sd.int.female[i]^2, scale=scale_g)
    spde_model_G_slope_female      <- RMgauss(var=sd.slope.female[i]^2, scale=scale_g)
    spde_model_G_beta_female       <- RMgauss(var=sd.beta.female[i]^2, scale=scale_g)
    int_kmeans_female = int.female[i] + RFsimulate(model = spde_model_G_int_female, x=loc_x[,1], y=loc_x[,2])@data[,1]
    slope_kmeans_female = slope.female[i] + RFsimulate(model = spde_model_G_slope_female, x=loc_x[,1], y=loc_x[,2])@data[,1]
    beta_kmeans_female = beta.female[i] + RFsimulate(model = spde_model_G_beta_female, x=loc_x[,1], y=loc_x[,2])@data[,1]
    growth_par[,1,i,2] = int_kmeans_female
    growth_par[,2,i,2] = slope_kmeans_female
    growth_par[,3,i,2] = beta_kmeans_female
  }
}
  
###################################################################################
# equilibrium state with no fishing mortality; M happens before growth
# n_g:  1-immature; 2-mature; 

Plot_equil = FALSE
n_t_equil = 10 # make sure n_t = n_t_equil otherwise growth and M don't match with the dimension
N_at_size_grid_male   = array (data = 0, dim = c(n_s,n_p,n_t_equil,n_g))
N_grid_total = array (data = 0, dim = c(n_s,n_t_equil,n_sex)) # 1-male; 2-female
N_at_size_grid_male[,,1,1]   = R_male_at_size_s[,,1]
N_grid_total[,1,1] = apply(N_at_size_grid_male[,,1,1],1,sum)

for (i in 2:n_t_equil)
{
  for (s in 1:n_s)
  {
    #growth_matrix_male = cal_GM(s_mid_points[1],s_mid_points[n_p],bin_width,growth_par[s,,i-1,1][1],growth_par[s,,i-1,1][2],growth_par[s,,i-1,1][3],growth_par[s,,i-1,1][4],growth_par[s,,i-1,1][5])
    if(!G_spatial){
      growth_matrix_male = GM
    }else{
      growth_matrix_male = growth_trans(growth_par[s,,i-1,1][1],growth_par[s,,i-1,1][2],growth_par[s,,i-1,1][3],low_s,up_s,bin_width)
    }
    N_at_size_grid_male[s,,i,1] = R_male_at_size_s[s,,i] + ((N_at_size_grid_male[s,,i-1,1]*exp(-M_at_size_s[s,,i-1]))%*%growth_matrix_male)*(1-mat_at_size_male)
    N_at_size_grid_male[s,,i,2] = ((N_at_size_grid_male[s,,i-1,1]*exp(-M_at_size_s[s,,i-1]))%*%growth_matrix_male)*mat_at_size_male+N_at_size_grid_male[s,,i-1,2]*exp(-M_at_size_s[s,,i-1])
    N_grid_total[s,i,1] = sum((N_at_size_grid_male[s,,i,1] + N_at_size_grid_male[s,,i,2]))
  }
}


############################################################

# initial condition
# n_g: 1-immature; 2-mature; 
N_at_size_ns_male   = array (data = 0, dim = c(n_s,n_p,n_t,n_g))
N_ns_total          = array (data = 0, dim = c(n_s,n_t,n_sex))
C_at_size_ns_male   = array (data = 0, dim = c(n_s,n_p,n_t)) # total catch; retained catch; bycatch
N_at_size_ns_male[,,1,]    = N_at_size_grid_male[,,n_t_equil,] 
N_at_size_ns_male[,,1,1]  = N_at_size_ns_male[,,1,1] - R_male_at_size_s[,,n_t_equil] + R_male_at_size_s[,,1]

# population dynamics - July 1st, survey -> fishery -> growth
alpha             = 0.62 # fishery time point  F_at_size_s

for (i in 2:n_t)
{   
  if (movement == TRUE){
    for (p in 1:n_p){
      N_at_size_ns_male[,,i-1,1][,p] = as.vector(Movement_m %*% N_at_size_ns_male[,,i-1,1][,p])
    }
  }
  
  for (s in 1:n_s)
  {
    
    # male
    #growth_matrix_male = cal_GM(s_mid_points[1],s_mid_points[n_p],bin_width,growth_par[s,,i-1,1][1],growth_par[s,,i-1,1][2],growth_par[s,,i-1,1][3],growth_par[s,,i-1,1][4],growth_par[s,,i-1,1][5])
    if(!G_spatial){
      growth_matrix_male = GM
    }else{
      growth_matrix_male = growth_trans(growth_par[s,,i-1,1][1],growth_par[s,,i-1,1][2],growth_par[s,,i-1,1][3],low_s,up_s,bin_width)
    }
    N_male_temp = (N_at_size_ns_male[s,,i-1,1]*exp(-M_at_size_s[s,,i-1]-F_at_size_s[s,,i-1,1]))%*%growth_matrix_male
    N_at_size_ns_male[s,,i,1] = R_male_at_size_s[s,,i] + N_male_temp*(1-mat_at_size_male)
    N_at_size_ns_male[s,,i,2] = N_male_temp*mat_at_size_male+N_at_size_ns_male[s,,i-1,2]*exp(-M_at_size_s[s,,i-1]-F_at_size_s[s,,i-1,1])
    N_grid_total[s,i,1] = sum((N_at_size_ns_male[s,,i,1] + N_at_size_ns_male[s,,i,2]))
    
    # catch
      C_at_size_ns_male[s,,i-1] = (1-exp(-F_at_size_s[s,,i-1,1]))*(N_at_size_ns_male[s,,i-1,1]+N_at_size_ns_male[s,,i-1,2])*exp(-alpha*M_at_size_s[s,,i-1])
    if (i == n_t)
    {
      C_at_size_ns_male[s,,i] = (1-exp(-F_at_size_s[s,,i,1]))*(N_at_size_ns_male[s,,i,1]+N_at_size_ns_male[s,,i,2])*exp(-alpha*M_at_size_s[s,,i])
    }
  }
}

# save results

Sim = list('N_at_size_ns_male'=N_at_size_ns_male,
           'C_at_size_ns_male'=C_at_size_ns_male,
           'F_at_size_s'=F_at_size_s,'n_s'=n_s,'n_t'=n_t,'n_p'=n_p)
save(Sim, file=paste0(SimFile,'/',"Sim.RData"))

# produce images

Dir = SimFile
P= P1
data = N_at_size_ns_male[,,,1] + N_at_size_ns_male[,,,2]

# data = C_at_size_ns_male[,,,1]
name = 'den'
# Plot_annimation(P,data,n_t,Dir,name)  

# Plot_annimation = function(P,data,n_t,Dir,name,breaks)
# {
#   setwd(paste(Dir))
#   for (i in 1:n_t){
#     P$data = apply(data[,,i],1,sum)
#     #min=min(P$data);max=max(P$data)
#     col=colorRampPalette(c("blue", "cyan", "yellow", "red"), bias=1)
#     png(paste(name,'-Year',i,'.png',sep=''), height = 6, width = 6, units = 'in', res=600)
#     rast.temp <- rasterize(P, rast, P$data, fun = mean)
#     plot(ocean,col="dark gray",axes=F,xlim=c(-180,-158),ylim=c(56,58))
#     image(rast.temp,col=col(12),axes=TRUE,breaks=breaks,add=T,xlim=c(-180,-150),ylim=c(55,65))
#     mtext(paste('Density of all size classes','','-','','Year',i,sep=''),1,-2,adj=0.2)
#     dev.off()
#     img_temp = image_read(paste(name,'-Year',i,'.png',sep=''))
#     if (i==1) img = img_temp
#     img = c(img,img_temp)
#   }
#   img <- image_scale(img, "600")
#   animation = image_animate(img, fps = 0.5, dispose = "previous")
#   image_write(animation, "abundance.gif")
# }

# rast <- raster(ncol=200,nrow=50)
# extent(rast) <- extent(coords)
# rast2 <- rasterize(P, rast, P$data, fun = mean)
# 
# par(mfrow=c(4,3))
# for (i in 1:n_t){
#   P$data = data[,,i]
#   min=min(P$data);max=max(P$data)
#   col=colorRampPalette(c("blue", "cyan", "yellow", "red"), bias=1)
#   #png(paste(name,'-Year',i,'.png',sep=''), height = 6, width = 6, units = 'in', res=600)
#   plot(ocean,col="dark gray",axes=T,xlim=c(-180,-155),ylim=c(55,60))
#   image(SpatialPixelsDataFrame(P,tolerance=0.51,P@data),col=col(10),
#         axes=TRUE,zlim = c(min,max),add=T,pch=16,
#         xlim=c(-180,-150),ylim=c(55,65))
#   mtext(paste('Year',i,sep=''),3,-2,adj=0.2)
#   #dev.off()
# }

############################################################################################
library(raster)
# density plots
area = 13.71962
Dir = SimFile
name = 'den_sim' 
P= P1
den.data = N_at_size_ns_male[,,,1] + N_at_size_ns_male[,,,2]
rast <- raster(ncol=200,nrow=100)
extent(rast) <- extent(coords)
col=colorRampPalette(c("blue", "cyan", "yellow", "red"), bias=1)
breakpoints <- c(0,0.05,0.1,0.2,0.3,0.5,0.7,0.9,1.1,1.3,1.5,2,2.5,3,3.5,4,4.5,5.5,6,6.5,10)


break_tot <- c(1,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,15)
#Plot_annimation(P,den.data/area,n_t,Dir,name,breaks=break_tot)

setwd(Dir)
png(paste(name,'-',Sys.Date(),'.png',sep=''), height = 9, width = 18, units = 'in', res=600)
par(mfrow=c(n_p,n_t))
par(mar = c(0, 0, 0, 0), oma = c(4, 4, 4, 0.5))
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))
for(p in 1:n_p){
  for(t in 1:n_t){
    P$data = den.data[,p,t]/area
    rast.temp <- rasterize(P, rast, P$data, fun = mean)
    plot(ocean,col="dark gray",axes=F,xlim=c(-180,-158),ylim=c(56,58))
    image(rast.temp,col=col(20),axes=TRUE,breaks=breakpoints,add=T,xlim=c(-180,-150),ylim=c(55,65),zlim = c(min,max))
    mtext(paste('Size class',p,',','Year',t,sep=''),1,-2,adj=1,cex=0.5)
    #filled.contour(x=x,y=y,z=den.data[,p,t]/area,levels=pretty(zlim,nlevels),nelevels=20,color.palette = cm.colors,col=color.palette(length(levels)-1))
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
    P$data = den.data[,p,t]/area
    rast.temp <- rasterize(P, rast, P$data, fun = mean)
    plot(ocean,col="dark gray",axes=F,xlim=c(-180,-158),ylim=c(56,58))
    image(rast.temp,col=col(20),axes=TRUE,breaks=breakpoints,add=T,xlim=c(-180,-150),ylim=c(55,65),zlim = c(min,max))
    mtext(paste('Year',t,sep=''),1,-2,adj=0.95,cex=0.8)
  }
  dev.off()
}

# recruitment
png(paste('recruitment.png'), height = 8, width = 6, units = 'in', res=600)
par(mfrow=c(4,3))
par(mar = c(0, 0, 0, 0), oma = c(4, 4, 4, 0.5))
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0)) 

for(t in 1:n_t){
  P$data = R_male_at_size_s[,1,t]/area
  rast.temp <- rasterize(P, rast, P$data, fun = mean)
  plot(ocean,col="dark gray",axes=F,xlim=c(-180,-158),ylim=c(56,58))
  image(rast.temp,col=col(20),axes=TRUE,breaks=breakpoints,add=T,xlim=c(-180,-150),ylim=c(55,65),zlim = c(min,max))
  mtext(paste('Year',t,sep=''),1,-2,adj=0.95,cex=0.8)
}
dev.off()

# fishing mortality
breakpoints_f=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.5,2)

png(paste('F_sim','.png',sep=''), height = 4, width = 8, units = 'in', res=600)
par(mfrow=c(2,5))
par(mar = c(0, 0, 0, 0), oma = c(1, 1, 1, 0.5))
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))

for(t in 1:n_t){
  P$data = F_kmeans[,t]
  rast.temp <- rasterize(P, rast, P$data, fun = mean)
  plot(ocean,col="dark gray",axes=F,xlim=c(-180,-158),ylim=c(56,58))
  image(rast.temp,col=col(14),axes=TRUE,breaks=breakpoints_f,add=T,xlim=c(-180,-150),ylim=c(55,65))
  mtext(paste('Year',t,sep=''),1,-2,adj=0.95,cex=0.8)
}

dev.off()

# animation 
# for (i in 1:n_t){
#   P$data = F_kmeans[,i]
#   col=colorRampPalette(c("blue", "cyan", "yellow", "red"), bias=1)
#   png(paste('F_sim','_year',i,'.png',sep=''), height = 6, width = 6, units = 'in', res=900)
#   par(mar = c(0, 0, 0, 0), oma = c(1, 1, 1, 0.5))
#   par(tcl = -0.25)
#   par(mgp = c(2, 0.6, 0))
#   rast.temp <- rasterize(P, rast, P$data, fun = mean)
#   plot(ocean,col="dark gray",axes=F,xlim=c(-180,-158),ylim=c(56,58))
#   image(rast.temp,col=col(14),axes=TRUE,breaks=breakpoints_f,add=T,xlim=c(-180,-150),ylim=c(55,65))
#   mtext(paste('Year',i,sep=''),1,-2,adj=0.95,cex=0.8)
#   dev.off()
#   img_temp = image_read(paste('F_sim','_year',i,'.png',sep=''))
#   if (i==1) img = img_temp
#   img = c(img,img_temp)
# }
# img <- image_scale(img, "900")
# animation = image_animate(img, fps = 0.5, dispose = "previous")
# image_write(animation, "F_sim.gif")



if(Plot_equil)
{
  setwd(SimFile)
  col=colorRampPalette(c("blue", "cyan", "yellow", "red"), bias=1)
  min=min(N_grid_total);max=max(N_grid_total)
  for (i in 1:n_t_equil)
  {
    P_male = P1; P_female = P1
    P_male$male    = N_grid_total[,i,1]
    P_female$femal = N_grid_total[,i,2]
    
    png(paste('Year',i,'.png',sep=''), height = 6, width = 12, units = 'in', res=900)
    par(mfrow=c(1,2))
    plot(ocean,col="dark gray",axes=T,xlim=c(-180,-165),ylim=c(55,60),main='Male')
    image(SpatialPixelsDataFrame(P_male,tolerance=0.3,P_male@data),col=col(10),
          axes=TRUE,zlim = c(min,max),add=T,
          xlim=c(-180,-150),ylim=c(55,65))
    mtext(paste('Year',i,sep=''),3,-2,adj=0.2)
    
    plot(ocean,col="dark gray",axes=T,xlim=c(-180,-165),ylim=c(55,60),main='Female')
    image(SpatialPixelsDataFrame(P_female,tolerance=0.3,P_female@data),col=col(10),
          axes=TRUE,zlim = c(min,max),add=T,
          xlim=c(-180,-150),ylim=c(55,65))
    mtext(paste('Year',i,sep=''),3,-2,adj=0.2)
    
    dev.off()
    
    img_temp = image_read(paste('Year',i,'.png',sep=''))
    
    if (i==1) img = img_temp
    img = c(img,img_temp)
  }
  
  img <- image_scale(img, "1200")
  animation = image_animate(img, fps = 2, dispose = "previous")
  image_write(animation, "abundance.gif")
}











