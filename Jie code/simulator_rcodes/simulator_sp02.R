

### parallel 
# library(doParallel)
# cl = makeCluster(8)
# registerDoParallel(cl)
# stopCluster(cl)  

# example spatial points
library(sp);library(INLA);library(RandomFields)
library(magick); library(maptools)
ocean<-readShapeSpatial("/Users/jiecao/Desktop/UW_work/model_devp/simulator/50m_land.shp")
setwd('/Users/jiecao/Desktop/UW_work/model_devp/simulator')
SimFile = ('/Users/jiecao/Desktop/UW_work/model_devp/simulator/results')
SimFile = paste0(SimFile,'/',format(Sys.time(), "%d-%b-%Y %H.%M"))
dir.create(SimFile)

#source('Calc_Kmeans.r'); source('Calc_Anisotropic_Mesh.r'); 
#source('/Users/jiecao/Google Drive/UW_work/snow_crab/simulator/simulator_rcodes/grow_tran.r')
source('/Users/jiecao/Desktop/UW_work/model_devp/simulator/cal_GM.r')
LP = FALSE # Lognormal-Poisson for generating correlated counts

library(SpatialDeltaGLMM)
strata.limits <- data.frame('STRATA'="All_areas")
Region = "user"
load('/Users/jiecao/Desktop/UW_work/model_devp/simulator/shrimp_example/shrimp_grid.rda')

Extrapolation_List = SpatialDeltaGLMM::Prepare_Extrapolation_Data_Fn( Region=Region, strata.limits=strata.limits, input_grid = input_grid )
Data_Extrap = Extrapolation_List$Data_Extrap[Extrapolation_List$Data_Extrap$Area_km2>0,]
coords = cbind(Data_Extrap$Lon,Data_Extrap$Lat)
coords = coords[coords[,1]<0,]
area = Extrapolation_List$Area_km2_x
P1 = SpatialPoints(coords)

input_grid$ind = c(1:nrow(input_grid))
high_ind = subset(input_grid,stratum < 4070)$ind
low_ind = subset(input_grid,stratum >= 4070)$ind

# data input 
  set.seed(10)
  
  sp.points = P1
  n_s = nrow(sp.points@coords)
  loc_x = P1@coords
  
  #if(movement==TRUE) source('/Users/jiecao/Desktop/UW_work/model_devp/simulator/simulator_rcodes/movement.r')

  # define population structure
    #n_p = 20 # numnber of size classes  
    n_sex = 1 # number of sexes
    n_g = 4  # 1-immature; 2-mature; 3-new shell; 4-old shell
    #n_v = 2  # number of shell condition/stages
    
    low_s = 9  # lower boundary of size range
    up_s  = 34 # upper boundary of size range
    bin_width = 5
    
    n_p = (up_s-low_s)/bin_width # number of bins
    s_mid_points = seq(low_s+bin_width/2, up_s-bin_width/2, bin_width)
    
  # define time priod
    n_t = 20
    mat_at_size_female = c(3.50736e-08, 3.364911e-05, 0.03127495, 0.9687251, 0.9999664)
    weight_at_size = exp(-7.4302+3.007*log(s_mid_points))
    
    # scale pars controling decorrelation distance
    scale_r = 0.5
    scale_f = 0.5
    scale_m = 0.5
    scale_g = 0.5
    scale_Nint = 0.5
    
  # R
    
    ENV=c(-0.25,-0.3,-0.1,-0.95,0.3,-0.6,0.25,0.15,-0.5,0.1,0.7,-0.4,-0.55,-0.4,-0.7,-0.9,
          -0.75,-1.1,-0.55,-1.3,0.3,0.85,-0.2,-1.2,-0.2,-0.8,-0.45,-1.35,-1.65,-1.95)[1:n_t]
    
    envRelationSD=0.1
    meant=mean(ENV)
    sdtemp=(sum((ENV-meant)^2)/30)^0.5
    zscore1=(ENV-mean(ENV))/sdtemp
    #zscore=temp-mean(temp)
    set.seed(220)
    zscore=zscore1+rnorm(length(ENV),0,envRelationSD)
    R_mean_v = 1000000*exp(zscore)/n_s
    
    #R_mean_v = rpois(n_t,lambda=900*1e6/n_s); # per grid 
    if (LP==TRUE){
      R_mean_v = log(R_mean_v)
    }
    SD_R_v = rep(0.5,n_t) # spatial variation for each year 
    R_size_pert_v = c(1,0,rep(0,n_p-2)) # percentages for each size bin
    R_sex_ratio = 1 # male/total
    R_at_size_s = array (data = NA, dim = c(n_s,n_p,n_t))
    
    for (i in 1:n_t){
      spde_model_R <- RMgauss(var=SD_R_v[i]^2, scale=scale_r)
      R_kmeans = R_mean_v[i] * exp(RFsimulate(model = spde_model_R, x=loc_x[,1], y=loc_x[,2])@data[,1])
      if (LP==TRUE){
        R_kmeans = exp(R_kmeans)
        for (ii in 1:(length(R_kmeans))){
          R_kmeans[ii] = ifelse(R_kmeans[ii]<0.00001, 0, rpois(1, R_kmeans[ii]))
        }
      }
      R_at_size_s[,,i] = as.matrix(R_kmeans)%*%R_size_pert_v
    }

    R_male_at_size_s = R_sex_ratio*R_at_size_s
    R_female_at_size_s = (1-R_sex_ratio)*R_at_size_s
    
  # specify life history parameters
    # M
    M_spatial = FALSE # spatial or non-spatial
    M_mean_v = rep(0.25,n_t); 
    SD_M_v = rep(0.1,n_t) # spatial variation of M for each year
    M_size_v = rep(1,n_p) # size-specific M
  
    M_at_size_s = array (data = NA, dim = c(n_s,n_p,n_t))
    if (!M_spatial) SD_M_v = rep(0,n_t)
    for(i in 1:n_t){
      spde_model_M <- RMgauss(var=SD_M_v[i]^2, scale=scale_m)
      M_kmeans = M_mean_v[i] * exp(RFsimulate(model = spde_model_M, x=loc_x[,1], y=loc_x[,2])@data[,1])
      M_at_size_s[,,i] = as.matrix(M_kmeans)%*%M_size_v
    }
    
    # F
    # selectivity logistic
    ############# selectivities ###############
    # cal_sels <- function(pars_v,sizem_v)
    # {
    #     sels <- 1.0/(1.0+exp((pars_v[1]- sizem_v)* pars_v[2]));
    #     sels <- sels/max(sels)
    #   return(sels)
    #  }
    
    n_f = 1 # total catch; retained catch; bycatch
    pars_male = c(25,0.5);
    sels_male <- 1.0/(1.0+exp((pars_male[1]- s_mid_points)* pars_male[2]));
    
    #F_mean_v =  rnorm(n_t,mean=0.5,sd=0.1); 
    F_mean_v = c(0.2375996,0.2709,0.3267,0.1917,0.1979997,0.2781,0.2141997,0.1988997,0.1134,
      0.1376997,0.3060001,0.5354995,0.5805,0.4509,0.2070002,0.2592002,0.1314001,
      0.03420006,0.0783,0.1071,0.1233,0.0801,0.1917003,0.1980002,0.2780997,0.207,
      0.2591995,0.5355,0.5805,0.4509)[1:n_t] *2
    F_spatial = TRUE # spatial or non-spatial
    SD_F_v = rep(0.5,n_t) # spatial variation by year
    if (!F_spatial)  SD_F_v = rep(0.0,n_t)
    
    F_at_size_s = array (data = NA, dim = c(n_s,n_p,n_t,n_sex))
    F_kmeans = matrix(NA, nrow=n_s, ncol=n_t)
    
    for(i in 1:n_t){
      spde_model_F <- RMgauss(var=SD_F_v[i]^2, scale=scale_f)
      #F_kmeans[,i] = F_mean_v[i] * exp(RFsimulate(model = spde_model_F, x=loc_x[,1], y=loc_x[,2])@data[,1])
      F_kmeans[,i][high_ind] = 2*F_mean_v[i] * exp(RFsimulate(model = spde_model_F, x=loc_x[,1][high_ind], y=loc_x[,2][high_ind])@data[,1])
      F_kmeans[,i][low_ind] = 0.5*F_mean_v[i] * exp(RFsimulate(model = spde_model_F, x=loc_x[,1][low_ind], y=loc_x[,2][low_ind])@data[,1])
      F_at_size_s[,,i,1] = as.matrix(F_kmeans[,i])%*%sels_male
      #F_at_size_s[,,i,2] = as.matrix(F_kmeans[,i])%*%sels_female
    }
    
    # Growth
    G_spatial = FALSE # spatial or non-spatial
    
    n_grpar = 5
    L_inf_v = rep(36,n_t)
    SD_Linf_v = rep(0.1,n_t)

    K_v = rep(0.3,n_t)
    SD_K_v = rep(0.02,n_t)

    L_inf_sd_v = rep(1.0,n_t)
    SD_Linf_sd_v = rep(0.03,n_t)

    K_sd_v = rep(0.2,n_t)
    SD_K_sd_v = rep(0.01,n_t)

    rho_v = rep(0.9,n_t)
    SD_rho_v = rep(0.02,n_t)

    if(G_spatial){
      growth_par = array (data = NA, dim = c(n_s,n_grpar,n_t,n_sex)) 
      
      for(i in 1:n_t){
        spde_model_G_Linf   <- RMgauss(var=SD_Linf_v[i]^2, scale=scale_g)
        spde_model_G_K      <- RMgauss(var=SD_K_v[i]^2, scale=scale_r)
        spde_model_G_Linfsd <- RMgauss(var=SD_Linf_sd_v[i]^2, scale=scale_g)
        spde_model_G_Ksd    <- RMgauss(var=SD_K_sd_v[i]^2, scale=scale_g)
        spde_model_G_rho    <- RMgauss(var=SD_rho_v[i]^2, scale=scale_g)
        
        Linf_kmeans = L_inf_v[i]*exp(RFsimulate(model = spde_model_G_Linf, x=loc_x[,1], y=loc_x[,2])@data[,1])
        K_kmeans = K_v[i]*exp(RFsimulate(model = spde_model_G_K, x=loc_x[,1], y=loc_x[,2])@data[,1])
        L_inf_sd_kmeans = L_inf_sd_v[i]*exp(RFsimulate(model = spde_model_G_Linfsd, x=loc_x[,1], y=loc_x[,2])@data[,1])
        K_sd_kmeans = K_sd_v[i]*exp(RFsimulate(model = spde_model_G_Ksd, x=loc_x[,1], y=loc_x[,2])@data[,1])
        rho_kmeans = rho_v[i]*exp(RFsimulate(model = spde_model_G_rho, x=loc_x[,1], y=loc_x[,2])@data[,1])
        
        growth_par[,1,i,1] = Linf_kmeans
        growth_par[,2,i,1] = L_inf_sd_kmeans
        growth_par[,3,i,1] = K_kmeans
        growth_par[,4,i,1] = K_sd_kmeans
        growth_par[,5,i,1] = rho_kmeans
      }
    }
    
    #################################################################################
    
    N_at_size_ns = array (data = 0, dim = c(n_s,n_p,n_t))
    N_ns_total   = array (data = 0, dim = c(n_s,n_t))
    C_at_size_ns = array (data = 0, dim = c(n_s,n_p,n_t)) # total catch; retained catch; bycatch
    N_at_size_ns_female = array (data = 0, dim = c(n_s,n_p,n_t))
    
    
    # initial condition
    pia=c(0.1,0.2,0.3,0.25,0.15)
    pia=pia/sum(pia)
    
    Nini=3200000*pia/n_s
    
    SD_N_int <- rep(0.5,n_p)
    for (p in 1:n_p){
      spde_model_N_int   <- RMgauss(var=SD_N_int[p]^2, scale=scale_Nint)
      Nint_kmeans = Nini[p] * exp(RFsimulate(model = spde_model_N_int, x=loc_x[,1], y=loc_x[,2])@data[,1])
      N_at_size_ns[,p,1] = Nint_kmeans
      }
    
    N_at_size_ns[,,1] = R_male_at_size_s[,,1] + N_at_size_ns[,,1]    
    
    for (i in 2:n_t){
      if(!G_spatial){
        growth_matrix = cal_GM(s_mid_points[1],s_mid_points[n_p],bin_width,L_inf_v[i-1],L_inf_sd_v[i-1],K_v[i-1],K_sd_v[i-1],rho_v[i-1])
        #growth_matrix = cal_GM(low_s,up_s,bin_width,L_inf_v[i-1],L_inf_sd_v[i-1],K_v[i-1],K_sd_v[i-1],rho_v[i-1])
      }
      for (s in 1:n_s){
        if(G_spatial){
          growth_matrix = cal_GM(s_mid_points[1],s_mid_points[n_p],bin_width,growth_par[s,,i-1,1][1],growth_par[s,,i-1,1][2],growth_par[s,,i-1,1][3],growth_par[s,,i-1,1][4],growth_par[s,,i-1,1][5])
          if(nrow(growth_matrix)>n_p){
            growth_matrix = cal_GM(s_mid_points[1],s_mid_points[n_p],bin_width,growth_par[s,,i-1,1][1],growth_par[s,,i-1,1][2],growth_par[s,,i-1,1][3],growth_par[s,,i-1,1][4],growth_par[s,,i-1,1][5])[1:n_p,1:n_p]
          }
        }
        N_at_size_ns[s,,i] = R_male_at_size_s[s,,i]+(N_at_size_ns[s,,i-1]*exp(-M_at_size_s[s,,i-1]-F_at_size_s[s,,i-1,1]))%*%growth_matrix
        C_at_size_ns[s,,i-1] = area[s]*(1-exp(-F_at_size_s[s,,i-1,1]-M_at_size_s[s,,i-1]))*(N_at_size_ns[s,,i-1])*(F_at_size_s[s,,i-1,1]/(F_at_size_s[s,,i-1,1]+M_at_size_s[s,,i-1]))
        if (i==n_t)
          C_at_size_ns[s,,i] = area[s]*(1-exp(-F_at_size_s[s,,i,1]-M_at_size_s[s,,i]))*(N_at_size_ns[s,,i])*(F_at_size_s[s,,i,1]/(F_at_size_s[s,,i,1]+M_at_size_s[s,,i]))
      }
    }
    
    ssb_agg = c()
    
    for (i in 1:n_t){
      for (s in 1:n_s){
        N_at_size_ns_female[s,,i] = (area[s]*N_at_size_ns[s,,i]/1000)*mat_at_size_female*weight_at_size
      }
      ssb_agg[i] = sum(N_at_size_ns_female[,,i])
    }
    
    # system.time(
    # for (i in 1:n_t){
    #   foreach(s=1:n_s) %dopar% {
    #     N_at_size_ns_female[s,,i] = N_at_size_ns[s,,i]*mat_at_size_female
    #       }
    #   ssb_agg[i] = sum(N_at_size_ns_female[,,i])
    # }
    # )
    
  # population selection
    N_at_size_ns_area = array (data = 0, dim = c(n_s,n_p,n_t))
    catch_tp = matrix(NA,nrow=n_t,ncol=n_p); N_tp = matrix(NA,nrow=n_t,ncol=n_p)
    F_tp = matrix(NA,nrow=n_t,ncol=n_p)
    
    for(t in 1:n_t){
      for(l in 1:n_p){
        N_at_size_ns_area[,l,t] = N_at_size_ns[,l,t]*area
      }
    }
    
    for(t in 1:n_t){
      N_tp[t,] = apply(N_at_size_ns_area[,,t],2,sum,na.rm=TRUE)
      catch_tp[t,] = apply(C_at_size_ns[,,t],2,sum,na.rm=TRUE)
      for(l in 1:n_p){
        Catch <- function(F) (1-exp(-F-M_mean_v[t]))*N_tp[t,l]*(F/(F+M_mean_v[t]))
        Catchzero <- function(F,Catch0) Catch(F)-Catch0
        F_tp[t,l] = uniroot(Catchzero,c(0,20),Catch=catch_tp[t,l])$root
      }
      #F_tp[t,] = F_tp[t,] / max(F_tp[t,])
    }
    
    #sel_sim = sels_male/max(sels_male)

  # save results
    Sim = list('N_at_size_ns'=N_at_size_ns,'N_at_size_ns_female'=N_at_size_ns_female,'C_at_size_ns'=C_at_size_ns,'F_at_size_s'=F_at_size_s,'ssb_agg'=ssb_agg,'F_tp'=F_tp)
    save(Sim, file=paste0(SimFile,'/',"Sim.RData"))
    
# produce images
  
    Plot_annimation = function(P,data,n_t,Dir,name,breaks)
    {
      setwd(paste(Dir))
      for (i in 1:n_t){
        P$data = apply(data[,,i],1,sum)
        min=min(P$data);max=max(P$data)
        col=colorRampPalette(c("blue", "cyan", "yellow", "red"), bias=1)
        
        png(paste(name,'-Year',i,'.png',sep=''), height = 6, width = 6, units = 'in', res=600)
        rast.temp <- rasterize(P, rast, P$data, fun = mean)
        plot(ocean,col="dark gray",axes=F,xlim=xlim,ylim=ylim)
        image(rast.temp,col=col(11),axes=TRUE,add=T,xlim=xlim,ylim=ylim)
        #image.plot( add=TRUE, nlevel=11,legend.only=TRUE, horizontal=TRUE,
        #            col=col(11))
        mtext(paste('Density of all size classes','','-','','Year',i,sep=''),1,-2,adj=0.2)
        dev.off()
        img_temp = image_read(paste(name,'-Year',i,'.png',sep=''))
        if (i==1) img = img_temp
        img = c(img,img_temp)
      }
      img <- image_scale(img, "600")
      animation = image_animate(img, fps = 0.5, dispose = "previous")
      image_write(animation, "abundance.gif")
    }
    
    
############################################################################################
    library(raster)
    #library(fields)
# density plots
    area = 1
    Dir = SimFile
    name = 'den_sim' 
    P= P1
    den.data = log(N_at_size_ns)
    rast <- raster(ncol=80,nrow=80)
    extent(rast) <- extent(coords)
    
    xlim = c(-72,-67) ; ylim = c(40,46)
    
    col=colorRampPalette(c("blue", "cyan", "yellow", "red"), bias=1)
  
    break_tot <- c(50,51,53,55,56,57,58,59,60,61,62,65)
    breakpoints <- break_tot/13
    
    # Plot_annimation(P,den.data/area,n_t,Dir,name,breaks=break_tot)
    
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
        plot(ocean,col="dark gray",axes=F,xlim=xlim,ylim=ylim)
        image(rast.temp,col=col(11),axes=TRUE,add=T,xlim=xlim,ylim=ylim,zlim = c(min,max))
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
          P$data = den.data[,p,t]/area
          rast.temp <- rasterize(P, rast, P$data, fun = mean)
          plot(ocean,col="dark gray",axes=F,xlim=xlim,ylim=ylim)
          image(rast.temp,col=col(11),axes=TRUE,add=T,xlim=xlim,ylim=ylim,zlim = c(min,max))
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
      plot(ocean,col="dark gray",axes=F,xlim=xlim,ylim=ylim)
      image(rast.temp,col=col(10),axes=TRUE,add=T,xlim=xlim,ylim=ylim,zlim = c(min,max))
      mtext(paste('Year',t,sep=''),1,-2,adj=0.95,cex=0.8)
    }
    dev.off()
    
    # fishing mortality
    breakpoints_f=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.5,2)
    
    png(paste('F_sim','.png',sep=''), height = 4, width = 8, units = 'in', res=600)
    par(mfrow=c(4,5))
    par(mar = c(0, 0, 0, 0), oma = c(1, 1, 1, 0.5))
    par(tcl = -0.25)
    par(mgp = c(2, 0.6, 0))
    zlim = c(min(F_kmeans),max(F_kmeans))
    
      for(t in 1:n_t){
        P$data = F_kmeans[,t]
        rast.temp <- rasterize(P, rast, P$data, fun = mean)
        plot(ocean,col="dark gray",axes=F,xlim=c(-73,-67),ylim=c(40,45))
        image(rast.temp,col=col(14),axes=TRUE,add=T,xlim=c(-73,-67),ylim=c(40,45))
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
    #   plot(ocean,col="dark gray",axes=F,xlim=xlim,ylim=ylim)
    #   image(rast.temp,col=col(14),axes=TRUE,add=T,xlim=xlim,ylim=ylim)
    #   mtext(paste('Year',i,sep=''),1,-2,adj=0.95,cex=0.8)
    #   dev.off()
    #   img_temp = image_read(paste('F_sim','_year',i,'.png',sep=''))
    #   if (i==1) img = img_temp
    #   img = c(img,img_temp)
    # }
    # img <- image_scale(img, "900")
    # animation = image_animate(img, fps = 0.5, dispose = "previous")
    # image_write(animation, "F_sim.gif")
    
    
    dev.off()
    
    
    # recruitment and F plot
    png(paste('R_F.png'), height = 8, width = 8, units = 'in', res=600)
    par(mfrow=c(2,1))
    par(mar = c(0, 0, 0, 0), oma = c(4, 6, 4, 0.5))
    par(tcl = -0.25)
    par(mgp = c(2, 0.6, 0))
    
    plot(1:n_t,R_mean_v,xlab='Year',ylab='Recruitment',type='b',pch=16,xaxt='n')
    mtext('Recruitment',2,line=2.2)
    plot(1:n_t,2*F_mean_v,xlab='Year',ylab='Fishing mortality',type='b',pch=16)
    lines(1:n_t,0.5*F_mean_v,type='b',col='red',pch=16)
    mtext('Fishing mortality',2,line=2.2)
    legend(0.5,2.3,legend=c('A1','A2'),col=c('black','red'),lty=1,box.lty=0)
    
    dev.off()
    
    
    