
###################### deterministic data ######################
#simulate true population - northern shrimp in the GOM
source('/Users/jiecao/Desktop/UW_work/model_devp/simulator/simulator_rcodes/simulator_sp02.r')
#generate data
runNonspatial = TRUE
sampling_error=FALSE
logsigma_p = rep(0.25,n_p)
n_loc_sample = 100
source('/Users/jiecao/Desktop/UW_work/model_devp/simulator/simulator_rcodes/simulator_genData2.r')

# configurations in spatial_size_#.r
sampling_error = FALSE # catch
logcatch_p = rep(0.25,n_p) # catch
catchlogsd = 0.01 # used in EM
plot_figures = FALSE

source('/Users/jiecao/Desktop/UW_work/model_devp/SSST/Codes/spatial_size_2.r')


########################### run spatial assessment model  #################################
source('/Users/jiecao/Desktop/UW_work/model_devp/simulator/simulator_rcodes/simulator_sp02.r')

#generate data
sampling_error=TRUE

ObsModel_sim = c("PL", "Lognormal")[2]
logsigma_p = rep(0.3,n_p)
logcatch_p = rep(0.3,n_p) # catch
catchlogsd = 0.3 # used in EM
plot_figures = FALSE
runtime = 50 + 1
#n_loc_samples = c(100,200,300,500,1000,2000)[3]
n_loc_sample = 200
n_reps_sp=array(NA,dim=c(runtime,n_t,n_p))
determ=FALSE
runNonspatial=FALSE

iter=1
while (iter < runtime){
  unlink("Kmeans-50.RData")
  source('/Users/jiecao/Desktop/UW_work/model_devp/simulator/simulator_rcodes/simulator_genData2.r')
  source('/Users/jiecao/Desktop/UW_work/model_devp/SSST/Codes/spatial_size_2.r')
  if(all(abs(Opt$diagnostics[,'final_gradient'])<0.001)){
    n_reps_sp[iter+1,,]=Report$Index_tp
  }
  n_reps_sp[1,,]=index
  iter=iter+1
}

save(n_reps_sp, file=paste0(DateFile,'/',"n_reps_sp.RData"))

for (p in 1:n_p){
  temp_p = c(n_reps_sp[2:(runtime),,p])
  temp_true = rep(n_reps_sp[1,,p],each=runtime-1)
  if(p==1) {
    abun = temp_p
    true_abun = temp_true
  }else{
    abun = c(abun,temp_p)
    true_abun = c(true_abun,temp_true)
  }
}

re2 = ((abun-true_abun)/true_abun)*100
year = rep(rep(1:n_t,each=runtime-1),n_p)
size_bin = rep(1:n_p,each=n_t*(runtime-1))
spatial_data = data.frame(year,size_bin,true_abun,abun,re2)

png(paste('RE','.png',sep=''), height = 9, width = 18, units = 'in', res=600)
par(mfrow=c(1,5))
par(mar = c(0, 0, 0, 0), oma = c(4, 6, 4, 0.5))
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))
ylimit=c(-1,2)

for(p in 1:n_p){
  temp = subset(spatial_data,size_bin==p)
  boxplot(re2~year,data=temp,ylim=c(-30,20))
}
dev.off()

png(paste('total.abun','.png',sep=''), height = 9, width = 18, units = 'in', res=600)
total.abun <- apply(n_reps_sp, c(1,2), sum)
total.abun.true = total.abun[1,]
total.abun.est = c(total.abun[2:(runtime),])
year2 = rep(1:n_t,each=runtime-1)
tatal.abun.data = data.frame(year2,total.abun.est)

par(mfrow=c(1,1))
boxplot(total.abun.est~year2,tatal.abun.data,ylab='Year')
points(1:n_t,total.abun.true,col='red',pch=16)
dev.off()

########################### run non-spatial assessment model #################################
runNonspatial = TRUE
runtime=1
sampling_error=FALSE
ObsModel_sim = c("PL", "Lognormal")[2] # survey
logsigma_p = rep(0.2,n_p) # survey
#totcatch.cv = 0.1 # used in non-spatial EM
#survey.cv = 0.1 # used in non-spatial EM
ind.cv = sd(survey.data[,4])/mean(survey.data[,4])
totcatch.cv = sd(catch.data[,4])/mean(catch.data[,4])

logcatch_p = rep(0.25,n_p)
catchlogsd = 0.01 # used in spatial EM
n.eff = 100 # for catch composition
n_loc_sample = 100

source('/Users/jiecao/Desktop/UW_work/model_devp/SSST/Codes/non_spatial_size_2.r')

########################### run comparing spatial to non-spatial ##################################
#simulate true population - northern shrimp in the GOM
source('/Users/jiecao/Desktop/UW_work/model_devp/simulator/simulator_rcodes/simulator_sp02.r')

runNonspatial = TRUE
runtime=1
sampling_error=　TRUE
ObsModel_sim = c("PL", "Lognormal")[2] # survey
logsigma_p = rep(0.3,n_p) # survey; used in OM
logcatch_p = rep(0.6,n_p) # catch; used in OM

#ind.cv = 0.3 # used in non-spatial EM
#totcatch.cv = 0.5 # used in non-spatial EM
#survey.cv = cv(survey.data)
#totcatch.cv = cv(catch.data)
catchlogsd = 0.6 # used in spatial EM
n.eff = 60 # for catch composition
n_loc_sample = 60
plot_figures = FALSE

require(PBSadmb)
require(TMB)
setwd('/Users/jiecao/Desktop/UW_work/model_devp/comparison')

n_reps_nonsp=array(NA,dim=c(runtime+1,n_t,n_p))
n_reps_sp=array(NA,dim=c(runtime+1,n_t,n_p))
n_reps_sel = matrix(NA,nrow=runtime+1,ncol=n_p)
n_reps_ssb_nonsp = matrix(NA,nrow=runtime+1,ncol=n_t)
n_reps_ssb_sp = matrix(NA,nrow=runtime+1,ncol=n_t)

dorun=0
N_non_conv=0
scenario='non-spatial'

while(dorun<runtime)
{
  unlink("Kmeans-50.RData")
  source('/Users/jiecao/Desktop/UW_work/model_devp/simulator/simulator_rcodes/simulator_genData2.r')
  if (dorun == 0){
    for(p in 1:n_p){
      n_reps_nonsp[1,,p]=index[,p]/1000
      n_reps_sp[1,,p]=index[,p]/1000
    }
    n_reps_sel[1,] = sels_male
    n_reps_ssb_nonsp[1,]= ssb_agg
    n_reps_ssb_sp[1,]=  ssb_agg 
  }
  
  ### spatial model
  source('/Users/jiecao/Desktop/UW_work/model_devp/SSST/Codes/spatial_size_2_comparison2.r')
  
  if(all(abs(Opt$diagnostics[,'final_gradient'])<0.001)){
    n_reps_sp[dorun+2,,]=Report$Index_tp/1000
    
    for(t in 1:n_t){
      n_reps_ssb_sp[dorun+2,t] = sum((Report$Index_tp[t,]/1000)*mat_at_size_female*weight_at_size)
    }
 
    unlink("SSAP_M.STD")
    source('/Users/jiecao/Desktop/UW_work/model_devp/Non-spatial_model/SSAP-SIM/write-inputs.r')
    
    runAD("SSAP_M",verbose=F)
    fileNames=c("SSAP_M.STD")
    
    if (file.exists(fileNames)==TRUE)
    {
      file.copy("SSAP_M.REP", paste(scenario,'-',dorun,".","rep",sep=""),overwrite = T)
      file.copy("SSAP_M.STD", paste(scenario,'-',dorun,".","std",sep=""),overwrite = T)
      file.copy("SSAP_M.PAR", paste(scenario,'-',dorun,".","par",sep=""),overwrite = T)
      file.copy("SSAP_M.COR", paste(scenario,'-',dorun,".","cor",sep=""),overwrite = T)
      
      report<-read.admb(paste(scenario,'-',dorun,sep=""))
      for(p in 1:n_p){
        n_reps_nonsp[dorun+2,,p]  = report$Abundance_at_Size[,p]
        n_reps_sel[dorun+2,] =  report$Fleet_selectivity[1,]
        n_reps_ssb_nonsp[dorun+2,] =  report$Spawning_stock_Biomass_input
      }
      
      dorun=dorun+1
    }
  }else{
    N_non_conv = N_non_conv +1
  }
}

save.results = list(n_reps_nonsp=n_reps_nonsp,n_reps_sp=n_reps_sp,n_reps_sel=n_reps_sel,n_reps_ssb_nonsp=n_reps_ssb_nonsp,n_reps_ssb_sp=n_reps_ssb_sp)
save(save.results, file=paste0('/Users/jiecao/Desktop/UW_work/model_devp/comparison','/',"Sim.RData"))

attach(save.results);n_p=5;n_t=20;runtime=50

for (p in 1:n_p){
  temp_p = c(n_reps_nonsp[2:(runtime+1),,p])
  temp_true = rep(n_reps_nonsp[1,,p],each=runtime)
  if(p==1) {
    abun_non = temp_p
    true_abun = temp_true
  }else{
    abun_non = c(abun_non,temp_p)
    true_abun = c(true_abun,temp_true)
  }
}
year = rep(rep(1:n_t,each=runtime),n_p)
size_bin = rep(1:n_p,each=n_t*runtime)
re1 = ((abun_non-true_abun)/true_abun)*100
nonspatial_data = data.frame(year,size_bin,true_abun,abun_non,re1)

for (p in 1:n_p){
  temp_p = c(n_reps_sp[2:(runtime+1),,p])
  temp_true = rep(n_reps_sp[1,,p],each=runtime)
  if(p==1) {
    abun = temp_p
    true_abun = temp_true
  }else{
    abun = c(abun,temp_p)
    true_abun = c(true_abun,temp_true)
  }
}

re2 = ((abun-true_abun)/true_abun)*100
spatial_data = data.frame(year,size_bin,true_abun,abun,re2)

png(paste('comparison1','.png',sep=''), height = 7, width = 12, units = 'in', res=600)
par(mfrow=c(2,5))
par(mar = c(0, 0, 0, 0), oma = c(4, 6, 4, 0.5))
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))
ylimit=c(-1,2)

for(p in 1:n_p){
  temp = subset(nonspatial_data,size_bin==p)
  #RMSE_temp = (sum((aggregate(temp[-1], list(temp$year),mean)$re1/100)^2)/20)^0.5
  RMSE_temp = ((sum((aggregate(temp[-1], list(temp$year),mean)$re1/100)^2)/20)^0.5)*100
  RB_temp = mean(aggregate(temp[-1], list(temp$year),mean)$re1)
    
  boxplot(re1~year,data=temp,ylim=c(-30,30),axes=F)
  abline(h=0,col='gray',lwd=1.5)
  text(10,27,paste('RMSE=',round(RMSE_temp,2),'%','  ','RB=',round(RB_temp,2),'%',sep=''))
  
  if(p==1) {
    mtext('Spatially-implicit model',side=2,line=2.2)
    axis(2,at=seq(-30,30,10))
  }
  box()
  mtext(paste('Size class','',p),3)
}

for(p in 1:n_p){
  temp = subset(spatial_data,size_bin==p)
  #RMSE_temp = (sum((aggregate(temp[-1], list(temp$year),mean)$re2/100)^2)/20)^0.5
  RMSE_temp = ((sum((aggregate(temp[-1], list(temp$year),mean)$re2/100)^2)/20)^0.5)*100
  
  #RMSE_temp = mean(abs(aggregate(temp[-1], list(temp$year),mean)$re2))
  RB_temp = mean(aggregate(temp[-1], list(temp$year),mean)$re2)
  
  boxplot(re2~year,data=temp,ylim=c(-30,30),axes=F)
  abline(h=0,col='gray',lwd=1.5)
  text(10,27,paste('RMSE=',round(RMSE_temp,2),'%','  ','RB=',round(RB_temp,2),'%',sep=''))
  axis(1,at=seq(5,20,5))
  if(p==1) {
    axis(2,at=seq(-30,30,10))
    mtext('Spatiotemporal model',side=2,line=2.2)
  }
  box()
}
mtext('Year',side=1,outer=T,line=2.2)
dev.off()

####################################################################

total.abun.non <- apply(n_reps_nonsp, c(1,2), sum)
total.abun.true = total.abun.non[1,]
total.abun.estnon = c(total.abun.non[2:(runtime+1),])
year2 = rep(1:n_t,each=runtime)
re = 100*(total.abun.estnon - rep(total.abun.true,each=runtime))/rep(total.abun.true,each=runtime)
tatal.abun.data.non = data.frame(year2,total.abun.estnon,re)

total.abun <- apply(n_reps_sp, c(1,2), sum)
total.abun.true = total.abun[1,]
total.abun.est = c(total.abun[2:(runtime+1),])
re = 100*(total.abun.est - rep(total.abun.true,each=runtime))/rep(total.abun.true,each=runtime)
year2 = rep(1:n_t,each=runtime)
tatal.abun.data = data.frame(year2,total.abun.est,re)

ssb_nonsp = c(n_reps_ssb_nonsp[2:(runtime+1),])
ssb_re_non = 100*(ssb_nonsp - rep(n_reps_ssb_nonsp[1,],each=runtime))/rep(n_reps_ssb_nonsp[1,],each=runtime)

ssb_sp = c(n_reps_ssb_sp[2:(runtime+1),])
ssb_re_sp = 100*(ssb_sp - rep(n_reps_ssb_sp[1,],each=runtime))/rep(n_reps_ssb_sp[1,],each=runtime)

ssb.nansp.data = data.frame(year2, ssb_nonsp, ssb_re_non)
ssb.sp.data = data.frame(year2, ssb_sp, ssb_re_sp)

png(paste('comparison2','.png',sep=''), height = 6, width = 8, units = 'in', res=600)
par(mfrow=c(2,2))
par(mar = c(0, 0, 0, 0), oma = c(4, 6, 4, 0.5))
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))

boxplot(re~year2,tatal.abun.data.non,ylim=c(-20,20),axes=F)
abline(h=0,col='gray',lwd=1.5)
axis(2,at=seq(-30,30,10))
box()
mtext(paste('Spatially-implicit model'),3)
mtext('Total abundance',side=2,line=2.2)

boxplot(re~year2,tatal.abun.data,ylim=c(-20,20),axes=F)
abline(h=0,col='gray',lwd=1.5)
box()
mtext(paste('Spatiotemporal model'),3)

boxplot(ssb_re_non~year2,ssb.nansp.data,ylim=c(-40,20),axes=F)
abline(h=0,col='gray',lwd=1.5)
axis(2,at=seq(-40,10,10))
axis(1,at=seq(1,20,1))
box()
mtext('Spawning stock biomass',side=2,line=2.2)

boxplot(ssb_re_sp~year2,ssb.sp.data,ylim=c(-40,20),axes=F)
abline(h=0,col='gray',lwd=1.5)
axis(1,at=seq(1,20,1))
box()
mtext('Year',side=1,outer=T,line=2.2)
dev.off()

# selectivity plot
par(mfrow=c(1,1))
library(fanplot)
plot(seq(1,n_p,1),F_tp[3,],xlab="Size bins",ylab="Selectivity",
     type='l',xlim=c(1,5),ylim=c(0,1.1))
fan(data = n_reps_sel[2:(runtime+1),], start = 1)

###########################################
# population selection knot-level OM
# F_knot = array (data = NA, dim = c(n_x,n_p,n_t))
# N_at_size_knot = array (data = NA, dim = c(n_x,n_p,n_t))
# loc_n = SpatialDeltaGLMM::Convert_LL_to_UTM_Fn(Lon = loc_x[,1], Lat = loc_x[,2], zone = Extrapolation_List$zone, flip_around_dateline = Extrapolation_List$flip_around_dateline)
# NN_n  = RANN::nn2(data = Spatial_List$Kmeans$centers[, c("E_km", "N_km")], query = loc_n[, c("X", "Y")], k = 1)
# 
# for (l in 1:n_p){
#   for (t in 1:n_t){
#     N_at_size_knot[,l,t] = tapply(N_at_size_ns_area[,l,t], INDEX = factor(NN_n$nn.idx, levels = 1:nrow(Spatial_List$Kmeans$centers)), FUN = sum)
#       for (x in 1:n_x){
#         Catch <- function(F) (1-exp(-F-M_mean_v[t]))*N_at_size_knot[x,l,t]*(F/(F+M_mean_v[t]))
#         Catchzero <- function(F,Catch0) Catch(F)-Catch0
#         F_knot[x,l,t] = uniroot(Catchzero,c(0,10),Catch=c_pkt[l,x,t])$root
#       }
#   }
# }
# 
# Fatsize_knot_est = array (data = NA, dim = c(n_x,n_p,n_t))
# sels_est <- 1.0/(1.0+exp((ParHat$select[1]- s_mid_points)* ParHat$select[2]));
# F_knot_est = exp(Save$Report$logF_male_kt)[1:n_x,]
# 
# for (t in 1:n_t){
#   for (x in 1:n_x){
#     Fatsize_knot_est[x,,t] = F_knot_est[x,t]*sels_est
#   }
# }
# 
# par(mfrow=c(5,10))
# par(mar = c(0, 0, 0, 0), oma = c(4, 6, 4, 0.5))
# par(tcl = -0.25)
# par(mgp = c(2, 0.6, 0))
# 
# for(l in 1:n_p){
#   for(t in 1:n_t){
#     plot(Fatsize_knot_est[,l,t],F_knot[,l,t],pch=16,xaxt='n',yaxt='n')
#   }
# }

######## population selection

F_nonsp = array(NA, dim = c(runtime,n_p,n_t))
F_mean2 = matrix(NA, nrow = n_t, ncol = n_p)
F_5 = matrix(NA, nrow = n_t, ncol = n_p)
F_95 = matrix(NA, nrow = n_t, ncol = n_p)


for (iter in 1:runtime){
  report<-read.admb(paste(scenario,'-',iter-1,sep=""))
  for(t in 1:n_t){
    F_nonsp[iter,,t] = report$Fishing_Mortality[t]*report$Fleet_selectivity[t,]
  }
}

for(t in 1:n_t){
  for(p in 1:n_p){
    F_5[t,p] = quantile(F_nonsp[,p,t],0.05)
    F_95[t,p] = quantile(F_nonsp[,p,t],0.95)
    F_mean2[t,p] = quantile(F_nonsp[,p,t],0.5)
  }
}

png(paste('comparison3','.png',sep=''), height = 6, width = 8, units = 'in', res=600)
par(mfrow=c(4,5))
par(mar = c(0, 0, 0, 0), oma = c(4, 6, 4, 0.5))
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))

for(t in 1:n_t){
  plot(1:n_p,F_mean2[t,],type="l",ylim=c(0,0.8),axes=FALSE,ylab='',xlab='')
  text(2, 0.6, paste('Year','',t))
  if(t %in% c(1,6,11,16)) axis(2,at=seq(0.1,0.8,0.1))
  if(t %in% c(16,17,18,19,20)) axis(1,at=seq(1,5,1))
  
  polygon(c(1:n_p,n_p:1),c(F_5[t,],rev(F_95[t,])),col="gray",border=NA)
  lines(1:n_p,F_mean2[t,],col='black',lwd=1.2)
  lines(1:n_p,F_tp[t,],col='red',lwd=1.2)
  box()
}

mtext('Size class',side=1,outer=T,line=2.2)
mtext('Fishing mortality',side=2,outer=T,line=2.2)
dev.off()

















