
###################### deterministic data ######################
#simulate true population - snow crab in the East Bering Sea
movement = FALSE
system.time(source('Jie code/simulator_rcodes/simulator_sp01.r'))

#generate data
sampling_error=FALSE
n_loc_sample = 200
source('simulator_genData.r')

# configurations in spatial_size_#.r
determ = TRUE
sampling_error = FALSE # catch
logcatch_p = c(0.25,0.25,0.25,0.25,0.25) # catch
logsigma_p = c(0.25,0.25,0.25,0.25,0.25)
catchlogsd = 0.01 # used in EM
plot_figures = TRUE

system.time(source('/Users/jiecao/Desktop/UW_work/model_devp/SSST/Codes/spatial_size_1.r'))

####################### stochastic data #########################
movement = FALSE
source('/Users/jiecao/Desktop/UW_work/model_devp/simulator/simulator_rcodes/simulator_sp01.r')

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

iter=1
while (iter < runtime){
  unlink("Kmeans-50.RData")
  source('/Users/jiecao/Desktop/UW_work/model_devp/simulator/simulator_rcodes/simulator_genData.r')
  source('/Users/jiecao/Desktop/UW_work/model_devp/SSST/Codes/spatial_size_1.r')
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


