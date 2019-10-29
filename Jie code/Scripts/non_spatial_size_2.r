
require(PBSadmb)
setwd('/Users/jiecao/Desktop/UW_work/model_devp/Non-spatial_model/SSAP-EM')

n_reps=array(NA,dim=c(runtime+1,n_t,n_p))
#r_reps=matrix(NA,ncol=n_t,nrow=runtime+1)

#r_reps[1,]=pop.dyn.obj$Recruits

dorun=0
N_non_conv=0
scenario='non-spatial'

while(dorun<runtime)
{
  if (runtime>1)
  {
    pb <- txtProgressBar(min = 0, max = runtime-1, style = 3, char = "*")
    Sys.sleep(0.1)
    setTxtProgressBar(pb, dorun)
  }
  
  source('/Users/jiecao/Desktop/UW_work/model_devp/simulator/simulator_rcodes/simulator_genData2.r')
  
  if (dorun == 0){
    for(p in 1:n_p){
      n_reps[1,,p]=index[,p]/1000
    }
  }
  unlink("SSAP_M.STD")
  source('/Users/jiecao/Desktop/UW_work/model_devp/Non-spatial_model/SSAP-SIM/write-inputs.r')
  
  runAD("SSAP_M",verbose=T)
  fileNames=c("SSAP_M.STD")
  
  if (file.exists(fileNames)==TRUE)
  {
    file.copy("SSAP_M.REP", paste(scenario,'-',dorun,".","rep",sep=""),overwrite = T)
    file.copy("SSAP_M.STD", paste(scenario,'-',dorun,".","std",sep=""),overwrite = T)
    file.copy("SSAP_M.PAR", paste(scenario,'-',dorun,".","par",sep=""),overwrite = T)
    file.copy("SSAP_M.COR", paste(scenario,'-',dorun,".","cor",sep=""),overwrite = T)
    
    report<-read.admb(paste(scenario,'-',dorun,sep=""))
    for(p in 1:n_p){
      n_reps[dorun+2,,p]  = report$Abundance_at_Size[,p]
    }
    dorun=dorun+1
  }else
  {
    N_non_conv=N_non_conv+1
    #dorun=dorun+1
  }
}
save(n_reps, file=paste0('/Users/jiecao/Desktop/UW_work/model_devp/Non-spatial_model/SSAP-EM','/',"Sim.RData"))

for (p in 1:n_p){
  temp_p = c(n_reps[2:(runtime+1),,p])
  temp_true = rep(n_reps[1,,p],each=runtime)
  if(p==1) {
    abun = temp_p
    true_abun = temp_true
  }else{
    abun = c(abun,temp_p)
    true_abun = c(true_abun,temp_true)
  }
}
year = rep(rep(1:n_t,each=runtime),n_p)
size_bin = rep(1:n_p,each=n_t*runtime)
re = ((abun-true_abun)/true_abun)*100
nonspatial_data = data.frame(year,size_bin,true_abun,abun,re)

par(mfrow=c(1,5))
par(mar = c(0, 0, 0, 0), oma = c(4, 6, 4, 0.5))
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))
ylimit=c(-1,2)

for(p in 1:n_p){
  temp = subset(nonspatial_data,size_bin==p)
  boxplot(re~year,data=temp,ylim=c(-150,150))
}

total.abun <- apply(n_reps, c(1,2), sum)
total.abun.true = total.abun[1,]
total.abun.est = c(total.abun[2:(runtime+1),])
year2 = rep(1:n_t,each=runtime)
tatal.abun.data = data.frame(year2,total.abun.est)

par(mfrow=c(1,1))
boxplot(total.abun.est~year2,tatal.abun.data)
points(1:n_t,total.abun.true,col='red')

