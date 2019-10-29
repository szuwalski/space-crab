
get_F_pkt <- function(dir,n_p,n_t,iter){
  setwd(dir)
  results = array(NA, dim=c(n_t,n_p,iter))
  for(i in 1:iter){
    load(paste('Save',i-1,'.RData',sep=''))
    catch_tp = t(apply(Save$Report$chat_pkt,c(1,3),sum))
    N_tp = Save$Report$Index_tp
    
    for(t in 1:n_t){
      for(l in 1:n_p){
        Catch <- function(F) (1-exp(-F-M_mean_v[t]))*N_tp[t,l]*(F/(F+M_mean_v[t]))
        Catchzero <- function(F,Catch0) Catch(F)-Catch0
        results[t,l,i] = uniroot(Catchzero,c(0,20),Catch=catch_tp[t,l])$root
      }
    }
  }
  return(results)
}
# population selection spatial EM
runtime=50
n_p=5
n_t=20
M_mean_v = rep(0.25,20)
dir = c('/Users/jiecao/Desktop/UW_work/model_devp/results_comparison')
F_tp = get_F_pkt(dir=dir,n_p=5,n_t=20,iter=runtime)

for (i in 1:runtime){
  for (t in 1:n_t){
    F_tp[t,,i] <- F_tp[t,,i]/max(F_tp[t,,i])
  }
}

F_mean_sp = matrix(NA, nrow = n_t, ncol = n_p)
F_5_sp = matrix(NA, nrow = n_t, ncol = n_p)
F_95_sp = matrix(NA, nrow = n_t, ncol = n_p)

for(t in 1:n_t){
  for(p in 1:n_p){
    F_5_sp[t,p] = quantile(F_tp[t,p,],0.05)
    F_95_sp[t,p] = quantile(F_tp[t,p,],0.95)
    F_mean_sp[t,p] = quantile(F_tp[t,p,],0.5)
  }
}

######## population selection nonspatial EM
source('/Users/jiecao/Desktop/UW_work/model_devp/Non-spatial_model/SSAP-SIM/reptoRlist.r')
F_nonsp = array(NA, dim = c(runtime,n_p,n_t))
F_mean2 = matrix(NA, nrow = n_t, ncol = n_p)
F_5 = matrix(NA, nrow = n_t, ncol = n_p)
F_95 = matrix(NA, nrow = n_t, ncol = n_p)

F_sp = array(NA, dim = c(runtime,n_p,n_t))

for (iter in 1:runtime){
  report<-read.admb(paste(dir,'/non-spatial','-',iter-1,sep=""))
  for(t in 1:n_t){
    #F_nonsp[iter,,t] = report$Fishing_Mortality[t]*report$Fleet_selectivity[t,]
    F_nonsp[iter,,t] = report$Fleet_selectivity[t,]/max(report$Fleet_selectivity[t,])
  }
}
for(t in 1:n_t){
  for(p in 1:n_p){
    F_5[t,p] = quantile(F_nonsp[,p,t],0.05)
    F_95[t,p] = quantile(F_nonsp[,p,t],0.95)
    F_mean2[t,p] = quantile(F_nonsp[,p,t],0.5)
  }
}

# F_tp_true
load("/Users/jiecao/Desktop/UW_work/model_devp/simulator/results/15-May-2018 14.55/Sim.RData")
F_tp_true = Sim$F_tp

for (t in 1:n_t){
  F_tp_true[t,] <- F_tp_true[t,]/max(F_tp_true[t,])
}

png(paste('figure6','.png',sep=''), height = 7, width = 8, units = 'in', res=600)
par(mfrow=c(4,5))
par(mar = c(0, 0, 0, 0), oma = c(4, 6, 4, 0.5))
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))

for(t in 1:n_t){
  plot(1:n_p,F_mean2[t,],type="l",ylim=c(0,1.05),axes=FALSE,ylab='',xlab='')
  text(2, 0.95, paste('Year','',t))
  if(t %in% c(1,6,11,16)) axis(2,at=seq(0.1,1,0.1))
  if(t %in% c(16,17,18,19,20)) axis(1,at=seq(1,5,1))
  
  polygon(c(1:n_p,n_p:1),c(F_5[t,],rev(F_95[t,])),col="red",border=NA)
  lines(1:n_p,F_mean2[t,],col='red',lwd=1.2)
  
  polygon(c(1:n_p,n_p:1),c(F_5_sp[t,],rev(F_95_sp[t,])),col="yellow",border=NA)
  lines(1:n_p,F_mean_sp[t,],col='yellow',lwd=1.2)
  lines(1:n_p,F_tp_true[t,],col='black',lwd=1.2,lty=3)
  
  box()
}

mtext('Size class',side=1,outer=T,line=2.2)
mtext('Selectivity',side=2,outer=T,line=2.2)
dev.off()

####################################################################
####################################################################

#png(paste('comparison3','.png',sep=''), height = 6, width = 8, units = 'in', res=600)
par(mfrow=c(4,5))
par(mar = c(0, 0, 0, 0), oma = c(4, 6, 4, 0.5))
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))

for(t in 1:n_t){
  plot(1:n_p,F_mean_sp[t,],type="l",ylim=c(0,0.8),axes=FALSE,ylab='',xlab='')
  text(2, 0.6, paste('Year','',t))
  if(t %in% c(1,6,11,16)) axis(2,at=seq(0.1,0.8,0.1))
  if(t %in% c(16,17,18,19,20)) axis(1,at=seq(1,5,1))
  
  polygon(c(1:n_p,n_p:1),c(F_5_sp[t,],rev(F_95_sp[t,])),col="gray",border=NA)
  lines(1:n_p,F_mean_sp[t,],col='gray',lwd=1.2)
  lines(1:n_p,F_tp_true[t,],col='red',lwd=1.2)
  box()
}

mtext('Size class',side=1,outer=T,line=2.2)
mtext('Fishing mortality',side=2,outer=T,line=2.2)
#dev.off()

