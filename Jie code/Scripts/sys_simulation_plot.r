
# systematic simulation experiment 

get_rmse <- function(dir,true_index,n_p,n_t,iter){
  setwd(dir)
  results = matrix(NA, nrow=iter, ncol=n_p)
  for(i in 1:iter){
    load(paste('Save',i,'.RData',sep=''))
    for (p in 1:n_p){
      results[i,p] = ((sum(((Save$Report$Index_tp[,p] - true_index[,p])/true_index[,p])^2)/n_t)^0.5)*100
    }
  }
  return(results)
}

iter = 50; n_p=5
results = array(NA, dim=c(iter,6,n_p)) # 1-crab 50; 2-crab 100; 3-crab 200; 4-shrimp 50; 5-shrimp 100; 6-shrimp 200

dir1 = '/Users/jiecao/Desktop/UW_work/model_devp/SSST/SSST_output/crab_K50'
dir2 = '/Users/jiecao/Desktop/UW_work/model_devp/SSST/SSST_output/crab_K100'
dir3 = '/Users/jiecao/Desktop/UW_work/model_devp/SSST/SSST_output/crab_K200'

dir4 = '/Users/jiecao/Desktop/UW_work/model_devp/SSST/SSST_output/iterations50'
dir5 = '/Users/jiecao/Desktop/UW_work/model_devp/SSST/SSST_output/iterations100'
dir6 = '/Users/jiecao/Desktop/UW_work/model_devp/SSST/SSST_output/iterations200'

# true index
# source('/Users/jiecao/Desktop/UW_work/model_devp/simulator/simulator_rcodes/simulator_sp02.r')
# N_at_size_ns_tot=array (data = 0, dim = c(n_s,n_p,n_t))
# index = matrix(0,nrow=n_t,ncol=n_p);area = Extrapolation_List$Area_km2_x
# for(t in 1:n_t){
#   for(p in 1:n_p){
#     for(s in 1:n_s){
#       N_at_size_ns_tot[s,p,t]=N_at_size_ns[s,p,t]*area[s]
#     }
#     index[t,p] = sum(N_at_size_ns_tot[,p,t])
#   }
# }
# true_index_shrimp = index
# 
# ###### crab example
# movement=FALSE
# source('/Users/jiecao/Desktop/UW_work/model_devp/simulator/simulator_rcodes/simulator_sp01.r')
# index = matrix(0,nrow=10,ncol=5)
# for(t in 1:n_t){
#   for(p in 1:n_p){
#     index[t,p] = sum(N_at_size_ns_male[,p,t,1] + N_at_size_ns_male[,p,t,2])
#   }
# }
# true_index_crab = index
# 

iter = 50
n_p = 5

results[,1,] = get_rmse(dir1,true_index_crab,n_p,10,iter)
results[,2,] = get_rmse(dir2,true_index_crab,n_p,10,iter)
results[,3,] = get_rmse(dir3,true_index_crab,n_p,10,iter)
  
results[,4,] = get_rmse(dir4,true_index_shrimp,n_p,20,iter)
results[,5,] = get_rmse(dir5,true_index_shrimp,n_p,20,iter)
results[,6,] = get_rmse(dir6,true_index_shrimp,n_p,20,iter)


  RE = c(c(results[,,1]),c(results[,,2]),c(results[,,3]),c(results[,,4]),c(results[,,5]))
  species = rep(c(rep('Snow crab',150),rep('Northern shrimp',150)),5)
  samples.n = rep(rep(c(50,100,200,50,100,200), each=50),5)
  size = rep(c(1,2,3,4,5),each=300)
  data.re = data.frame(RE=RE,species=species,samples.n=samples.n,size=size)
  
  data.re$species <- factor(data.re$species)
  data.re$samples.n <- factor(data.re$samples.n)
  data.re$size <- factor(data.re$size)

  ##############################################
  library(ggplot2)
  setwd('/Users/jiecao/Desktop')
  png(paste('figure7','.png',sep=''), height = 7, width = 12, units = 'in', res=600)
  version.labs <- c(`1`="Size class 1", `2`="Size class 2", `3`="Size class 3", `4`="Size class 4", `5`="Size class 5")
  ggplot(data.re, aes(x=factor(samples.n),y=RE)) + 
    #geom_jitter(position=position_jitter(width=0.3, height=0.2), aes(colour=factor(species)), alpha=0.6) +
    geom_boxplot(alpha = 0.5, show.legend = TRUE, aes(fill=factor(species))) + 
    facet_grid(.~size, labeller = as_labeller(version.labs)) +
    labs(x="Number of sampling locations", y="Root-mean-square error ") + 
    theme(axis.title.x = element_text(size=13,hjust=0.5),
          axis.title.y = element_text(size=13,vjust=1),
          axis.text.x = element_text(size=11,color='black'),
          axis.text.y = element_text(size=11,color='black'),
          legend.title = element_text(face="bold", color="black", size=14),
          legend.text = element_text(face="bold", color="black", size=12),
          strip.text.x = element_text(size=9, color="black", face="bold"),
          panel.background = element_blank()) +
    guides(color=guide_legend("Species"))
  dev.off()
  
  #########
  
  ###########################################################################################

  get_rb <- function(dir,true_index,n_p,n_t,iter){
    setwd(dir)
    results = matrix(NA, nrow=iter, ncol=n_p)
    for(i in 1:iter){
      load(paste('Save',i,'.RData',sep=''))
      for (p in 1:n_p){
        results[i,p] = (sum((Save$Report$Index_tp[,p] - true_index[,p])/true_index[,p])/n_t)*100
      }
    }
    return(results)
  }
  
  iter = 50; n_p=5
  results = array(NA, dim=c(iter,6,n_p)) # 1-crab 50; 2-crab 100; 3-crab 200; 4-shrimp 50; 5-shrimp 100; 6-shrimp 200
  
  iter = 50
  n_p = 5
  
  results[,1,] = get_rb(dir1,true_index_crab,n_p,10,iter)
  results[,2,] = get_rb(dir2,true_index_crab,n_p,10,iter)
  results[,3,] = get_rb(dir3,true_index_crab,n_p,10,iter)
  
  results[,4,] = get_rb(dir4,true_index_shrimp,n_p,20,iter)
  results[,5,] = get_rb(dir5,true_index_shrimp,n_p,20,iter)
  results[,6,] = get_rb(dir6,true_index_shrimp,n_p,20,iter)
  
  
  RE = c(c(results[,,1]),c(results[,,2]),c(results[,,3]),c(results[,,4]),c(results[,,5]))
  species = rep(c(rep('Snow crab',150),rep('Northern shrimp',150)),5)
  samples.n = rep(rep(c(50,100,200,50,100,200), each=50),5)
  size = rep(c(1,2,3,4,5),each=300)
  data.rb = data.frame(RE=RE,species=species,samples.n=samples.n,size=size)
  
  data.rb$species <- factor(data.rb$species)
  data.rb$samples.n <- factor(data.rb$samples.n)
  data.rb$size <- factor(data.rb$size)
  
  ##############################################
  library(ggplot2)
  setwd('/Users/jiecao/Desktop')
  png(paste('figure8','.png',sep=''), height = 7, width = 12, units = 'in', res=600)
  version.labs <- c(`1`="Size class 1", `2`="Size class 2", `3`="Size class 3", `4`="Size class 4", `5`="Size class 5")
  ggplot(data.re, aes(x=factor(samples.n),y=RE)) + 
    #geom_jitter(position=position_jitter(width=0.3, height=0.2), aes(colour=factor(species)), alpha=0.15) +
    geom_hline(yintercept = 0,col='gray') + 
    geom_boxplot(alpha = 0.5, show.legend = TRUE, aes(fill=factor(species))) + 
    facet_grid(.~size, labeller = as_labeller(version.labs)) +
    labs(x="Number of sampling locations", y="Relative bias") + 
    theme(axis.title.x = element_text(size=13,hjust=0.5),
          axis.title.y = element_text(size=13,vjust=1),
          axis.text.x = element_text(size=11,color='black'),
          axis.text.y = element_text(size=11,color='black'),
          legend.title = element_text(face="bold", color="black", size=14),
          legend.text = element_text(face="bold", color="black", size=12),
          strip.text.x = element_text(size=9, color="black", face="bold"),
          panel.background = element_blank()) +
    guides(color=guide_legend("Species"))
  dev.off()
  

  ########################
  png(paste('figure7','.png',sep=''), height = 7, width = 12, units = 'in', res=600)
  # boxplots.triple = boxplot(RE~species + samples.n + size, data = data.re, 
  #                           col=c('white','gray'), xaxt='n', ylim = c(0,10),
  #                           at = c(seq(1,6,1),  seq(8,13,1), seq(15,20,1), seq(22,27,1), seq(29,34,1)))
  # abline(v=7, col="gray");abline(v=14, col="gray");abline(v=21, col="gray");abline(v=28, col="gray")
  
  par(mfrow=c(2,5))
  par(mar = c(0, 0, 0, 0), oma = c(4, 6, 4, 0.5))
  par(tcl = -0.25)
  par(mgp = c(2, 0.6, 0))
  ylimit=c(-1,2)
  
  for(p in 1:5){
    temp = subset(data.re,size==p)
    boxplot(RE~species + samples.n, data = temp, col=c('white','gray'), xaxt='n', yaxt='n', ylim = c(0,10), at = c(1,2,4,5,7,8))
    
    if(p==1) {
      mtext('Root-mean-square-error (%)',side=2,line=2.2)
      axis(2,at=seq(0,10,2))
    }
    mtext(paste('Size class','',p),3)
    if(p==5) {
      legend(x=9.5,bty='n',xjust=-0.5,legend=c('Northern shrimp','Snow crab'),fill=c('white','gray'))
    }
  }
  
  for(p in 1:5){
    temp = subset(data.rb,size==p)
    boxplot(RE~species + samples.n, data = temp, col=c('white','gray'), xaxt='n', yaxt='n', ylim = c(-5,5), at = c(1,2,4,5,7,8))
    
    if(p==1) {
      mtext('Relative bias (%)',side=2,line=2.2)
      axis(2,at=seq(-4,4,1))
    }
    axis(1,at=c(1.5,4.5,7.5),labels=c(50,100,200))
    abline(h=0,col='gray')
  }

  mtext('Number of sampling locations',side=1,outer=T,line=2.2)
  dev.off()


