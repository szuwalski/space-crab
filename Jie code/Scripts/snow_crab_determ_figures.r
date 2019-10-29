
# figures 
# figure 1 side by side comparison of simulated and estimated distribution
setwd('/Users/jiecao/Desktop/UW_work/model_devp/SSST/SSST_output/Crab_example/crab_determ_move_true_N200_K100')
# density on the map

name = 'den_est' 
den.est = Save$Report$d_pkt

PlotDF1 = subset(MapDetails_List$PlotDF,Include==TRUE)

coords2 = cbind(PlotDF1$Lon,PlotDF1$Lat)
coords2 = coords2[coords2[,1]<0,]
P2 = SpatialPoints(coords2)
rast <- raster(ncol=200,nrow=100)
extent(rast) <- extent(coords2)
col=colorRampPalette(c("blue", "cyan", "yellow", "red"), bias=1)
#breakpoints <- c(0,0.1,0.5,1,1.5,2,3,5,10,50,150)
breakpoints <- c(1.2,1.5,2,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.8,4.1,4.4,4.7,5,8)

png(paste('figure1','.png',sep=''), height = 9, width = 12, units = 'in', res=600)

layout(mat = matrix(c(seq(1,10,1),rep(0,5),seq(11,20,1),rep(0,5),seq(21,30,1)), 
                    nrow = 5, 
                    ncol = 8),
       heights = rep(1,40),    # Heights of the two rows
       widths = c(1,1,0.1,1,1,0.1,1,1) )     # Widths of the two columns

layout.show(30)

par(mar = c(0, 0, 0, 0), oma = c(4, 4, 4, 1))
#par(tcl = -0.25)
#par(mgp = c(2, 0.6, 0))
for(p in c(1,3,5)){
  for(t in c(1,3,5,7,9)){
    P2$data = den.est[p,,t][PlotDF1$x2i]
    rast.temp <- rasterize(P2, rast, P2$data, fun = mean)
    plot(ocean,col="dark gray",axes=F,xlim=c(-180,-158),ylim=c(56,58))
    image(rast.temp,col=col(15),axes=TRUE,breaks=breakpoints,add=T,xlim=c(-180,-150),ylim=c(55,65))
    box()
    if(p==1){
      mtext(paste('Year ',t,sep=''),2,0.5,cex=0.9)
    }
    if(t==1){
      mtext(paste('Est'),3,1)
    }
    }
  for(t in c(1,3,5,7,9)){
    P$data = log(den.data[,p,t]/area)
    rast.temp <- rasterize(P, rast, P$data, fun = mean)
    plot(ocean,col="dark gray",axes=F,xlim=c(-180,-158),ylim=c(56,58))
    image(rast.temp,col=col(15),axes=TRUE,breaks=breakpoints,add=T,xlim=c(-180,-150),ylim=c(55,65))
    box()
    if(t==1){
      mtext(paste('Sim'),3,1)
    }
    if(t==9){
      if(p==1){
        mtext('Size class 1',1,1,adj=-0.6)
      }
      if(p==3){
        mtext('Size class 3',1,1,adj=-0.6)
      }
      if(p==5){
        mtext('Size class 5',1,1,adj=-0.6)
      }
    }
  }
}

dev.off()

# fishing mortality

F_est = exp(Save$Report$logF_male_kt)
breakpoints=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.5,2)

png(paste('F_est','.png',sep=''), height = 8, width = 10, units = 'in', res=600)

layout(mat = matrix(c(seq(1,10,1),rep(0,5),seq(11,20,1)), 
                    nrow = 5, 
                    ncol = 5),
       heights = rep(1,40),    # Heights of the two rows
       widths = c(1,1,0.3,1,1) )     # Widths of the two columns

layout.show(20)

par(mar = c(0, 0, 0, 0), oma = c(4, 4, 4, 1))

for (t in 1:5){
  P2$data = F_est[,t][PlotDF1$x2i]
  rast.temp <- rasterize(P2, rast, P2$data, fun = mean)
  plot(ocean,col="dark gray",axes=F,xlim=c(-180,-158),ylim=c(56,58))
  image(rast.temp,col=col(14),axes=TRUE,breaks=breakpoints,add=T,xlim=c(-180,-150),ylim=c(55,65))
  box()
  mtext(paste('Year ',t,sep=''),2,1,cex=0.9)
  if(t==1){
    mtext('Est',3,1)
  }
 }
  for(t in 1:5){
    P$data = F_kmeans[,t]
    rast.temp <- rasterize(P, rast, P$data, fun = mean)
    plot(ocean,col="dark gray",axes=F,xlim=c(-180,-158),ylim=c(56,58))
    image(rast.temp,col=col(14),axes=TRUE,breaks=breakpoints_f,add=T,xlim=c(-180,-150),ylim=c(55,65))
    box()
    if(t==1){
      mtext('Sim',3,1)
    }
  }
  
for (t in 6:10){
  P2$data = F_est[,t][PlotDF1$x2i]
  rast.temp <- rasterize(P2, rast, P2$data, fun = mean)
  plot(ocean,col="dark gray",axes=F,xlim=c(-180,-158),ylim=c(56,58))
  image(rast.temp,col=col(14),axes=TRUE,breaks=breakpoints,add=T,xlim=c(-180,-150),ylim=c(55,65))
  box()
  mtext(paste('Year ',t,sep=''),2,1,cex=0.9)
  if(t==6){
    mtext('Est',3,1)
  }
}
for(t in 6:10){
  P$data = F_kmeans[,t]
  rast.temp <- rasterize(P, rast, P$data, fun = mean)
  plot(ocean,col="dark gray",axes=F,xlim=c(-180,-158),ylim=c(56,58))
  image(rast.temp,col=col(14),axes=TRUE,breaks=breakpoints_f,add=T,xlim=c(-180,-150),ylim=c(55,65))
  box()
  if(t==6){
    mtext('Sim',3,1)
  }
}

dev.off()
































