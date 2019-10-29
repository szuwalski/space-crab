
#mr_north = 0.1
#mr_south = 0.1
#mr_west  = 0.1
#mr_east  = 0.1
mr_dif = 0.4

# library(VAST);library(Matrix)
# strata.limits <- data.frame('STRATA'="All_areas")
# Region = "Eastern_Bering_Sea"
# Extrapolation_List = SpatialDeltaGLMM::Prepare_Extrapolation_Data_Fn( Region=Region, strata.limits=strata.limits )
# Data_Extrap = Extrapolation_List$Data_Extrap[Extrapolation_List$Data_Extrap$Area_in_survey_km2>0,]
# coords = cbind(Data_Extrap$Lon,Data_Extrap$Lat)
# coords = coords[coords[,1]<0,]

x=coords[,1] # x & y are unique
y=coords[,2]
grid_dataframe = data.frame(x=x,y=y)


DeltaX = range( diff(sort(grid_dataframe[,'x'])) )
DeltaY = range( diff(sort(grid_dataframe[,'y'])) )

n.loc = length(x)
# M_north = Matrix(0,nrow=n.loc,ncol=n.loc,sparse=TRUE)
# M_south = Matrix(0,nrow=n.loc,ncol=n.loc,sparse=TRUE)
# M_west = Matrix(0,nrow=n.loc,ncol=n.loc,sparse=TRUE)
# M_east = Matrix(0,nrow=n.loc,ncol=n.loc,sparse=TRUE)
Ins_Movement_m = Matrix(0,nrow=n.loc,ncol=n.loc,sparse=TRUE)

for (i in 1:n.loc){
  ind.x = which(abs(x-x[i])<3.5*DeltaX[2])
  ind.y = which(abs(y-y[i])<3.5*DeltaY[2])
  nb.ind = intersect(ind.x,ind.y)
  
  temp.x = x[nb.ind]-x[i]
  temp.y = y[nb.ind]-y[i]
  
  Ins_Movement_m[nb.ind[which(temp.x!=0)],i] = mr_dif
  
  # if (min(temp.x)!=0&max(temp.x)!=0&length(nb.ind)>=5){
  #   M_east[nb.ind[which(temp.x==max(temp.x))],i] = mr_east  #east
  #   M_west[nb.ind[which(temp.x==min(temp.x))],i] = mr_west  #west
  #   M_north[nb.ind[which(temp.y==max(temp.y))],i] = mr_north #north
  #   M_south[nb.ind[which(temp.y==min(temp.y))],i] = mr_south #south
  # }else{
  #   NULL
  # }
}

#Ins_Movement_m = M_east + M_west + M_north + M_south
diag(Ins_Movement_m) = -colSums(Ins_Movement_m)

temp = Matrix(0,nrow=n.loc,ncol=n.loc,sparse=TRUE)
diag(temp)=1

Movement_m = (temp + Ins_Movement_m/10) %*% (temp + Ins_Movement_m/10) %*% (temp + Ins_Movement_m/10) %*% (temp + Ins_Movement_m/10) %*% (temp + Ins_Movement_m/10) %*% 
             (temp + Ins_Movement_m/10) %*% (temp + Ins_Movement_m/10) %*% (temp + Ins_Movement_m/10) %*% (temp + Ins_Movement_m/10) %*% (temp + Ins_Movement_m/10)

#colSums(Movement_m)







