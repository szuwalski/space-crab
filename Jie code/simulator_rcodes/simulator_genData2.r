
# generate survey data
library(Rlab)
n_loc_sample = n_loc_sample
#sampling_error=FALSE
if (sampling_error==TRUE){
  #ObsModel_sim = c("PL", "Lognormal")[1]
  obs = switch(ObsModel_sim,"PL"=1, "Lognormal"=2)
}
#logsigma_p = c(0.1,0.1,0.1,0.1,0.1)
#logsigma_p = c(0.25,0.25,0.25,0.25,0.25)

data_male = c()
data_female = c()
catch_data =c()

strata.limits <- data.frame('STRATA'="All_areas")
Region = "user"
load('/Users/jiecao/Google Drive/UW_work/snow_crab/simulator/shrimp_example/shrimp_grid.rda')
Extrapolation_List = SpatialDeltaGLMM::Prepare_Extrapolation_Data_Fn( Region=Region, strata.limits=strata.limits, input_grid = input_grid )
coords = cbind(Extrapolation_List$Data_Extrap$Lon,Extrapolation_List$Data_Extrap$Lat)
coords = coords[coords[,1]<0,]

area = Extrapolation_List$Area_km2_x

Dist = stats::dist(loc_x)
range = 10*min(Dist)

for (i in 1:n_t){
  
  # pb <- txtProgressBar(min = 0, max = n_t-1, style = 3, char = "*")
  # Sys.sleep(0.1)
  # setTxtProgressBar(pb, i)
  
  sample_id = sample(1:n_s, n_loc_sample, replace=FALSE)
  min_dist = min(stats::dist(loc_x[sample_id,]))
  
  while(min_dist<range){
    sample_id = sample(1:n_s, n_loc_sample, replace=FALSE)
    min_dist = min(stats::dist(loc_x[sample_id,]))
  }
  
  for (p in 1:n_p){
    temp_male = N_at_size_ns[,,i][,p]

    for (k in 1:n_loc_sample){
      temp_n = temp_male[sample_id[k]]
      if(sampling_error){
        if(obs==2){
          temp_n=rlnorm(1,log(temp_n)-(logsigma_p[p]^2)/2,logsigma_p[p])
          #print(temp_n)
        }
        if(obs==1){
          encout = 1-exp(-temp_n)
          bp = rbern(1,encout)
          #print(bp)
            if(bp==0){
              temp_n = 0
            }else{
              temp_n = rlnorm(1,log(temp_n)-(logsigma_p[p]^2)/2,logsigma_p[p])
            }
        }
      }
      temp1 = c(p,i,temp_n,coords[sample_id[k],],1)
      #temp1 = c(p,i,temp_n,coords[sample_id[k],],area[sample_id[k]])
      #temp1 = c(p,i,round(temp_male[sample_id[k]]/area,0),coords[sample_id[k],],area)
      data_male = rbind(data_male,temp1)
    }
  }
}

N_at_size_ns_tot=array (data = 0, dim = c(n_s,n_p,n_t))
index = matrix(0,nrow=n_t,ncol=n_p)
for(t in 1:n_t){
  for(p in 1:n_p){
    for(s in 1:n_s){
      N_at_size_ns_tot[s,p,t]=N_at_size_ns[s,p,t]*area[s]
    }
    index[t,p] = sum(N_at_size_ns_tot[,p,t])
  }
}

#survey_data = rbind(data_male,data_female)
#if(R_sex_ratio == 1) 
survey_data = data_male
colnames(survey_data)=c('size_class','year','count','lon','lat','area_swept')

write.csv(survey_data,file = '/Users/jiecao/Desktop/UW_work/model_devp/SSST/Simulation_Data/shrimp/survey_data.csv',row.names = FALSE)

catch_data_temp = C_at_size_ns

for (i in 1:n_t){
  catch_data = rbind(catch_data,catch_data_temp[,,i])
}
year = rep(seq(1,n_t,1),each=n_s)
lat = rep(coords[,2],n_t)
lon = rep(coords[,1],n_t)
catch_data = cbind(lat,lon,year,catch_data)

if(sampling_error){
  for (i in 1:nrow(catch_data)){
    for (p in 1:n_p){
      catch_data[i,(3+p)] = rlnorm(1,log(catch_data[i,(3+p)])-logcatch_p[p]^2/2,logcatch_p[p])
    }
  }
}

write.csv(catch_data,file = '/Users/jiecao/Desktop/UW_work/model_devp/SSST/Simulation_Data/shrimp/catch_data.csv')

survey_data = read.csv('/Users/jiecao/Desktop/UW_work/model_devp/SSST/Simulation_Data/shrimp/survey_data.csv')
catch_data = read.csv('/Users/jiecao/Desktop/UW_work/model_devp/SSST/Simulation_Data/shrimp/catch_data.csv')

# data for non-spatial assessment model
if (runNonspatial == TRUE){

tot_c = apply(catch_data[,5:(4+n_p)],1,sum)
catch_data = cbind(catch_data,tot_c)

survey_data_nonspatial = matrix(NA,nrow=n_t,ncol=(n_p+6))
catch_data_nonspatial = matrix(NA,nrow=n_t,ncol=(n_p+9))
#nsize_at_loc_temp = matrix(NA,nrow=n_loc_sample,ncol=n_p)
temp3=c();temp5=c();temp_tot_catch=c();temp_tot_catch_s=c()
sample.comp = matrix(NA,nrow=n.eff,ncol=n_p)
sample.comp = array(NA, dim = c(n.eff,n_p+2,n_t))
catch_tot_error_n = c()

for (i in 1:n_t){
  temp2=subset(survey_data,year==i)
  temp4=subset(catch_data,year==i)
      for(p in 1:n_p){
        temp3[p]=sum(subset(temp2,size_class==p)$count)/(n_loc_sample*1000)
        temp5[p]=sum(temp4[,4+p])/1000
      }
  
  # if (sampling_error==TRUE){
  #   for (n in 1:nrow(temp4)){
  #     temp_tot_catch_s[n] = rlnorm(1,log(temp4[n,4+n_p+1]/1000)-logcatch_p^2,logcatch_p)
  #   }
  #   catch_tot_error_n=c(catch_tot_error_n,temp_tot_catch_s)
  #   temp_tot_catch[i] = sum(temp_tot_catch_s)
  #   sample.id = sample(1:nrow(temp4),n.eff,replace=F)
  #   sample.comp[,1+n_p,i] = sample.id
  #   temp.comp = cbind(temp4[sample.id,],temp_tot_catch_s[sample.id]*1000)
  #   for (nn in 1:nrow(temp.comp)){
  #     sample.comp[nn,(1:n_p),i]=rmultinom(1, temp.comp[nn,4+n_p+2], temp.comp[nn,5:(4+n_p)]/sum(temp.comp[nn,5:(4+n_p)]))
  #   }
  #   temp5 = apply(sample.comp[,(1:n_p),i],2,sum)/sum(sample.comp[,(1:n_p),i])
  #   
  # }else{
  #   temp_tot_catch[i] = sum(temp4[,4+n_p+1])/1000
  # }
  temp_tot_catch[i] = sum(temp4[,4+n_p+1])/1000
  survey_data_nonspatial[i,]=c(i,1,0,sum(temp2$count)/n_loc_sample,sd(temp2$count)/mean(temp2$count),n_loc_sample,temp3)
  catch_data_nonspatial[i,]=c(i,1,1,temp_tot_catch[i],sd(temp4[,4+n_p+1])/mean(temp4[,4+n_p+1]),-1,1,-1,n.eff,temp5)
}

survey.data = survey_data_nonspatial
catch.data = catch_data_nonspatial

# if (sampling_error)
#   catch_data_error = cbind(catch_data,catch_tot_error_n*1000)
}





