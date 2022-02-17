###===============================###===============================###
### function for computing F-madogram
fun_fmado <- function(season){
  ### compute F-madogram for all runs
  ### load global event set (0 lag)
  mado_dist<-list()
  for(r in 1:length(vect_sim)){
    setwd(dir_data)
    load(file=paste(set_name,'_',vect_sim[r],"_global_events_",day_lag,".RData",sep="")) ### loads events_global_list
    events_global_list_all <- events_global_list
    
    events_global_list <- list()
    ### allow for seasonal F-madograms to be constructed
    ### extract seasonal data
    ### all seasons
    if(season=='all'){
      events_global_list <- events_global_list_all
    }
    if(season=='winter'){
      for(l in 1:length(events_global_list_all)){
        ### months
        months <- format(events_global_list_all[[l]]$Date,'%m')
        ### only events happening in winter
        events_global_list[[l]] <- events_global_list_all[[l]][which(months%in%c('12','01','02')),]
      }
    }
    if(season=='spring'){
      for(l in 1:length(events_global_list_all)){
      ### months
      months <- format(events_global_list_all[[l]]$Date,'%m')
      ### only events happening in winter
      events_global_list[[l]] <- events_global_list_all[[l]][which(months%in%c('03','04','05')),]
      }
    }
    if(season=='summer'){
      for(l in 1:length(events_global_list_all)){
      ### months
      months <- format(events_global_list_all[[l]]$Date,'%m')
      ### only events happening in winter
      events_global_list[[l]] <- events_global_list_all[[l]][which(months%in%c('06','07','08')),]
      }
    }
    if(season=='fall'){
      for(l in 1:length(events_global_list_all)){
      ### months
      months <- format(events_global_list_all[[l]]$Date,'%m')
      ### only events happening in winter
      events_global_list[[l]] <- events_global_list_all[[l]][which(months%in%c('09','10','11')),]
      }
    }
    ### spring and summer
    if(season=='spring_summer'){
      for(l in 1:length(events_global_list_all)){
      ### months
      months <- format(events_global_list_all[[l]]$Date,'%m')
      ### only events happening in winter
      events_global_list[[l]] <- events_global_list_all[[l]][which(months%in%c('03','04','05','06','07','08')),]
      }
    }
    
    ### remove columns representing temperature data (available only for a subset of stations)
    events_values <- events_dates <- list()
    for(i in 1:length(events_global_list)){
      ### extract only values, exclude all other information (will be needed for correlation analysis)
      events_values[[i]] <- events_global_list[[i]]$Q
      ### extract only dates
      events_dates[[i]] <- events_global_list[[i]]$Date
    }
    
    ### reformat data: put all the stations into the same data frame
    ### rows: events, columns: stations
    ### data: A matrix representing the data. Each column corresponds to one location.
    data <- do.call(cbind, events_values)
    
    ### coord: A matrix that gives the coordinates of each location. Each row corresponds to one location.
    coords <- cbind(coordo1[,1], coordo1[,2])
    
    ### compute F-madogram for pairs of stations: use as dissimilarity matrix
    ### package SpatialExtremes by Ribatet.
    ### The F-madogram is a rank-based distance measure. It is related to the extremal coefficient.
    if(length(data)>0){
      mado <- fmadogram(data,coords,plot=FALSE)
      ### extract F-madogram values
      mado_dist[[r]] <- mado[,2]
    }else{ ### set NA in the case of no events
      mado_dist[[r]] <- NA
    }
  }
  
  ### compute Euclidean distance (will later be replaced by river distance)
  #dist_mat_euc <- dist(cbind(camels_points@data$lon_cen,camels_points@data$lat_cen)[set_cascades,],upper=FALSE,diag=FALSE)
  dist_mat_euc <-dist(cbind(dfmeta$DX, dfmeta$DY))
  
  setwd(dir_analysis)
  # pdf('F_madogram_obs_vs_sim.pdf',width=10,height=4)
  # par(mfrow=c(1,2),mar=c(5,5,2,0))
  pdf(paste(season,'_F_madogram_obs_vs_sim.pdf',sep=''),width=5,height=4)
  par(mfrow=c(1,1),mar=c(5,5,2,2))
  for(r in 1:length(vect_sim)){
    ### euclidean distance
    if(r==1){
      # plot(dist_mat_river,mado_dist,col=col_Q,
      # ylab="F-madogram floods",xlab="River distance [X/Y]", ylim=c(0,0.25))
      ### fit a smoothing spline 
      spl_floods <- smooth.spline(dist_mat_euc,mado_dist[[r]],spar=0.95)
      plot(spl_floods,ylab="F-madogram Outaouais",
           xlab="Distance Euclidien  [X/Y en m]",type='l',col='black',ylim=c(0.05,0.15))
    }else{
      if(!is.na(mado_dist[[r]])){
        spl_floods <- smooth.spline(dist_mat_euc,mado_dist[[r]],spar=0.95)
        lines(spl_floods,col='orange')
      }
    }
  }
  #   for(r in 1:length(vect_sim)){
  #     ### river distance
  #     if(r==1){
  #       spl_floods <- smooth.spline(dist_mat_river[which(!is.na(dist_mat_river))],mado_dist[[r]][which(!is.na(dist_mat_river))],spar=0.95)
  #       plot(spl_floods,col=col_obs,lwd=2,type='l',ylab="F-madogram floods",
  #            xlab="River distance [X/Y]",ylim=c(0,0.16))
  #     }else{
  #       spl_floods <- smooth.spline(dist_mat_river[which(!is.na(dist_mat_river))],mado_dist[[r]][which(!is.na(dist_mat_river))],spar=0.95)
  #       lines(spl_floods,col=col_sim,lwd=2)
  #     }
  # }
  
  dev.off()
}

### plot F-madogram for certain season
fun_fmado(season='spring_summer')

### (g) extremal coefficient/ tail dependence chibar
### see spatialExtremes by Gilleland
###===============================###===============================###
### tail dependence: estimates tail dependence parameters chi and chibar
### chibar as introduced by Coles 1999
### estimators given by Reiss and Thomas 2007.

### write function for seasonal tail dependence computation
fun_taildep <- function(season){
  ### load data of the individual stations
  set_name<-'outaouais'
  data_stations_all <- list()
  data_stations <- list()
  for(l in station_numbers){
    setwd(dir_data)
    #load(file=paste(l,'_',set_name,"_stoch_sim_phase_rand_",'_wave_',n_wave,".RData",sep=""))
    #data_stations_all[[l]] <- data_stoch
    bassin<- BassinsOutaouais[l]
    load(file=paste(dir_data, bassin,"_stoch_sim_outaouais_Kappa_9997_LD",".Rdata",sep=""))
    data$simulation$timestamp<-NULL
    data_stations_all[[l]] <- data$simulation
    ### compute specific discharge to compare discharge at different sites
    data_stations_all[[l]][,4:length(data_stations_all[[l]])]<- data_stations_all[[l]][,4:length(data_stations_all[[l]])] *(1/(aire[l])) * 1000
    ### append a timestamp
    data_stations_all[[l]]$Date <-   paste(data_stations_all[[l]]$YYYY,"-",data_stations_all[[l]]$MM,"-",data_stations_all[[l]]$DD,sep="")
    data_stations_all[[l]]$Date <- as.POSIXct(strptime(data_stations_all[[l]]$Date,format="%Y-%m-%d",tz="GMT"))
    # ### remove first three rows
    data_stations_all[[l]] <- data_stations_all[[l]][,-c(1:3)]
    ### fill in missing dates if necessary
    if(length(is.na(data_stations_all[[l]]$Date))>0){
      ### create date sequence
      dates <- seq.Date(as.Date(data_stations_all[[l]]$Date[1]), as.Date(data_stations_all[[l]]$Date[length(data_stations_all[[l]]$Date)]), by=1)
      ### remove 29.2
      dates <- dates[format(dates, "%m %d") != "02 29"]
      data_stations_all[[l]]$Date <- dates
    }
    
    ### only seasonal data
    if(season=='all'){
      data_stations[[l]] <- data_stations_all[[l]]
    }
    ### winter
    if(season=='winter'){
      ### extract months
      months <- format(data_stations_all[[l]]$Date,'%m')
      data_stations[[l]] <- data_stations_all[[l]][which(months%in%c('12','01','02')),]
    }
    ### spring
    if(season=='spring'){
      ### extract months
      months <- format(data_stations_all[[l]]$Date,'%m')
      data_stations[[l]] <- data_stations_all[[l]][which(months%in%c('03','04','05')),]
    }
    ### summer
    if(season=='summer'){
      ### extract months
      months <- format(data_stations_all[[l]]$Date,'%m')
      data_stations[[l]] <- data_stations_all[[l]][which(months%in%c('06','07','08')),]
    }
    ### fall
    if(season=='fall'){
      ### extract months
      months <- format(data_stations_all[[l]]$Date,'%m')
      data_stations[[l]] <- data_stations_all[[l]][which(months%in%c('09','10','11')),]
    }
    ### spring and summer
    if(season=='spring_summer'){
      ### extract months
      months <- format(data_stations_all[[l]]$Date,'%m')
      data_stations[[l]] <- data_stations_all[[l]][which(months%in%c('03','04','05','06','07','08')),]
    }
  }
  ### remove empty list entries
  data_stations <- data_stations[lengths(data_stations) != 0]
  
  ### define threshold
  thresh <-0.95
  ### run for different thresholds
  for(t in c(0.8,0.95)){
    thresh <- t
    ### run through all simulations runs plus observations
    # for(r in 1:length(vect_sim)){
    for(r in 1:2){
      ### run tail dependence test for pairs of stations.
      #chi_mat <-chibar_mat <-tail_test_mat<- matrix(nrow = length(set_cascades),ncol=length(set_cascades))
      chi_mat <-chibar_mat <-tail_test_mat<- matrix(nrow = length(station_numbers),ncol=length(station_numbers))
      for(j in 1:length(station_numbers)){
        for(i in 1:length(station_numbers)){
          ### estimators given by Reiss and Thomas 2007.
          ### compute chi
          chi_mat[j,i] <- taildep(x=data_stations[[j]][,r],y=data_stations[[i]][,r], u=thresh,type='chi',na.rm=T)
          ### chi bar
          chibar_mat[j,i] <- taildep(x=data_stations[[j]][,r],y=data_stations[[i]][,r], u=thresh,na.rm=T,type='chibar')
          ### tail dependence test: Reiss and Thomas 2007
          tail_test_mat[j,i] <- taildep.test(x=cbind(data_stations[[j]][,r],data_stations[[i]][,r]),cthresh=-thresh,trans='relative.rank',div='n+1',na.action=na.omit)$p.value
        }
      }
      ### save results for each run
      setwd(dir_data)
      save(file=paste(season,'_',vect_sim[r],'_thresh_',thresh,'_tail_dep.RData',sep=''),chi_mat,chibar_mat,tail_test_mat)
    }
  }
  
  ### Note: H0 of tail dependence test is: random variables X and Y are dependent. 
  ### H0 rejected if p-value <0.05. I.e. if p-value smaller than 0.05 independence.
  ### determine colors
  #col_dep <- adjustcolor('#238443',0.5)
  #col_indep <- adjustcolor('grey',0.5)
  #col_obs_dep <- adjustcolor('black',0.5)
  #col_sim_dep <- adjustcolor(col_sim, 0.7)
# change color ldLD
  col_dep <- adjustcolor('blue',0.7)
  col_indep <- adjustcolor('red',0.3)
  col_obs_dep <- adjustcolor('black',0.5)
  col_sim_dep <- adjustcolor(col_sim, 0.7)
  ### run through two thresholds
  setwd(dir_analysis)
  # pdf('tail_dependence_obs_vs_sim.pdf',width=10,height=5)
  pdf(paste(season,'_tail_dependence_obs_vs_sim_chi.pdf',sep=''),width=10,height=5)
  
  par(mfrow=c(1,2))
  for(t in c(0.8,0.95)){
    thresh<-t
    ### look at results for observations
    setwd(dir_data)
    load(file=paste(season,'_',vect_sim[1],'_thresh_',thresh,'_tail_dep.RData',sep='')) ### loads chi_mat,chibar_mat,tail_test_mat)
    head(tail_test_mat)
    tail_test_mat_obs <- tail_test_mat
    chi_mat_obs <- chi_mat
    chibar_mat_obs <- chibar_mat
    ### compare to simulations
    load(file=paste(vect_sim[2],'_thresh_',thresh,'_tail_dep.RData',sep='')) ### loads chi_mat,chibar_mat,tail_test_mat)
    head(tail_test_mat)
    tail_test_mat_sim <- tail_test_mat
    chi_mat_sim <- chi_mat
    chibar_mat_sim <- chibar_mat
    
    ### determine agreement of obs vs. sim tail dependence
    ### focus on upper triangular matrix
    tail_test_mat_obs <- tail_test_mat_obs[upper.tri(tail_test_mat_obs)]
    tail_test_mat_sim <- tail_test_mat_sim[upper.tri(tail_test_mat_sim)]
    ### determine tail dependent pairs
    tail_dep_obs <- which(tail_test_mat_obs>0.05)
    tail_dep_sim <- which(tail_test_mat_sim>0.05)
    ### non-tail dependent pairs
    no_td_obs <- which(tail_test_mat_obs<0.05)
    no_td_sim <- which(tail_test_mat_sim<0.05)
    ### determine agreement/disagreement
    ### i) both simulations and observations indicate dependence
    dep <- tail_dep_obs[which(tail_dep_obs%in%tail_dep_sim)]
    ### ii) both simulations and observations indicate independence
    indep <- no_td_obs[which(no_td_obs%in%no_td_sim)]
    ### iii) observations dependence and simulations independence
    dep_obs <- tail_dep_obs[which(tail_dep_obs%in%no_td_sim)]
    ### iv) observations independence and simulations dependence
    dep_sim <- no_td_obs[which(no_td_obs%in%tail_dep_sim)]
    
    ### plot the chi estimator
    # plot(chi_mat_obs[upper.tri(chi_mat_obs)],chi_mat_sim[upper.tri(chi_mat_sim)],xlab='Observations',ylab='Simulations')
    
    ### plots for chi
    ### indicate pairs of stations where observations and simulations agree on independence
    plot(chi_mat_obs[upper.tri(chi_mat_obs)][indep],chi_mat_sim[upper.tri(chi_mat_sim)][indep],col=col_indep,
         xlab='Observations',ylab='Simulations',main=paste(thresh,'tail dependence chi',sep=' '),ylim=c(0,1),xlim=c(0,1))
    ### indicate pairs of stations where observations and simulations agree on dependence
    points(chi_mat_obs[upper.tri(chi_mat_obs)][dep],chi_mat_sim[upper.tri(chi_mat_sim)][dep],
           col=col_dep,xlab='Observations',ylab='Simulations',main='Agreement on dependence',ylim=c(0,1),xlim=c(0,1))
    ### pairs of stations where observations say dependence and simulations independence
    points(chi_mat_obs[upper.tri(chi_mat_obs)][dep_obs],chi_mat_sim[upper.tri(chi_mat_sim)][dep_obs],
           col=col_obs_dep,xlab='Observations',ylab='Simulations',main='Obs dep, sim indep',ylim=c(0,1),xlim=c(0,1))
    ### pairs of stations where simulations say dependence and observations independence
    points(chi_mat_obs[upper.tri(chi_mat_obs)][dep_sim],chi_mat_sim[upper.tri(chi_mat_sim)][dep_sim],
           col=col_sim_dep,xlab='Observations',ylab='Simulations',main='Sim dep, obs indep',ylim=c(0,1),xlim=c(0,1))
    ### add a legend
    legend('bottomright',legend=c('Independence','Dependence','Obs dependence','Sim dependence'),
           pch=1,col=c(col_indep,col_dep,col_obs_dep,col_sim_dep))
    abline(0,1)
    
    # ### plot for chibar
    # ### indicate pairs of stations where observations and simulations agree on independence
    # plot(chibar_mat_obs[upper.tri(chibar_mat_obs)][indep],chibar_mat_sim[upper.tri(chibar_mat_sim)][indep],col=col_indep,
    #      xlab='Observations',ylab='Simulations',main=paste(thresh,'tail dependence chibar',sep=''),ylim=c(-1,1),xlim=c(-1,1))
    # ### indicate pairs of stations where observations and simulations agree on dependence
    # points(chibar_mat_obs[upper.tri(chibar_mat_obs)][dep],chibar_mat_sim[upper.tri(chibar_mat_sim)][dep],
    #        col=col_dep,xlab='Observations',ylab='Simulations',main='Agreement on dependence',ylim=c(-1,1),xlim=c(-1,1))
    # ### pairs of stations where observations say dependence and simulations independence
    # points(chibar_mat_obs[upper.tri(chibar_mat_obs)][dep_obs],chibar_mat_sim[upper.tri(chibar_mat_sim)][dep_obs],
    #        col=col_obs_dep,xlab='Observations',ylab='Simulations',main='Obs dep, sim indep',ylim=c(-1,1),xlim=c(-1,1))
    # ### pairs of stations where simulations say dependence and observations independence
    # points(chibar_mat_obs[upper.tri(chibar_mat_obs)][dep_sim],chibar_mat_sim[upper.tri(chibar_mat_sim)][dep_sim],
    #        col=col_sim_dep,xlab='Observations',ylab='Simulations',main='Sim dep, obs indep',ylim=c(-1,1),xlim=c(-1,1))
    # ### add a legend
    # legend('bottomright',legend=c('Independence','Dependence','Obs dependence','Sim dependence'),
    #        pch=1,col=c(col_indep,col_dep,col_obs_dep,col_sim_dep))
    # abline(0,1)
  }
  
  dev.off()
}

### apply to winter
dir_data="H:/Projets_communs/2020/Outaouais PRsim/01_Intrants/AnalyseSpatial60Series"
fun_taildep(season='winter')
fun_taildep(season='spring_summer')
fun_fmado(season='winter')
fun_fmado(season='spring_summer')
fun_fmado(season='all')
fun_fmado(season='summer')
