### extract binary matrix of flood events
### plot spatial flood events: for each event (simulated within the 671 joint simulations), indicate on map, which stations were affected
### extract above threshold events
#ldLD
### wavelet analysis
library(WaveletComp)
### shapefiles
#library(rgdal)
### data.frame extensions
library(data.table)
### f-madogram
library(SpatialExtremes)
### tail dependence
library(extRemes)
# lecture excel pour madogram
library("xlsx")
# Cleaning memory
rm(list=ls())
gc()



### define directories
#repbase<-"H:\\Projets_communs\\2020\\Outaouais PRSIM\\01_Intrants\\Example Rdata pour Leslie\\"
repbase<-"H:\\Projets_communs\\2020\\Outaouais PRsim\\01_Intrants\\AnalyseSpatial60Series\\"

repdata<-"H:\\Projets_communs\\2020\\Outaouais PRSIM\\02_Calculs\\madogram\\"

excelfilecoord <- stringr::str_c(repdata,"HSAMI_centroides.xlsx")
dfmeta <- read.xlsx(excelfilecoord, "BassinsOutaouais",as.data.frame=TRUE, header=TRUE)

set.seed(1)
alt <- dfmeta$ALTITUDE
lon <- (dfmeta$DX)
lat <- dfmeta$DY
aire<-dfmeta$SUPERFICIE
coordo <- as.matrix(data.frame(lon=lon,lat=lat, alt=alt, aire=aire))
coordo1 <- coordo




#dir_main <- "~/Documents/Travail de MB/" ### path needs to be pointing to the hydro_Quebec folder I sent you.
dir_main <- "H:/Projets_communs/2020/Outaouais PRsim/02_Calculs/madogram/"

dir_data <- repbase
#ldLD
dir_stoch_sim_multi_site <- dir_data
dir_analysis <- paste(dir_main,"results",sep='')

### Define colors for plotting
###===============================###===============================###
col_sim <- adjustcolor('#d95f0e',0.3)
col_sim_tran <- adjustcolor("#d95f0e",alpha=0.2)
col_obs <- adjustcolor( "black", alpha.f = 0.2)
### define colors for plotting wavelets
col_wave <- c('#a6cee3','#b2df8a','#fb9a99','#fdbf6f')
### darker for observations
col_wave_obs <- c('#3182bd','#31a354','#de2d26')

### define helper functions
### function of identifying peaks: local maxima with m=15 points on either side being smaller than it
find_peaks <- function (x, m = 15){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks
}



n_series = 60
#set_cascades<-c(668, 661,590,608)
Bassins =c("Bark Lake",                    "Baskatong",                    "Cabonga",                     
           "Carillon et Hull",             "Chelsea",                      "Chenaux et Noire",            
           "Chute-des-Chats" ,             "Des Joachims" ,                "Dozois" ,                     
           "High Falls",                   "Kamaniskeg",                   "Kiamika" ,                    
           "Kipawa",                       "Lac des Quinze",               "Lac du poisson blanc",        
           "Lac Temiscamingue a Angliers", "Lac Victoria et lac Granet",   "Lady Evelyn",                 
           "Lower Notch et Indian Chute",  "Madawaska-Arnprior",           "Maniwaki",                    
           "Masson",                       "Mistinikon",                   "Mitchinamecus",               
           "Mont-Laurier",                 "Mountain Chute",               "Otto Holden",                 
           "Paugan",                       "Rabbit Lake",                  "Rapide-2" ,                   
           "Rapide-7",                     "Riviere Blanche",              "Riviere Bonnechere",          
           "Riviere Coulonge",             "Riviere Dumoine",              "Riviere Kinojevis",           
           "Riviere Mattawa",              "Riviere Mississippi",          "Riviere Petawawa" ,           
           "Riviere Petite Nation",        "Riviere Rideau" ,              "Riviere Rouge",               
            "Riviere South Nation",         "Du Loup" ,                     "Du Nord",                     
           "LAssomption" ,                 "Maskinonge" ,                  "Nicolet" ,                    
           "Richelieu",                    "Riviere Saint Francois",       "Yamaska")                
  
BassinsOutaouais =c("Bark Lake",                    "Baskatong",                    "Cabonga",                     
           "Carillon et Hull",             "Chelsea",                      "Chenaux et Noire",            
           "Chute-des-Chats" ,             "Des Joachims" ,                "Dozois" ,                     
           "High Falls",                   "Kamaniskeg",                   "Kiamika" ,                    
           "Kipawa",                       "Lac des Quinze",               "Lac du poisson blanc",        
           "Lac Temiscamingue a Angliers", "Lac Victoria et lac Granet",   "Lady Evelyn",                 
           "Lower Notch et Indian Chute",  "Madawaska-Arnprior",           "Maniwaki",                    
           "Masson",                       "Mistinikon",                   "Mitchinamecus",               
           "Mont-Laurier",                 "Mountain Chute",               "Otto Holden",                 
           "Paugan",                       "Rabbit Lake",                  "Rapide-2" ,                   
           "Rapide-7",                     "Riviere Blanche",              "Riviere Bonnechere",          
           "Riviere Coulonge",             "Riviere Dumoine",              "Riviere Kinojevis",           
           "Riviere Mattawa",              "Riviere Mississippi",          "Riviere Petawawa" ,           
           "Riviere Petite Nation",        "Riviere Rideau" ,              "Riviere Rouge",               
           "Riviere South Nation")        

set_name="Outaouais"

#nb = length(Bassins)
nb = length(BassinsOutaouais)
# #only once
# #fichstochsim<-paste(repbase,"stoch_sim_10_outaouais_Kappa_5_9997_LD.Rdata",sep='')
# fichstochsim<-paste("H:\\Projets_communs\\2020\\Outaouais PRsim\\01_Intrants\\Example Rdata pour Leslie\\0.9997\\","stoch_sim_r60.Rdata",sep='')
# load(fichstochsim)
# for(b in 1: nb) {
#    noms<-names(stoch_sim_concat_r60)
#    bassin=noms[b]
#    print(bassin)
#    data<-stoch_sim_concat_r60[[b]]
#    #plot(data$simulation$timestamp,  data$simulation$r20)
#    nomfich=paste(dir_data, bassin, "_stoch_sim_outaouais_Kappa_9997_LD",".Rdata",sep="")
#    save(file=nomfich,data)
#  }


station_numbers= 1:length(BassinsOutaouais)

fun_joint_floods_extraction <- function(station_numbers,set_name,n_series){
  setwd(dir_data)
  data_stations <- list()
  
  for(l in station_numbers){
    ### load data of the individual stations
    #load(file=paste(l,'_',set_name,"_stoch_sim_phase_rand_",'_wave_',n_wave,".RData",sep=""))
    #POUR OUTAOUAIS
    bassin<- BassinsOutaouais[l]
    load(file=paste(dir_data, bassin,"_stoch_sim_outaouais_Kappa_9997_LD",".Rdata",sep=""))
    data$simulation$timestamp<-NULL
    data_stations[[l]] <- data$simulation
    #data_stations$Nom<- bassin
    ### compute specific discharge to compare discharge at different sites
    data_stations[[l]][,4:length(data_stations[[l]])]<- data_stations[[l]][,4:length(data_stations[[l]])] *(1/(aire[l])) * 1000
    ### append a timestamp
    data_stations[[l]]$timestamp <-   paste(data_stations[[l]]$YYYY,"-",data_stations[[l]]$MM,"-",data_stations[[l]]$DD,sep="")
    data_stations[[l]]$timestamp <- as.POSIXct(strptime(data_stations[[l]]$timestamp,format="%Y-%m-%d",tz="GMT"))
  }
  ### remove empty list entries
  data_stations <- data_stations[lengths(data_stations) != 0]
  
  ### define simulation of interest
  
  vect_sim <- c('Qobs',paste('r',1:n_series,sep=''))
  events_run <- list()
  ### run through the different simulation runs
  for(r in 1:length(vect_sim)){
    bt_ids<-events_station <- list()
    event_dates <- list()
    ### extraction of below threshold values
    for(l in 1:length(station_numbers)){
      plot(data_stations[[l]]$timestamp,data_stations[[l]]$Qobs,typ='l')
      ### determine threshold: annual maxima based
      ### compute AMs values
      ams <- data.frame()
      # ams <- do.call("rbind", by(data_stations[[i]][,which(names(data_stations[[i]])==vect_sim[r])], data_stations[[i]]$YYYY,
      #                            function(x) x[which.max(x)]))
      ams <- aggregate(data_stations[[l]][,which(names(data_stations[[l]])==vect_sim[r])],by=list(data_stations[[l]]$YYYY) ,max)
      ### compute mean AMs
      # threshold <- mean(ams$Flow) ### produces only very few events
      # threshold <- quantile(ams$Flow,0.25)
      threshold <- quantile(ams$x,0.25,na.rm=T)
      ### plot threshold
      abline(h=threshold,col="grey",lty=2)
      ### define below threshold values
      # if(upper==FALSE){
      #   # bt_values[[i]] <- data_in[,i][which(data_in[,i]<quant_ex[[i]])]
      #   bt_ids[[i]]<-which(data_in[,i]<quant_ex[[i]])
      # }else{
      ### extract above threshold values
      # bt_values[[i]] <- data_in[,i][which(data_in[,i]>quant_ex[[i]])]
      bt_ids[[l]]<-which(data_stations[[l]][,which(names(data_stations[[l]])==vect_sim[r])]>threshold)
      ### identify corresponding dates
      event_dates[[l]] <- data_stations[[l]]$timestamp[bt_ids[[l]]]
      ### determine values over threshold
      pot <- data.frame('Date'=event_dates[[l]],'Flow'=data_stations[[l]][,which(names(data_stations[[l]])==vect_sim[r])][bt_ids[[l]]])
      
      ### choose only one event per month instead
      pot$timestamp <- format(pot$Date,format="%Y %m")
      
      ### go out of station if there is no event (ts too short)
      if(length(pot$Date)==0){
        return()
      }
      
      ### choose unique timestamps
      events <- unique(pot$timestamp)
      ### assign event ids
      pot$event_id <- NA
      event_df <- NA
      for(i in 1:length(events)){
        pot$event_id[which(pot$timestamp%in%events[i])] <- i
        ### determine peak discharge of individual events
        max_d <- max(pot$Flow[pot$event_id==i])
        ### set up data frame with time of occurrence and maximum dischcharge
        event <- pot[pot$event_id==i,][which.max(pot$Flow[pot$event_id==i]),]
        ### append to dataframe
        event_df <- rbind(event_df,event)
      }
      ### omit first row
      event_df <- event_df[-1,]
      event_df_orig <- event_df
      
      ### remove non-independent events: prescribe a minimum time lag between independent events
      time_lag <- 10
      ### run through all events
      for(i in 1:length(row.names(event_df))){
        ### determine date of occurrence of event i
        date_occ <- event_df[i,]$Date
        ### compute time difference to other events and determine thos events closer than time lag
        pos_close <- which(abs(event_df$Date-date_occ)/(60*60*24)<time_lag)
        ### if an event is not separated by the prescribed time lag, remove the smaller of the two
        if(length(pos_close)>1){
          for(p in 2:length(pos_close)){
            if(event_df[i,]$Flow>event_df[pos_close[p],]$Flow)
            {
              ### remove event post_close[p]
              event_df[pos_close[p],]<-NA
            }else{
              ### remove tested event
              event_df[i,]<-NA
            }
          }
        }else{
          event_df <-event_df ### if event is independent already
        }
      }
      
      ### remove dependent events
      event_df <- event_df[apply(event_df, 1, function(y) !all(is.na(y))),]
      ### reassign event ids
      event_df$event_id <- seq(1:length(event_df$event_id))
      ### add to plot
      points(event_df$Date,event_df$Flow,col="red")
      events_station[[l]] <- event_df
      print(paste(l,'done',sep=' '))
      ### go to the next station
    }
    ### set up joint event set 
    ### extract data for individual stations 
    event_dates_ind <- event_mag_ind <- list()
    for(l in 1:length(station_numbers)){
      event_dates_ind[[l]] <- events_station[[l]]$Date
      ### determine event magnitudes
      event_mag_ind[[l]] <- events_station[[l]]$Flow
    }
    
    ### combine pot events of all stations
    pot_all <- do.call("c",event_dates_ind)
    pot_all <- as.POSIXct(strptime(pot_all,format="%Y-%m-%d",tz="GMT"))
    ### convert to week of occurrence
    # pot_all <- format(pot_all,format="%Y %W")
    
    ### identify unique events according to day of occurrence
    events <- unique(pot_all)
    ### sort events by date of occurrence
    events <- sort(events)
    
    ### total number of events
    length(events)
    
    ### look at overall events
    ###===============================###===============================###
    ### count number of stations affected by a certain event
    ### allow for a short time lag between event occurrences at different stations
    day_lag <- 2
    
    ### run through all events identified over all stations
    occ_count <-rep(list(rep(list(NA),times=length(station_numbers))),times=length(events))
    stations_occ_event<-list()
    for(i in 1:length(events)){
      ### choose event of interest
      event <- events[i]
      ### determine dates +/- day_lag around this event
      start <- event-(day_lag*60*60*24)
      end <- event+(day_lag*60*60*24)
      ### event dates to search for
      dates <- seq(start,end,by=60*60*24)
      
      ### run through all stations and check whether the same event was also identified at the station under consideration
      
      for(l in 1:length(station_numbers)){
        ### distinguis between occurrence and non-occurrence
        occ_count[[i]][l]<- length(event_dates_ind[[l]][which(event_dates_ind[[l]]%in%dates)])
      }
      ### count stations with event occurrence
      ### store stations of co-occurrence for each event
      stations_occ_event[[i]] <- which(unlist(occ_count[[i]])>0)
    }
    
    ### count stations of occurrence per event
    number_occ<-list()
    for(i in 1:length(events)){
      number_occ[[i]] <- length(stations_occ_event[[i]])
    }
    number_occ <- unlist(number_occ)
    
    ### remove non-independent events
    ### remove entries of 'double' counts, i.e. events listed twice because they were slightly shifted in space (e.g. two successive days)
    ### For events happening within the same time window (7 days), only retain the date when most stations were affected
    time_lag <- 7
    ### for each event, compute time difference to other events
    ### if there is an event date closer than seven days, remove the event, which affected less stations
    ### test each individual event
    for(i in 1:length(events)){
      pos_close <- which((abs(events[i]-events)/(60*60*24))<time_lag)
      events[pos_close]
      ### compare number of stations affected
      ### remove tested event if it is not the one affecting the most stations
      larger_events <- which(number_occ[pos_close]>number_occ[i])
      ### if at least one event has a larger extent, remove tested event
      if(length(larger_events)>0){
        events[i]<-NA ### set to NA
      }
    }
    
    ### there are still some dependent events in case the same number of stations was affected on two successive days.
    ### in this case, remove the second date
    for(i in 1:length(events)){
      pos_close <- which((abs(events[i]-events)/(60*60*24))<time_lag)
      ### if there are still non-independent events,
      if(length(events[pos_close])>1){
        ### remove the second event
        events[i+1]<-NA
      }
    }
    
    ### remove the non-independent events
    pos_dep <- which(is.na(events))
    events <- events[-pos_dep]
    length(events)
    number_occ <- number_occ[-pos_dep]
    stations_occ_event <- stations_occ_event [-pos_dep]
    ### set up a data frame where each overall events gets an entry.
    ### set stations with non-occurrence to NA
    ### data frame with one event per row
    df_all_events <- data.frame("event"=events)
    row.names(df_all_events) <- events
    ### for each station, add relevant events
    for(l in 1:length(station_numbers)){
      ### create new row
      df_all_events[,l] <-events
      ### run through alls events
      for(i in 1:length(events)){
        event <- events[i]
        day_lag <- 2
        ### determine dates +/- day_lag around this event
        start <- event-(day_lag*60*60*24)
        end <- event+(day_lag*60*60*24)
        ### event dates to search for
        dates <- seq(start,end,by=60*60*24)
        ### does this event occur in this catchment?
        occ <- length(which(event_dates_ind[[l]]%in%dates))
        ### if it does not,
        if(occ==0){
          ### set data frame entry to NA
          df_all_events[,l][i]<-NA
        }else{
          ### if it does, add date to data frame
          df_all_events[,l][i]<-event_dates_ind[[l]][which(event_dates_ind[[l]]%in%dates)]
        }
      }
    }
    
    ### attribute names to columns
    names(df_all_events)<-station_numbers
    ### save as RData
    setwd(dir_analysis)
    save(file=paste(set_name,'_',vect_sim[r],'_df_event_dates_all_stations.Rdata',sep=''),df_all_events)
    print(paste('run',r,'done',sep=' '))
    ### to to next run
  }
  
}

### apply to whole data set
fun_joint_floods_extraction(station_numbers=station_numbers,set_name='outaouais',n_series=60)

#===========================================================================================================================
### requires joint event extraction

### define simulation of interest
vect_sim <- c('Qobs',paste('r',1:n_series,sep=''))
fun_joint_event_extraction<- function(day_lag,n_wave,set_name){
  for(r in 1:length(vect_sim)){
    setwd(dir_analysis)
    ### load event set corresponding to run r
    load(file=paste(set_name,'_',vect_sim[r],'_df_event_dates_all_stations.Rdata',sep=''))### loads df_all_events
    events <- row.names(df_all_events) ### currently refers to lag 7
    events <- as.POSIXct(strptime(events,format="%Y-%m-%d",tz="GMT"))
    
    events_list <- list()
    data_stations<- list()
    for(l in 1:nb){
      ### load data of the individual stations
      setwd(dir_data)
      #load(file=paste(l,'_',set_name,"_stoch_sim_phase_rand_",'_wave_',n_wave,".RData",sep=""))
      #data_stations[[l]] <- data_stoch
      
      bassin<- BassinsOutaouais[l]
      load(file=paste(dir_data, bassin,"_stoch_sim_outaouais_Kappa_9997_LD",".Rdata",sep=""))
      data$simulation$timestamp<-NULL
      data_stations[[l]] <- data$simulation
      
      
      ### compute specific discharge to compare discharge at different sites
      data_stations[[l]][,4:length(data_stations[[l]])]<- data_stations[[l]][,4:length(data_stations[[l]])] *(1/(aire[l])) * 1000
      ### append a timestamp
      data_stations[[l]]$Date <-   paste(data_stations[[l]]$YYYY,"-",data_stations[[l]]$MM,"-",data_stations[[l]]$DD,sep="")
      data_stations[[l]]$Date <- as.POSIXct(strptime(data_stations[[l]]$Date,format="%Y-%m-%d",tz="GMT"))
      data <- data_stations[[l]]
      ### remove first three rows
      data <- data[,-c(1:3)]
      ### plot time series
      plot(data$Date,data[,r],ylab = expression(bold(paste("Discharge [m"^3,"/s]"))),
           xlab=expression(bold("Time [y]")),type="l")
      
      ### for each overall event, determine station maxima =/- x days around the event
      
      events_station <- list()
      for(i in 1:length(events)){
        start <- which(data$Date==events[i])-day_lag ### define start of search window
        end <- which(data$Date==events[i])+day_lag ### define end of search window
        if(length(start)==0){
          data_max <- data.frame(data[1,])
        }else{
          if(start<1){
            start<-1
          }
          ### determine values within this window
          data_i <- data[start:end,]
          if(is.na(data_i[,r])){
            data_max <- data_i[1,]
          }else{
            ### determine max within this window
            data_max <- data_i[which.max(data_i[,r]),]
          }
        }
        
        points(data_max$Date,data_max[,r],col="red")
        events_station[[i]] <- data.frame('Date'=data_max$Date,'Q'=data_max[,r])
      }
      ### combine events to datafrmae
      # events_station <- do.call("rbind.data.frame", events_station)
      events_station <- data.frame(rbindlist(events_station, fill = TRUE))
      ### store in list
      events_list[[l]] <- events_station
      print(paste('set',r,'station',l,'done',sep=' '))
    }
    
    ### save events_list
    setwd(dir_analysis)
    save(events_list,file=paste(set_name,'_',vect_sim[r],"_joint_events_lag_",day_lag,".RData",sep=""))
    print(paste(r,'done',sep=' '))
  }
}

### apply for different time lags
### 0: exact date of occurrence as for original station
#--------------------ATT TIME CONSUMING-----------------------------
fun_joint_event_extraction(day_lag=0,n_wave=10,set_name='outaouais')


### (iii) identify events relevant at a regional scale
### global event extraction
###===============================###===============================###
day_lag<-0
fun_global_events <- function(day_lag,set_name){
  setwd(dir_analysis)
  for(r in 1:length(vect_sim)){
    ### load event set corresponding to run r
    load(file=paste(set_name,'_',vect_sim[r],'_df_event_dates_all_stations.Rdata',sep=''))### loads df_all_events
    
    events <- row.names(df_all_events) ### currently refers to lag 7
    events <- as.POSIXct(strptime(events,format="%Y-%m-%d",tz="GMT"))
    
    ### load joint events
    setwd(dir_analysis)
    load(file=paste(set_name,'_',vect_sim[r],"_joint_events_lag_",day_lag,".RData",sep=""))
    ### remove empty list entries
    events_list <- events_list[lengths(events_list) != 0]
    ### set up data frame
    df.ranks <- data.frame(matrix(data=NA,ncol=length(station_numbers),nrow=length(events_list[[2]]$Q)))
    ### extract ranks
    for(l in 1:length(station_numbers)){
      df.ranks[l] = rank(events_list[[l]]$Q,na.last="keep")
    }
    colnames(df.ranks) <- station_numbers
    # row.names(df.ranks) <- format(events_list[[2]]$Date,format="%Y %m %d")
    row.names(df.ranks) <- events
    
    ### global events: identify those events where all stations show high ranks
    ### compute standard deviation over ranks: if small: global event
    sd_ranks <- apply(df.ranks,MARGIN=1,FUN=sd,na.rm=T)
    ### compute rank sum over all stations
    sum_ranks <- apply(df.ranks,MARGIN=1,FUN=sum,na.rm=T)
    #sd(df.ranks[1,],na.rm=T)
    plot(sort(sd_ranks))
    ### determine threshold between local and global events
    # abline(h=quantile(sort(sd_ranks),0.25)) ### everyting below threshold would be global
    abline(h=quantile(sort(sd_ranks),0.25)) ### everyting below threshold would be global
    
    plot(sort(sum_ranks))
    ### determine threshold
    # abline(h=quantile(sort(sum_ranks),0.25))
    abline(h=quantile(sort(sum_ranks),0.25))
    
    ### define global events according to rank variability
    # global_events_ranks <- df.ranks[which(sd_ranks<quantile(sort(sd_ranks),0.25)),]
    global_events_ranks <- df.ranks[which(sd_ranks<quantile(sort(sd_ranks),0.25)),]
    
    ### define global events according to rank sum across all stations
    # global_events_sum <- df.ranks[which(sum_ranks>quantile(sort(sum_ranks),0.25)),]
    global_events_sum <- df.ranks[which(sum_ranks>quantile(sort(sum_ranks),0.25)),]
    
    ### schnittmenge of the two
    ### define date of occurrence of global events
    global_dates <- row.names(global_events_ranks)[which(row.names(global_events_ranks)%in%row.names(global_events_sum))]
    global_dates_formatted <- as.POSIXct(global_dates,format="%Y-%m-%d",tz="GMT")
    
    ### compare these dates to dates chosen when looking at number of stations where global event was detected
    #event_dates_global ### extracted via counting number of stations where event
    ### was detected
    global_dates ### extracted my strategy (after having composed synchronised event data set)
    
    ### extract globel events data
    events_global_list <- list()
    pos_global <- list()
    for(l in 1:length(station_numbers)){
      for(i in 1:length(global_dates)){
        ### try to identify exact correspondence
        if(length(which(format(events_list[[l]]$Date,format="%Y %m %d")==format(global_dates_formatted[i],format="%Y %m %d")))>0){
          pos_global[i] <- which(format(events_list[[l]]$Date,format="%Y %m %d")==format(global_dates_formatted[i],format="%Y %m %d"))
        } else {
          ### dates of occurrence might be shifted by one day
          ### test coincidence of week of occurrence
          pos_global[i] <- which(format(events_list[[l]]$Date,format="%Y %m")==format(global_dates_formatted[i],format="%Y %m"))[1]
        }
      }
      pos_global <- unlist(pos_global)
      #pos_global <- which(format(events_list[[l]]$Date,format="%Y %m %d")%in%format(global_dates_formatted,format="%Y %m %d"))
      events_global_list[[l]]  <- events_list[[l]][pos_global,]
      print(paste(l,'done',sep=' '))
    }
    
    ### save global event set
    setwd(dir_analysis)
    save(file=paste(set_name,'_',vect_sim[r],"_global_events_",day_lag,".RData",sep=""),events_global_list)
  }
}
fun_global_events(day_lag=0,set_name='outaouais')

### compute F-madogram for different simulation runs
###===============================###===============================###

### compute F-madogram for all runs
### load global event set (0 lag)
vect_sim <- c('Qobs',paste('r',1:n_series,sep=''))
set_name='outaouais'
day_lag<-0
mado_dist<-list()
for(r in 1:length(vect_sim)){
  setwd(dir_analysis)
  load(file=paste(set_name,'_',vect_sim[r],"_global_events_",day_lag,".RData",sep="")) ### loads events_global_list
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
  mdata = data
  mado <- fmadogram(mdata,coords,plot=FALSE)
  ### extract F-madogram values
  mado_dist[[r]] <- mado[,2]
}
setwd(dir_stoch_sim_multi_site)


### compute Euclidean distance (will later be replaced by river distance)
#dist_mat_euc <- dist(cbind(camels_points@data$lon_cen,camels_points@data$lat_cen)[set_cascades,],upper=FALSE,diag=FALSE)
dist_mat_euc <-dist(cbind(dfmeta$DX, dfmeta$DY))

setwd(dir_analysis)
# pdf('F_madogram_obs_vs_sim.pdf',width=10,height=4)
# par(mfrow=c(1,2),mar=c(5,5,2,0))
pdf('F_madogram_obs_vs_sim.pdf',width=5,height=4)
par(mfrow=c(1,1),mar=c(5,5,2,1))
for(r in 1:length(vect_sim)){
  ### euclidean distance
  if(r==1){
    # plot(dist_mat_river,mado_dist,col=col_Q,
    #ylab="F-madogram floods",xlab="River distance [lon/lat]", ylim=c(0.05,0.15))
    ### fit a smoothing spline 
    spl_floods <- smooth.spline(dist_mat_euc,mado_dist[[r]],spar=0.95)
    plot(spl_floods,ylab="F-madogram floods", ylim=c(0.05,0.15),
         xlab="Euclidean Distance [X/Y]",type='l',col='black')
  }else{
    spl_floods <- smooth.spline(dist_mat_euc,mado_dist[[r]],spar=0.95)
    lines(spl_floods,col='orange')
  }
}
#   for(r in 1:length(vect_sim)){
#     ### river distance
#     if(r==1){
#       spl_floods <- smooth.spline(dist_mat_river[which(!is.na(dist_mat_river))],mado_dist[[r]][which(!is.na(dist_mat_river))],spar=0.95)
#       plot(spl_floods,col=col_obs,lwd=2,type='l',ylab="F-madogram floods",
#            xlab="River distance [lon/lat]",ylim=c(0,0.16))
#     }else{
#       spl_floods <- smooth.spline(dist_mat_river[which(!is.na(dist_mat_river))],mado_dist[[r]][which(!is.na(dist_mat_river))],spar=0.95)
#       lines(spl_floods,col=col_sim,lwd=2)
#     }
# }

dev.off()

### (g) extremal coefficient/ tail dependence chibar
### see spatialExtremes by Gilleland
###===============================###===============================###
### tail dependence: estimates tail dependence parameters chi and chibar
### chibar as introduced by Coles 1999
### estimators given by Reiss and Thomas 2007.

### load data of the individual stations
set_name<-'outaouais'
data_stations <- list()
for(l in station_numbers){
  setwd(dir_data)
  # load(file=paste(l,'_',set_name,"_stoch_sim_phase_rand_",'_wave_',n_wave,".RData",sep=""))
  # data_stations[[l]] <- data_stoch
  bassin<- BassinsOutaouais[l]
  load(file=paste(dir_data, bassin,"_stoch_sim_outaouais_Kappa_9997_LD",".Rdata",sep=""))
  data$simulation$timestamp<-NULL
  data_stations[[l]] <- data$simulation
  ### compute specific discharge to compare discharge at different sites
  data_stations[[l]][,4:length(data_stations[[l]])]<- data_stations[[l]][,4:length(data_stations[[l]])] *(1/(aire[l])) * 1000
  ### append a timestamp
  data_stations[[l]]$Date <-   paste(data_stations[[l]]$YYYY,"-",data_stations[[l]]$MM,"-",data_stations[[l]]$DD,sep="")
  data_stations[[l]]$Date <- as.POSIXct(strptime(data_stations[[l]]$Date,format="%Y-%m-%d",tz="GMT"))
  # ### remove first three rows
  data_stations[[l]] <- data_stations[[l]][,-c(1:3)]
}
### remove empty list entries
data_stations <- data_stations[lengths(data_stations) != 0]

### define threshold
thresh <-0.95
### run for different thresholds
for(t in c(0.8,0.95)){
  thresh <- t
  ### run through all simulations runs plus observations
  for(r in 1:length(vect_sim)){
  #for(r in 1:2){
    ### run tail dependence test for pairs of stations.
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
    setwd(dir_analysis)
    save(file=paste(vect_sim[r],'_thresh_',thresh,'_tail_dep.RData',sep=''),chi_mat,chibar_mat,tail_test_mat)
  }
}

### Note: H0 of tail dependence test is: random variables X and Y are dependent. 
### H0 rejected if p-value <0.05. I.e. if p-value smaller than 0.05 independence.
### determine colors
col_dep <- adjustcolor('blue',0.7)
col_indep <- adjustcolor('red',0.3)
col_obs_dep <- adjustcolor('black',0.5)
col_sim_dep <- adjustcolor('green',0.8)
### run through two thresholds
setwd(dir_analysis)
# pdf('tail_dependence_obs_vs_sim.pdf',width=10,height=5)
pdf('tail_dependence_obs_vs_sim_chi.pdf',width=10,height=5)

par(mfrow=c(1,2))
for(t in c(0.8,0.95)){
  thresh<-t
  ### look at results for observations
  load(file=paste(vect_sim[1],'_thresh_',thresh,'_tail_dep.RData',sep='')) ### loads chi_mat,chibar_mat,tail_test_mat)
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



