#Fait par: Fabian Tito Arandia Martinez
#Date: 2020-02-18
#Objectif: Prendre les intrants hydrologiques et les mettre dans une liste.

library(PRSim)
library(tidyverse)
library(lubridate)
library(sparklyr)
library(foreach)
library(doParallel)
rm(list=ls())
setwd("/home/tito/Desktop/PRsim entry")

tests<-list()
fichiers<-list.files()
for(fichier in fichiers){
  
  nom<-strsplit(fichier,'_')[[1]][4]
  date_debut<-strsplit(fichier,'_')[[1]][3]
  date_fin<-substr(strsplit(fichier,'_')[[1]][5],1,8)
  
  
  Qobs<-read.delim(fichier,header = FALSE,col.names = 'Qobs')
  Qobs<-Qobs$Qobs
  
  dates_complete<-seq(as.Date(date_debut,'%Y%m%d'),as.Date("20161231",'%Y%m%d'),by='day')
  dates<-seq(as.Date(date_debut,'%Y%m%d'),as.Date(date_fin,'%Y%m%d'),by='day')
  
  Qobs<-Qobs[dates==dates_complete]
  
  YYYY<-year(dates[dates==dates_complete])
  MM<-month(dates[dates==dates_complete])
  DD<-day(dates[dates==dates_complete])
  
  df<-data.frame(YYYY,MM,DD,Qobs)
  
  
  
  tests[[nom]]<-df
  
}

bvs<-names(tests)
filename2<-paste("/home/tito/Desktop/obs_outaouais.Rdata")
save(tests, file = filename2)

# cores<-detectCores()
# cl<-makeCluster(cores[1]-8)
# registerDoParallel(cl)
# 
# foreach(i = 1:10) %dopar% {
#   require(PRSim)
#   stoch_sim<-prsim.wave(tests, "Qobs", 10, suppWarn=TRUE)
#   names(stoch_sim)<-bvs
#   
#   filename<-paste("/home/tito/Desktop/stoch_sim_100_outaouais_",as.character(i),".Rdata")
#   
#   save(stoch_sim, file = filename)
# 
# }
# 
# stopCluster(cl)

trace(prsim.wave,edit=TRUE)#enlever les set seed.

#for(i in 1:9999) {
  #require(PRSim)
  stoch_sim_70<-prsim.wave(tests, "Qobs", 10, suppWarn=FALSE,win_h_length = 35,GoFtest = 'AD')
  stoch_sim_30<-prsim.wave(tests, "Qobs", 10, suppWarn=FALSE,win_h_length = 15,GoFtest = c('AD','KS'))
  #les sorties des tests de fit sont souvent nulles.
  
  names(stoch_sim)<-bvs
  
  #filename<-paste("/home/tito/Desktop/stoch_sim_10_outaouais_kappa_35dayshalfwindow_",as.character(i),".Rdata",sep='')
  
  #save(stoch_sim, file = filename)
  
#}

##################################
#Analyse des s?ries stochastiques#            
##################################

# load(filename)
 sim <- stoch_sim_30[[1]]$simulation
 par(mai=c(.9,.9,.1,.1))
# 
 plot(sim$timestamp, sim$Qobs, type="l", 
     xlab="Time [d]", ylab=expression(paste("Discharge [m"^3,"/s]")))
# 
 matlines(sim$timestamp, sim[, grep("r", names(sim))],
          lty=1, col="gray")
# 
# #deseasonalized
#  require(deseasonalize)
#  out<-ds(sim$Qobs)
#  des<-out$z
# # # periodogram of deseasonalized
#  kern <- kernel("modified.daniell",c(10,10))#verifier si la taille fait du sens
#  sp1 <- spec.pgram(sim$Qobs, k=kern, taper=0, log="no", plot=FALSE)
#  sp2 <- spec.pgram(des, k=kern, taper=0, log="no", plot=FALSE)
#  plot(sp1, xlim=c(0,.05))     
#  plot( sp2, add=T, col=2)
# # 
# # # Peaks correspond to the following cycles:
#  1/sp1$freq[head(order(sp1$spec, decreasing=TRUE))]
# # 
# # # compare periodogram of simulated series
#  plot(sp1, xlim=c(0,.05))     # would be nice to identify the peaks...
#  for (i in grep("r",names(sim))) {
#    spi <- spec.pgram(sim[,i], k=kern, taper=0, log="no", plot=FALSE)#probleme a elucider
#    plot( spi, add=T, col="gray")
#  }
# # 
#  sp3 <- spec.pgram(sim$Qobs, taper=0, log="no", plot=FALSE)
#  1/sp3$freq[head(order(sp3$spec, decreasing=TRUE))]
# # 
# # 
# # ### plot mean regime for each simulation run and compare to observed regime
# # ### define plotting colors
#  col_sim <- adjustcolor("#fd8d3c",alpha=0.8)
#  col_sim_tran <- adjustcolor("#fd8d3c",alpha=0.2)
#  col_obs <- adjustcolor( "black", alpha.f = 0.2)
# # 
#  year <- unique(sim$YYYY)
# 
# ### compute mean runoff hydrograph
 sim$day_id <- rep(seq(1:365),times=length(year))
 mean_hydrograph_obs <- aggregate(sim$Qobs, by=list(sim$day_id), FUN=mean,simplify=FALSE)
 plot(unlist(mean_hydrograph_obs[,2]), lty=1, lwd=1, col="black", ylab=expression(paste("Discharge [m"^3,"/s]")),
      xlab="Time [d]", main="Mean hydrographs", ylim=c(0,max(unlist(mean_hydrograph_obs[,2]))*1.5),type="l")
# 
# ### add mean runoff hydrographs
 for(r in 6:(length(names(sim))-1)){
   mean_hydrograph <- aggregate(sim[,r], by=list(sim$day_id), FUN=mean,simplify=FALSE)
   lines(mean_hydrograph, lty=1, lwd=1, col=col_sim)
 }
# ### redo observed mean
 lines(mean_hydrograph_obs, lty=1, lwd=1, col="black")
# 
# ### autocorrelation
# acf_mare <- list()
# acf_obs <- acf(sim$Qobs, plot=FALSE)
# plot(acf_obs$acf, type="l", xlab="Lag", main="Autocorrelation", ylab="ACF")
# for(r in 6:(length(names(sim))-2)){
#   acf_sim <- acf(sim[,r], plot=FALSE)#probleme de distribution de R?
#   lines(acf_sim$acf, col=col_sim, type="l")
#   ### compute mean relative error in the acf
#   acf_mare[[r]]<- mean(abs((acf_obs$acf-acf_sim$acf)/acf_obs$acf))
# }
# lines(acf_obs$acf)
# 
# ### partial autocorrelation function
# pacf_obs <- pacf(sim$Qobs, plot=FALSE)
# pacf_mare <- list()
# plot(pacf_obs$acf, type="l", xlab="Lag", main="Partial autocorrelation", ylab="PACF")
# for(r in 6:(length(names(sim))-2)){
#   pacf_sim <- pacf(sim[,r], plot=FALSE)
#   lines(pacf_sim$acf, col=col_sim, type="l")
#   ### compute mean relative error in the acf
#   pacf_mare[[r]] <- mean(abs((pacf_obs$acf-pacf_sim$acf)/pacf_obs$acf))
# }
# lines(pacf_obs$acf)
# 
# ### compute seasonal statistics
# ### Q50,Q05,Q95, boxplots
# ### define seasons: Winter:12,1,2; spring:3,4,5; summer: 6,7,8; fall: 9,10,11
sim$season <- "winter"
sim$season[which(sim$MM%in%c(3,4,5))] <- "spring"
sim$season[which(sim$MM%in%c(6,7,8))] <- "summer"
sim$season[which(sim$MM%in%c(9,10,11))] <- "fall"

### all simulated series show the same seasonal statistics. plot only one
#attention, le graphique n'est valide que pour une simulation.
boxplot(sim$Qobs[which(sim$season=="winter")], sim$r1[which(sim$season=="winter")],
        sim$Qobs[which(sim$season=="spring")], sim$r1[which(sim$season=="spring")],
        sim$Qobs[which(sim$season=="summer")], sim$r1[which(sim$season=="summer")],
        sim$Qobs[which(sim$season=="fall")], sim$r1[which(sim$season=="fall")],
        border=c("black", col_sim, "black", col_sim, "black", col_sim, "black", col_sim),
        xaxt="n", main="Seasonal statistics", outline=FALSE)
mtext(side=1, text=c("Winter", "Spring", "Summer", "Fall"), at=c(1.5,3.5,5.5,7.5))

# ###############################
# #Calcul des pointes mensuelles#            
# ###############################
# 
# df_max_prt<-df %>% na.pass() %>% select(DATE,DEBIT) %>%
#   gather(variable, value, -DATE) %>% mutate(MONTH = as.numeric(format(DATE, "%m")), YEAR = format(DATE, "%Y")) %>%
#   filter(YEAR>1910 & YEAR<2020) %>% 
#   filter(MONTH == 3 | MONTH == 4 | MONTH == 5| MONTH == 6) %>% group_by(variable,YEAR) %>%  summarise(max_year_prt = max(value) )
# 
#




