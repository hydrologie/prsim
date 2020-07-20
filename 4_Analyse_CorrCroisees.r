#Fait par: Fabian Tito Arandia Martinez
#Date de cr?ation: 2020-02-21
#Objectif: Analyse des corr?lations.


##################################
#Analyse des series stochastiques#            
##################################
rm(list=ls())

filename<-paste('/home/tito/Desktop/sims_kappaLD/stoch_sim_10_outaouais_Kappa_1000_LD.Rdata')

load(filename)
#sim <- stoch_sim[[1]]$simulation #Qobs, r1 

simulations_multi_sites<-list(stoch_sim[[1]]$simulation,stoch_sim[[2]]$simulation,stoch_sim[[3]]$simulation,stoch_sim[[4]]$simulation)

sim <- stoch_sim[[1]]$simulation
#sim<-simulations_multi_sites
#deseasonalized
require(deseasonalize)
out<-ds(sim$Qobs)
des<-out$z
# periodogram of deseasonalized
kern <- kernel("modified.daniell",c(10,10))#verifier si la taille fait du sens
sp1 <- spec.pgram(sim$Qobs, k=kern, taper=0, log="no", plot=FALSE)
sp2 <- spec.pgram(des, k=kern, taper=0, log="no", plot=FALSE)
plot(sp1, xlim=c(0,.05))     
plot( sp2, add=T, col=2)

# Peaks correspond to the following cycles:
1/sp1$freq[head(order(sp1$spec, decreasing=TRUE))]

# compare periodogram of simulated series
plot(sp1, xlim=c(0,.05))     # would be nice to identify the peaks...
for (i in grep("r",names(sim))) {
  spi <- spec.pgram(sim[,i], k=kern, taper=0, log="no", plot=FALSE)#probleme a elucider
  plot( spi, add=T, col="gray")
}

sp3 <- spec.pgram(sim$Qobs, taper=0, log="no", plot=FALSE)
1/sp3$freq[head(order(sp3$spec, decreasing=TRUE))]


### plot mean regime for each simulation run and compare to observed regime
### define plotting colors
col_sim <- adjustcolor("#fd8d3c",alpha=0.8)
col_sim_tran <- adjustcolor("#fd8d3c",alpha=0.2)
col_obs <- adjustcolor( "black", alpha.f = 0.2)

year <- unique(sim$YYYY)

### autocorrelation
acf_mare <- list()
acf_obs <- acf(sim$Qobs, plot=FALSE)
plot(acf_obs$acf, type="l", xlab="Lag", main="Autocorrelation", ylab="ACF")
for(r in 6:(length(names(sim))-2)){
  acf_sim <- acf(sim[,r], plot=FALSE)#probleme de distribution de R?
  lines(acf_sim$acf, col=col_sim, type="l")
  ### compute mean relative error in the acf
  acf_mare[[r]]<- mean(abs((acf_obs$acf-acf_sim$acf)/acf_obs$acf))
}
lines(acf_obs$acf)

### partial autocorrelation function
pacf_obs <- pacf(sim$Qobs, plot=FALSE)
pacf_mare <- list()
plot(pacf_obs$acf, type="l", xlab="Lag", main="Partial autocorrelation", ylab="PACF")
for(r in 6:(length(names(sim))-2)){
  pacf_sim <- pacf(sim[,r], plot=FALSE)
  lines(pacf_sim$acf, col=col_sim, type="l")
  ### compute mean relative error in the acf
  pacf_mare[[r]] <- mean(abs((pacf_obs$acf-pacf_sim$acf)/pacf_obs$acf))
}
lines(pacf_obs$acf)









