###===============================###===============================###
### Function PRSim.wave for varying marginals
### Manuela Brunner, NCAR
### 11/27/2020
###===============================###===============================###
rm(list=ls())
library(devtools)
library(PRSim)
# load_all()

# load example data
load("H:/Projets_communs/2020/Outaouais PRsim/02_Calculs/Travail de MB/runoff_multi_sites.rda")
data<- runoff_multi_sites

### define directories
#dir_main <- "~/Projects/stochastic_simulation/hydro_Quebec/" ### path needs to be pointing to the hydro_Quebec folder I sent you.
dir_main <- "H:/Projets_communs/2020/Outaouais PRsim/02_Calculs/Travail de MB/"
# dir_data <- paste(dir_main,"example_data",sep='')
dir_analysis <- paste(dir_main,"results",sep='')


### define marginal distributions
### GEV
require("evd")
require("ismev")
rGEV <- function(n, theta, max_prob)  rgev(n, theta[1], theta[2], theta[3])
pGEV <- function(x, theta)  pgev(x, theta[1], theta[2], theta[3])
GEV_fit <- function( xdat, ...)   gev.fit(xdat, show=FALSE, ...)$mle

### Normal
library(fitdistrplus)
rNORM <- function(n, theta, max_prob)  rnorm(n, theta[1], theta[2])
pNORM <- function(x, theta)  pnorm(x, theta[1], theta[2])
NORM_fit <- function( xdat, ...)   fitdistr( xdat, 'normal', show=FALSE, ...)$estimate

### Gumbel with varying max prob
rgumbelLD <- function(n, theta, max_prob){
  F <- runif(n, min = 1e-10, max = max_prob)
  x = qgumbel(F, theta[1],theta[2])
  return(x)
}
pgumbelLD <- function(x, theta)  pgumbel(x, theta[1], theta[2])
gumbelLD_fit <- function( xdat, ...)   gum.fit( xdat, show=FALSE, ...)$mle

### kappa in rCDF format
rKAPPA <- function (n, theta, max_prob) 
{
  F <- runif(n, min = 1e-10, max = max_prob)
  x <- invF.kappa(F, theta[1], theta[2], theta[3], theta[4])
  return(x)
}
KAPPA_fit <- function(xdat){
  ll <- homtest::Lmoments(xdat)
  test <- try(par.kappa(ll[1],ll[2],ll[4],ll[5]), silent = TRUE)
  if(length(test)>1){
    x <- test
    x <- unlist(x)
  }else{
    x <- NA
  }
  return(x)
}
pKAPPA <- function(x, theta){
  x <- F.kappa(x, theta[1], theta[2], theta[3], theta[4])
  return(x)
}  

### new arguments compared to original prsim.wave function:
### marginal_list instead of marginal: list with marginal distribution names,
### one name for each station in the same order as data is provided
### n_par_list instead of n_par: list with number of parameters per marginal distribution,
### one number for each station in the same order data is provided
marginal_list=list('kappa','kappa','kappa','kappa')
marginal_list <- list('gumbelLD','gumbelLD','gumbelLD','gumbelLD')
marginal_list <- list('KAPPA','KAPPA','KAPPA','KAPPA')

### Implementation of a time varying max probability parameter in rand.kappa
### This will allow to cap distribution if necessary (e.g. in spring)
### 03.02.2021: implement new parameter prob_list
### allow for choosing seasonally variable probabilities to be used as max parameter in rand.kappa function
prob_list<-rep(list(0.997),times=365) ### must contain one probability value per day
### requires redefining rand.kappa function to include an additional probabity parameter

### 24.02.2021
### implement seasonally varying marginals, use daily varying distribution choice
### daily distribution definition station 1, use one distribution per month
### in this example, I am using kappa for Nov-April (winter) and normal for May-Oct (summer)
### distributions station 1
dist_1 <- c(rep(list('kappa'),times=31),rep(list('kappa'),times=28),rep(list('kappa'),times=31), ### jan to march
                         rep(list('kappa'),times=30),rep(list('NORM'),times=31),rep(list('NORM'),times=30), ### april to june
                         rep(list('NORM'),times=31),rep(list('NORM'),times=31),rep(list('NORM'),times=30), ### july to september
                         rep(list('NORM'),times=31),rep(list('kappa'),times=30),rep(list('kappa'),times=31)) ### october to september
### distributions station 2
dist_2 <- c(rep(list('GEV'),times=31),rep(list('GEV'),times=28),rep(list('GEV'),times=31), ### jan to march
            rep(list('GEV'),times=30),rep(list('NORM'),times=31),rep(list('NORM'),times=30), ### april to june
            rep(list('NORM'),times=31),rep(list('NORM'),times=31),rep(list('NORM'),times=30), ### july to september
            rep(list('NORM'),times=31),rep(list('GEV'),times=30),rep(list('GEV'),times=31)) ### october to september
### station 3
dist_3 <- dist_2
### station 4
dist_4 <- dist_1

### store daily distribution choice for each catchment
marginal_list <- list(dist_1,dist_2,dist_3,dist_4)

### the number of parameter list also needs to be updated to a daily format
par_1 <- c(rep(list(4),times=31),rep(list(4),times=28),rep(list(4),times=31), ### jan to march
            rep(list(4),times=30),rep(list(2),times=31),rep(list(2),times=30), ### april to june
            rep(list(2),times=31),rep(list(2),times=31),rep(list(2),times=30), ### july to september
            rep(list(2),times=31),rep(list(4),times=30),rep(list(4),times=31)) ### october to september
par_2 <- c(rep(list(3),times=31),rep(list(3),times=28),rep(list(3),times=31), ### jan to march
            rep(list(3),times=30),rep(list(2),times=31),rep(list(2),times=30), ### april to june
            rep(list(2),times=31),rep(list(2),times=31),rep(list(2),times=30), ### july to september
            rep(list(2),times=31),rep(list(3),times=30),rep(list(3),times=31)) ### october to september
### station 3
par_3 <- par_2
### station 4
par_4 <- par_1
### number of paramters all stations and days
n_par_list <- list(par_1,par_2,par_3,par_4)

### define PRSim.wave function with varying margins depending on station
prsim.wave.marg <- function(data, station_id="Qobs", number_sim=1, win_h_length=15, 
                       marginal_list=list("kappa","empirical","GEV","NORM"), n_par_list=list(4,NA,3,2), n_wave=100, marginalpar=TRUE, 
                       GoFtest=NULL, verbose=TRUE, suppWarn=FALSE, prob_list=rep(list(0.997),times=365),...){  
  
  ### function for backtransformation of continuous wavelet transform
  ### inverse wavelet transform
  ### x is the input matrix
  fun_icwt<-function(x){
    wt.r<-Re(x)
    
    ### define number of scales
    J<-length(x[1,])
    # Reconstruct as in formula (11):
    dial<-2*2^(0:J*.125)
    rec<-rep(NA,(length(x[,1])))
    for(l in 1:(length(x[,1]))){
      rec[l]<-0.2144548*sum(wt.r[l,]/sqrt(dial)[1:length(wt.r[l,])])
    }
    return(rec)
  }
  
  ### redefine kappa distribution to include a new max_prob parameter
  ### will allow us to use a seasonally varying maximum probability
  rand.kappa <- function (numerosita, xi, alfa, k, h, max_prob) 
  {
    F <- runif(numerosita, min = 1e-10, max = max_prob)
    x <- invF.kappa(F, xi, alfa, k, h)
    return(x)
  }
  
  ## start preparing arguments.
  if (!is.null(GoFtest)) {
    GoFtest <- toupper(GoFtest)[1]
    if (!(GoFtest %in% c("AD","KS"))) stop("'GoFtest' should be either 'NULL', 'AD' or 'KS'.")
  } else  GoFtest <- "NULL"
  
  
  ### list for storing distributions
  rCDF <- rep(list(rep(list(NA),times=365)),length(marginal_list))
  CDF_fit <- rep(list(rep(list(NA),times=365)),length(marginal_list))
  ### store distribution for each catchment
  for(l in 1:length(marginal_list)){
    for(d in 1:365){
      marginal <- marginal_list[[l]][[d]]    # check all list elements
      if (!(marginal %in% c("kappa","empirical"))) {   # check if distributions exist
        if (!is.character(marginal)) stop("'marginal' should be a character string.")
        rCDF[[l]][[d]] <- get(paste0("r",marginal), mode = "function")
        CDF_fit[[l]][[d]] <- get(paste0(marginal,"_fit"), mode = "function")
        if (GoFtest=="AD")	  pCDF <- get(paste0("p",marginal), mode = "function")
      }
    }
  }
  
  
  op <- options("warn")$warn
  ### input data needs to be of the format year (four digits), month (two digits), day (one digit), input discharge time series
  
  ### check for correct input data labels and length
  ### run through all stations: list
  for(l in 1:length(data)){
    if (nrow(data[[l]])[1]<730) stop("At least one year of data required.")
    if (is.numeric(station_id)){
      station_id <- colnames(data[[l]])[station_id]
    }
    if (is.na(station_id)||!("Qobs" %in% colnames(data[[l]]))) stop("Wrong column (name) for observations selected.")
    
    # test for proper format:
    if (any(class(data[[l]][,1])%in%c("POSIXct","POSIXt"))){
      data <- data.frame(YYYY=as.integer(format(data[[l]][,1],'%Y')), 
                         MM=as.integer(format(data[[l]][,1],'%m')),
                         DD=as.integer(format(data[[l]][,1],'%d')),
                         Qobs=data[[l]][,station_id],
                         timestamp=data[[l]][,1])
    } else {
      if(!all(c("YYYY","MM","DD") %in% colnames(data[[l]]))) stop("Wrong time column names")
      
      data[[l]] <- data[[l]][,c("YYYY","MM","DD", station_id)]
      tmp <- paste(data[[l]]$YYYY,data[[l]]$MM,data[[l]]$DD,sep=" ")
      names(data[[l]]) <- c("YYYY","MM","DD","Qobs")
      data[[l]]$timestamp <- as.POSIXct(strptime(tmp, format="%Y %m %d", tz="GMT"))
    }
    
    ### remove February 29
    data[[l]] <- data[[l]][format(data[[l]]$timestamp, "%m %d") != "02 29",]
    ### remove incomplete years
    if(which(format(data[[l]]$timestamp,format='%j')=='001')[1]>1){
      data[[l]] <- data[[l]][-c(1:(which(format(data[[l]]$timestamp,format='%j')=='001')[1]-1)),]
    }
    if ((nrow(data[[l]]) %% 365)>0) stop("No missing values allowed. Some days are missing.")
    
    ### replace missing data by mean values
    if(length(which(is.na(data[[l]]$timestamp)))>0){
      ### replace days with missing data
      data[[l]][which(is.na(data[[l]]$timestamp)),]$Qobs <- mean(data[[l]]$Qobs,na.rm=T)
    }
    
    ### replace NA values
    if(length(which(is.na(data[[l]]$Qobs)))>0){
      data[[l]][which(is.na(data[[l]]$Qobs)),]$Qobs <- mean(data[[l]]$Qobs,na.rm=T)
    }
    
    ### generate a day index
    data[[l]]$index <- as.numeric(format(data[[l]]$timestamp,format='%j'))
    ### replace empty index positions
    if(length(which(is.na(data[[l]]$index))>0)){
      data[[l]]$index[which(is.na(data[[l]]$index))] <- rep(c(1:365), times=length(unique(data[[l]]$YYYY)))[which(is.na(data[[l]]$index))]
    }
  }
  
  if (verbose) cat(paste0("Detrending with (half-)length ",win_h_length,"...\n"))  
  
  ### (1) Generation of white noise for random phases generation
  ### generate random sample of indices for each simulation run
  set.seed(10)
  noise_mat_r <- list()
  for (r in 1:number_sim){
    ts_wn <- rnorm(n=length(data[[1]]$Qobs), mean = 0, sd = 1) ### iid time seris
    ### determine scale range
    scale.range = deltat(data[[l]]$Qobs) * c(1, length(data[[l]]$Qobs))
    ### sampling interval
    sampling.interval <- 1
    ### determine octave
    octave <- logb(scale.range, 2)
    ### determine wavelet scales
    scale <- ifelse1(n_wave > 1, 2^c(octave[1] + seq(0, n_wave -
                                                       2) * diff(octave)/(floor(n_wave) - 1), octave[2]), scale.range[1])
    scale <- unique(round(scale/sampling.interval) * sampling.interval)
    
    wt_morlet <- cwt_wst(signal=ts_wn,scales=scale,wname='MORLET',makefigure=FALSE,dt=1,powerscales=FALSE)
    noise_mat_r[[r]] <- as.matrix(wt_morlet$coefs)
  }
  
  ### fitting of kappa distribution to all stations for which simulations are to be derived  
  par_day_list <- rep(list(rep(list(NA),times=365)),length(marginal_list))
  for(l in 1:length(data)){
    ### daily fitting of Kappa distribution
    ### fit the parameters of the Kappa distribution for each day separately.
    ### To enable a sufficient sample size by using daily values in moving window around day i (i.e., reduce uncertainty due to fitting)
    p_vals <- numeric(365) 
    
    # density_kap <- list()
    ### define window length  
    win_length <- c(1:win_h_length)
    for(d in c(1:365)){
    ### data[[l]]$index is somehow overwritten
    if(marginal_list[[l]][[d]]=='empirical'){
      marginal_list[[l]][[d]]<-'empirical'
    }
    if(marginal_list[[l]][[d]]=="kappa"){
      marginal_list[[l]][[d]] <- 'kappa'
      par_day <- matrix(0, nrow=365, ncol=4)
        ### define start and end of window
        before <- data[[l]]$index[d+365-win_length]
        after <- data[[l]]$index[d+365+win_length-1]
        ### define days within window
        ids <- c(before, after)
        ### determine values in window around day i
        data_window <- data[[l]]$Qobs[which(data[[l]]$index%in%ids)]
        # par.kappa(data_monthly)
        ll<- homtest::Lmoments(data_window)
        
        ### test whether Kappa distribution can be fit
        if (suppWarn) {
          suppressWarnings( test <- try(par.kappa(ll[1],ll[2],ll[4],ll[5]), silent = TRUE) )
        } else {
          test <- try(par.kappa(ll[1],ll[2],ll[4],ll[5]), silent = TRUE)
        }
        
        if(length(test)>1){
          kap_par <- test
          par_day[d,] <- unlist(kap_par)
          ### define vector of quantiles
          quant <- sort(data_window)
          thresh <- kap_par$xi + kap_par$alfa*(1 - kap_par$h^(-kap_par$k))/kap_par$k
          if(!is.na(thresh)){
            ##        min(quant)>thresh
            ### only use quantiles larger than threshold (as in f.kappa function documentation)
            quant <- quant[which(quant>thresh)]
          }

          data_kap <- rand.kappa(length(data_window), xi=kap_par$xi,alfa=kap_par$alfa, k=kap_par$k, h=kap_par$h, max_prob=prob_list[[d]])
          ### test effect of varying maximum probabilities
          # data_kap_999 <- rand.kappa(length(data_window), xi=kap_par$xi,alfa=kap_par$alfa, k=kap_par$k, h=kap_par$h, max_prob=0.999)
          # data_kap_997 <- rand.kappa(length(data_window), xi=kap_par$xi,alfa=kap_par$alfa, k=kap_par$k, h=kap_par$h, max_prob=0.997)
          # data_kap_99 <- rand.kappa(length(data_window), xi=kap_par$xi,alfa=kap_par$alfa, k=kap_par$k, h=kap_par$h, max_prob=0.99)
          # data_kap_9 <- rand.kappa(length(data_window), xi=kap_par$xi,alfa=kap_par$alfa, k=kap_par$k, h=kap_par$h, max_prob=0.9)
          # 
          # par(mfrow=c(1,4))
          # plot(sort(data_window),sort(data_kap_999),xlab='obs',ylab='sim',main='max prob = 0.999')
          # abline(0,1)
          # plot(sort(data_window),sort(data_kap_997),xlab='obs',ylab='sim',main='max prob = 0.997')
          # abline(0,1)
          # plot(sort(data_window),sort(data_kap_99),xlab='obs',ylab='sim',main='max prob = 0.99')
          # abline(0,1)
          # plot(sort(data_window),sort(data_kap_9),xlab='obs',ylab='sim',main='max prob = 0.9')
          # abline(0,1)
          ### adjusting max_prob seems to be an efficient way of controlling the heaviness of the tail efficiently
          
          if (tolower(GoFtest)=="ks")
            p_vals[d] <- ks_test(data_window, data_kap) ### kappa distribution not rejected at alpha=0.05
          #        p_vals[d] <- ks.test(data_window, data_kap)$p.value ### kappa distribution not rejected at alpha=0.05
          if (tolower(GoFtest)=="ad") {
            
            try_ad_test <- try(ad.test(data_window,F.kappa,xi=kap_par$xi,alfa=kap_par$alfa,k=kap_par$k,h=kap_par$h), silent=TRUE) 
            if(length(try_ad_test)==1){
              p_vals[d]  <- NA
            }else{
              p_vals[d]  <- try_ad_test$p.value
            }
          }
        } else{
          if(d==1){
            p_vals[d] <- NA
            par_day[d,] <- NA
          }else{
            p_vals[d] <- p_vals[d-1]
            par_day[d,] <- par_day[d-1,]
          }
        }
      }
      
      
      ### Treatment for the case when Kappa distribution can not be fitted
      ### a) parameters can't be fitted for any of the days
      if(length(which(is.na(par_day[,1])))==365){
        ### use empirical distribution instead
        marginal_list[[l]][[d]]<-'empirical'
      } else{
        ### b) parameters can be fitted for some days
        ### replace NA entries by values estimated for subsequent day
        if(length(which(is.na(par_day[,1])))>0){
          indices <- rev(which(is.na(par_day[,1])))
          for(i in 1:length(indices)){
            par_day[indices[i],] <- par_day[indices[i]+1,]
          }
        }
      }
      par_day_list[[l]][[d]] <- par_day[d,]
    # }      
    
    ### use either a predefined distribution in R or define own function
    if(marginal_list[[l]][[d]]!="kappa" & marginal_list[[l]][[d]]!="empirical"){
      marginal_list[[l]][[d]] <- marginal_list[[l]][[d]]
      p_vals <- numeric(365) 
      par_day <- matrix(0, nrow=365, ncol=n_par_list[[l]][[d]])
      # for(d in c(1:365)){
        ### define window length
        win_length <- seq(1:15)
        ### define start and end of window
        before <- data[[l]]$index[d+365-win_length]
        after <- data[[l]]$index[d+365+win_length-1]
        ### define days within window
        ids <- c(before,after)
        ### determine values in window around day i
        data_window <- data[[l]]$Qobs[which(data[[l]]$index%in%ids)]
        theta <-  CDF_fit[[l]][[d]](xdat=data_window)
        
        # density_gengam[[d]] <- density(data_gengam)
        # hist(data_window)
        # hist(data_random,add=T,col="red")
        if (tolower(GoFtest)=="ks"){
          ### goodness of fit test
          data_random <- rCDF[[l]][[d]](n=length(data_window), theta, max_prob=prob_list[[d]])
          p_vals[d] <- ks_test(data_window,data_random)
          #            p_vals[d] <- ks.test(data_window,data_random)$p.value 
        }
        if (tolower(GoFtest)=="ad"){
          p_vals[d] <-  ad.test(data_window,pCDF,theta)$p.value
        }
        ### store parameters
        par_day[d,] <- theta
      } 
      par_day_list[[l]][[d]] <- par_day[d,]
    }
  }
  
  
  ### replace NA values by mean if necessary: otherwise, problems with transformation
  for(l in 1:length(data)){
    ### center_data: substract mean from values
    data[[l]]$norm <- data[[l]]$Qobs-mean(data[[l]]$Qobs,na.rm=T)
  }
  
  ### repeat stochastic simulation several times
  
  if(verbose) cat(paste0("Starting ",number_sim," simulations:\n"))
  
  
  ### run through all stations
  out_list<-list()
  for(l in 1:length(data)){
    ### list for storing results
    data_sim <- list()
    ### simulate n series
    for (r in c(1:number_sim)){
      ### use alternative R-package instead
      # wt_morlet <- WaveletTransform(x=data[[l]]$norm,dt=1,dj=1/8)
      ### determine scale range
      scale.range = deltat(data[[l]]$norm) * c(1, length(data[[l]]$norm))
      ### sampling interval
      sampling.interval <- 1
      ### determine octave
      octave <- logb(scale.range, 2)
      ### determine wavelet scales
      scale <- ifelse1(n_wave > 1, 2^c(octave[1] + seq(0, n_wave -
                                                         2) * diff(octave)/(floor(n_wave) - 1), octave[2]), scale.range[1])
      scale <- unique(round(scale/sampling.interval) * sampling.interval)
      ### these scales correspond to the scales originall used in wavCWT()
      
      ### apply continuous wavelet transform: use package wavScalogram
      wt_morlet <- cwt_wst(signal=data[[l]]$norm,scales=scale,wname='MORLET',
                           powerscales=FALSE,makefigure=FALSE,dt=1,wparam=5)
      ### return CWT coefficients as a complex matrix with rows and columns representing times and scales, respectively.
      morlet_mat <- as.matrix(wt_morlet$coefs)
      
      ### something is wrong with the scale of modulus
      
      ### derive modulus of complex numbers (radius)
      modulus <- Mod(morlet_mat)
      ### extract phases (argument)
      phases <- Arg(morlet_mat)
      
      ### use the noise matrix corresponding to this run
      noise_mat <- noise_mat_r[[r]]
      phases_random <- Arg(noise_mat)
      
      ### iv) combine this randomised phase and the WT modulus of the original signal to obtain a surrogate time-frequency distribution
      ### create a new matrix
      ### combine modulus of original series to randomised phase: create new matrix of complex values
      mat_new <- matrix(complex(modulus=modulus,argument=phases_random),ncol=ncol(phases_random))
      ### plug into the original time-frequency object
      ### wmtsa package does not allow for the inverser transform of a CWT object
      
      ### v) inverse wavelet transform
      ### apply inversion to CWT of original data
      rec_orig = fun_icwt(x=morlet_mat)+mean(data[[l]]$Qobs)
      ### apply wavelet reconstruction to randomized signal
      rec<- fun_icwt(x=mat_new)        
      ### add mean
      rec_random<-rec+mean(data[[l]]$Qobs)
      
      ### create new data frame
      data_new <- data.frame("random"=rec_random)
      
      ### add months and years
      data_new$MM <- data[[l]]$MM
      data_new$DD <- data[[l]]$DD
      data_new$YYYY <- data[[l]]$YYYY
      data_new$index <- data[[l]]$index
      
      ### use transformed data directly
      data_new$seasonal <- data_new$random
      # ### derive the ranks of the data
      data_new$rank <- rank(data_new$seasonal)
      
      ### vi) rescale the surrogate to the distribution of the original time series 
      ### apply daily backtransformation: ensures smoothness of regimes
      d<-1
      data_new$simulated_seasonal <- NA
      
      
      for(d in c(1:365)){    
        data_day <- data[[l]][which(data[[l]]$index%in%c(d)),]
        
        ### use kappa distribution for backtransformation
        if(marginal_list[[l]][[d]]=="kappa"){
          names(par_day_list[[l]][[d]]) <- names(kap_par)
          
          ### use monthly Kappa distribution for backtransformation
          ### simulate random sample of size n from Kappa disribution
          data_day$kappa <- rand.kappa(length(data_day$Qobs),
                            xi=par_day_list[[l]][[d]][1],alfa=par_day_list[[l]][[d]][2],
                            k=par_day_list[[l]][[d]][3],h=par_day_list[[l]][[d]][4],max_prob=prob_list[[d]])
          
          
          data_day$rank <- rank(data_day$kappa)
          
          data_new$rank <- rank(data_new$seasonal)
          data_new$rank[ which(data[[l]]$index%in%c(d)) ] <- rank(data_new[which(data[[l]]$index%in%c(d)), ]$seasonal)
          ### derive corresponding values from the kappa distribution
          ### identify value corresponding to rank in the kappa time series
          data_ordered <- data_day[order(data_day$rank),]
          data_new$simulated_seasonal[which(data_new$index%in%c(d))] <- data_ordered$kappa[data_new$rank[which(data[[l]]$index%in%c(d))]]
          ### if error was applied, replace negative values by 0 values
          ### in any case, replace negative values by 0. Corresponds to a bounded Kappa distribution
          
          if(length(which(data_new$simulated_seasonal<0))>0){
            ### do not use 0 as a replacement value directly
            # data_new$simulated_seasonal[which(data_new$simulated_seasonal<0)] <- 0
            ### sample value from a uniform distribution limited by 0 and the minimum observed value
            ### determine replacement value
            rep_value <- runif(n=1,min=0,max=min(data_day$Qobs))
            data_new$simulated_seasonal[which(data_new$simulated_seasonal<0)]<-rep_value
          } 
        }
        ### use empirical distribution for backtransformation
        if(marginal_list[[l]]=="empirical"){
          data_day$rank <- rank(data_day$Qobs)
          data_new$rank <- rank(data_new$seasonal)
          data_new$rank[which(data[[l]]$index%in%c(d))] <- rank(data_new[which(data[[l]]$index%in%c(d)),]$seasonal)
          ### derive corresponding values from the empirical distribution
          ### identify value corresponding to rank in the original time series
          data_ordered <- data_day[order(data_day$rank),]
          data_new$simulated_seasonal[which(data_new$index%in%c(d))] <- data_ordered$Qobs[data_new$rank[which(data[[l]]$index%in%c(d))]]
          # }
        }
        
        ### use any predefined distribution for backtransformation
        if(marginal_list[[l]][[d]]!="kappa" & marginal_list[[l]][[d]]!="empirical"){
          ### use monthly distribution for backtransformation
          ### simulate random sample of size n from disribution
          if(is.na(par_day_list[[l]][[d]])){
            ### use mean parameters if parameters could not be fitted (concerns KAPPA)
            data_day$cdf <-   rCDF[[l]][[d]](n=length(data_day$Qobs), 
                                             colMeans(par_day_list[[l]][[d]],na.rm=TRUE),max_prob=prob_list[[d]])
          }else{
            data_day$cdf <-   rCDF[[l]][[d]](n=length(data_day$Qobs), 
                                             par_day_list[[l]][[d]],max_prob=prob_list[[d]])
          }
          data_day$rank <- rank(data_day$cdf)
          
          data_new$rank <- rank(data_new$seasonal)
          
          # hist(data_day$Qobs)
          # hist(data_day$cdf,add=T,col="blue")
          # data_day$rank <- rank(data_day$cdf)
          data_new$rank[which(data[[l]]$index%in%c(d))] <- rank(data_new[which(data[[l]]$index%in%c(d)),]$seasonal)
          ### derive corresponding values from the kappa distribution
          ### identify value corresponding to rank in the kappa time series
          data_ordered <- data_day[order(data_day$rank),]
          data_new$simulated_seasonal[which(data_new$index%in%c(d))] <- data_ordered$cdf[data_new$rank[which(data[[l]]$index%in%c(d))]]
        }
      }  # end for loop
      data_sim[[r]] <- data_new$simulated_seasonal
      
      if(verbose) cat(".")
      ### next simulation run
    }
    if(verbose) cat("\nFinished.\n")
    ### put observed and simulated data into a data frame
    data_sim <- as.data.frame(data_sim)
    names(data_sim) <- paste("r",seq(1:number_sim),sep="")
    data_stoch <- data.frame(data[[l]][,c("YYYY", "MM", "DD", "timestamp", "Qobs")],
                             data_sim)
    
    if (GoFtest=="NULL") {  
      p_vals <- NULL 
    }
    
  
    ### store values in list
    if(marginal != "empirical"){
      if (marginalpar) {  # also return intermediate results
        # return(list(simulation=data_stoch, pars=par_day, p_val=p_vals))
        out_list[[l]] <- list(simulation=data_stoch, pars=par_day, p_val=p_vals)
      } else {
        # return(list(simulation=data_stoch, pars=NULL, p_val=p_vals)) 
        out_list[[l]] <- list(simulation=data_stoch, pars=NULL, p_val=p_vals)
      }
    }else{
      # return(list(simulation=data_stoch))
      out_list[[l]] <- list(simulation=data_stoch, pars=NULL, p_val=NULL)
    }
    # }
    ### to to next station
  }
  # if(is.na(out_dir)){
  return(out_list)
  # }
}

### application example
### simulate using four different distributions
out <- prsim.wave.marg(data=runoff_multi_sites, number_sim=1, 
                       marginal_list=list("kappa","empirical","GEV","NORM"),
                       n_par_list=list(4,NA,3,2), 
                       GoFtest="KS")

### visualize simulations for all test stations
### store stochastically simulated time series
setwd(dir_analysis)
pdf('test_diff_marginals.pdf',width=8,height=8)
par(mfrow = c(4,2),mar=c(2,3,2,1))
### station l
for(l in 1:length(data)){
  sim <- out[[l]]$simulation
  
  ### plot example of simulated time series
  par(mai=c(.9,.9,.1,.1))
  ### observed time series
  plot(sim$timestamp[1:1000], sim$Qobs[1:1000], type="l", main=marginal_list[[l]],
       xlab="Time [d]", ylab=expression(paste("Discharge [m"^3,"/s]")))
  ### add simulations
  matlines(sim$timestamp[1:1000], sim[1:1000, grep("r", names(sim))],
           lty=1, col="gray")
  
  ### compare distributions without outliers
  ### without outliers
  # boxplot(sim$Qobs,sim$r1,outline=F,col=c('black','grey'),names=c('Obs','Sim'))
  ### with outliers
  boxplot(sim$Qobs,sim$r1,col=c('black','grey'),names=c('Obs','Sim'))
}
dev.off()

### apply function using seasonally varying maximum probabilities for kappa distribution
### compile daily maximum probability list:
### lower values for winter/spring to cap distribution (0.997): Dec-may
### and higher values for summer/fall to allow for heavier tail (0.999): Jun-Nov
### list can be adjusted as desired.
prob_list <- c(rep(list(0.997),times=31),rep(list(0.997),times=28),rep(list(0.997),times=31), ### jan to march
  rep(list(0.997),times=30),rep(list(0.997),times=31),rep(list(0.999),times=30), ### april to june
  rep(list(0.999),times=31),rep(list(0.999),times=31),rep(list(0.999),times=30), ### july to september
  rep(list(0.999),times=31),rep(list(0.999),times=30),rep(list(0.997),times=31)) ### october to september

### run for Kappa
out <- prsim.wave.marg(data=runoff_multi_sites, number_sim=1, 
                       marginal_list=list("kappa","kappa","kappa","kappa"),
                       n_par_list=list(4,4,4,4), 
                       GoFtest=NULL,prob_list=prob_list)

### look at seasonal distributions
setwd(dir_analysis)
pdf('winter_vs_summer_dist.pdf',width=8,height=4)
par(mfrow=c(2,4),mar=c(4,4,2,1))
for(l in 1:4){
  sim <- out[[l]]$simulation
  ### determine winter sim
  sim_w <- sim[which(sim$MM%in%c(12,01,02,03,04,05)),]
  ### summer sim
  sim_s <- sim[which(sim$MM%in%c(06,07,08,09,10,11)),]
  ### compare observed to simulated distributions
  plot(sort(sim_w$Qobs),sort(sim_w$r1,na.last=FALSE),xlab='obs',ylab='sim',main=paste('winter',l,sep=' '))
  abline(0,1)
  plot(sort(sim_s$Qobs),sort(sim_s$r1,na.last=FALSE),xlab='obs',ylab='sim',main=paste('summer',l,sep=' '))
  abline(0,1)
}
dev.off()

### run for gumbelLD (i.e. specific rCDF)
out <- prsim.wave.marg(data=runoff_multi_sites, number_sim=1, 
                       marginal_list=list("gumbelLD","gumbelLD","gumbelLD","gumbelLD"),
                       n_par_list=list(4,4,4,4), 
                       GoFtest=NULL,prob_list=prob_list)
### run for external KAPPA
out <- prsim.wave.marg(data=runoff_multi_sites, number_sim=1, 
                       marginal_list=list("KAPPA","KAPPA","KAPPA","KAPPA"),
                       n_par_list=list(4,4,4,4), 
                       GoFtest=NULL,prob_list=prob_list)

### application using seasonal distribution choices
###===============================###===============================###
### define a daily distribution vector for each catchment
### you can either use the same vector for each catchment or a different combination of distributions for each catchment
### distributions station 1
dist_1 <- c(rep(list('kappa'),times=31),rep(list('kappa'),times=28),rep(list('kappa'),times=31), ### jan to march
            rep(list('kappa'),times=30),rep(list('NORM'),times=31),rep(list('NORM'),times=30), ### april to june
            rep(list('NORM'),times=31),rep(list('NORM'),times=31),rep(list('NORM'),times=30), ### july to september
            rep(list('NORM'),times=31),rep(list('kappa'),times=30),rep(list('kappa'),times=31)) ### october to september
### distributions station 2
dist_2 <- c(rep(list('GEV'),times=31),rep(list('GEV'),times=28),rep(list('GEV'),times=31), ### jan to march
            rep(list('GEV'),times=30),rep(list('NORM'),times=31),rep(list('NORM'),times=30), ### april to june
            rep(list('NORM'),times=31),rep(list('NORM'),times=31),rep(list('NORM'),times=30), ### july to september
            rep(list('NORM'),times=31),rep(list('GEV'),times=30),rep(list('GEV'),times=31)) ### october to september
### station 3
dist_3 <- dist_2
### station 4
dist_4 <- dist_1

### combine all distribution vectors in list
### store daily distribution choice for each catchment
marginal_list <- list(dist_1,dist_2,dist_3,dist_4)

### also create a daily vector of the number of distribution parameters
### the number of parameter list also needs to be updated to a daily format
par_1 <- c(rep(list(4),times=31),rep(list(4),times=28),rep(list(4),times=31), ### jan to march
           rep(list(4),times=30),rep(list(2),times=31),rep(list(2),times=30), ### april to june
           rep(list(2),times=31),rep(list(2),times=31),rep(list(2),times=30), ### july to september
           rep(list(2),times=31),rep(list(4),times=30),rep(list(4),times=31)) ### october to september
par_2 <- c(rep(list(3),times=31),rep(list(3),times=28),rep(list(3),times=31), ### jan to march
           rep(list(3),times=30),rep(list(2),times=31),rep(list(2),times=30), ### april to june
           rep(list(2),times=31),rep(list(2),times=31),rep(list(2),times=30), ### july to september
           rep(list(2),times=31),rep(list(3),times=30),rep(list(3),times=31)) ### october to september
### station 3
par_3 <- par_2
### station 4
par_4 <- par_1
### combine all parameters in one list
### number of paramters all stations and days
n_par_list <- list(par_1,par_2,par_3,par_4)

### apply function to four stations using seasonally varying distributions as defined above
out <- prsim.wave.marg(data=runoff_multi_sites, number_sim=1, 
                       marginal_list=marginal_list,
                       n_par_list=n_par_list, 
                       GoFtest=NULL,prob_list=prob_list)

plot(data[[4]]$Qobs,typ="l")
lines(out[[4]]$simulation$r1,col="red")

length(data[[1]]$Qobs)
length(out[[1]]$simulation$r1)
