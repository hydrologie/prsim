#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
###===============================###===============================###
### Function PRSim.wave for varying marginals
### Manuela Brunner, NCAR
### 11/27/2020
###===============================###===============================###

rm(list=ls())
library(devtools)
library(PRSim)
library(PearsonDS)


### define directories
#dir_main <- "/home/tito/Documents/Travail de MB/" ### path needs to be pointing to the hydro_Quebec folder I sent you.
#dir_main<-"H:\\Projets_communs\\2020\\Outaouais PRsim\\01_Intrants\\AnalyseSpatial60Series\\"

dir_main<-"H:\\Projets_communs\\2020\\Outaouais PRSIM\\02_Calculs\\Resultats\\"
dir_analysis <- paste(dir_main,"results",sep='')

### define marginal distributions
### GEV
require("evd")
require("ismev")
rGEV <- function(n, theta)  rgev(n, theta[1], theta[2], theta[3])
pGEV <- function(x, theta)  pgev(x, theta[1], theta[2], theta[3])
GEV_fit <- function( xdat, ...)   gev.fit(xdat, show=FALSE, ...)$mle

### Normal
library(fitdistrplus)
rNORM <- function(n, theta)  rnorm(n, theta[1], theta[2])
pNORM <- function(x, theta)  pnorm(x, theta[1], theta[2])
NORM_fit <- function( xdat, ...)   fitdistr( xdat, 'normal', show=FALSE, ...)$estimate

### Lognormal
rLNORM <- function(n, theta)  rlnorm(n, theta[1], theta[2])
pLNORM <- function(x, theta)  plnorm(x, theta[1], theta[2])
LNORM_fit <- function( xdat, ...)   fitdistr( xdat, "log-normal", show=FALSE, ...)$estimate


### GUMBEL
dgumbel <- function(x, a, b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
pgumbel <- function(q, a, b) exp(-exp((a-q)/b))
qgumbel <- function(p, a, b) a-b*log(-log(p))

rGUMBEL <- function(n, theta)  rgumbel(n, theta[1], theta[2])
pGUMBEL <- function(x, theta)  pgumbel(x, theta[1], theta[2])
GUMBEL_fit <- function( xdat, ...)   fitdist(xdat, "gumbel",start=list(a=10, b=10))$estimate

### GAMMA
rGAMMA <- function(n, theta)  rgamma(n, theta[1], theta[2])
pGAMMA <- function(x, theta)  pgamma(x, theta[1], theta[2])
GAMMA_fit <- function( xdat, ...)   fitdist(xdat, "gamma",lower=c(0,0),start=list(scale=1, shape=1))$estimate

### PIII  Pearson III
dPIII<-function(x, shape, location, scale) PearsonDS::dpearsonIII(x, shape, location, scale, log=FALSE)
pPIII<-function(q, shape, location, scale) PearsonDS::ppearsonIII(q, shape, location, scale, lower.tail = TRUE, log.p = FALSE)
qPIII<-function(p, shape, location, scale) PearsonDS::qpearsonIII(p, shape, location, scale, lower.tail = TRUE, log.p = FALSE)

rPIII <- function(n, theta)  rpearsonIII(n, theta[1], theta[2], theta[3])
PIII_fit <- function( xdat, ...)   fitdist(xdat, "PIII", method = "mse", start=list(shape=1, location=1, scale=1))$estimate




### new arguments compared to original prsim.wave function:
### marginal_list instead of marginal: list with marginal distribution names,
### one name for each station in the same order as data is provided
### n_par_list instead of n_par: list with number of parameters per marginal distribution,
### one number for each station in the same order data is provided
prsim.wave.marg <- function(data, station_id="Qobs", number_sim=1, win_h_length=15, 
                            marginal_list=list("kappa","empirical","GEV","NORM"), n_par_list=list(4,NA,3,2), n_wave=100, marginalpar=TRUE, 
                            GoFtest=NULL, verbose=TRUE, suppWarn=FALSE, ...){  
  
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
  
  ## start preparing arguments.
  if (!is.null(GoFtest)) {
    GoFtest <- toupper(GoFtest)[1]
    if (!(GoFtest %in% c("AD","KS"))) stop("'GoFtest' should be either 'NULL', 'AD' or 'KS'.")
  } else  GoFtest <- "NULL"
  
  
  ### list for storing distributions
  rCDF <- list()
  CDF_fit <- list()
  ### store distribution for each catchment
  for(l in 1:length(marginal_list)){
    marginal <- marginal_list[[l]]    # check all list elements
    if (!(marginal %in% c("kappa","empirical"))) {   # check if distributions exist
      if (!is.character(marginal)) stop("'marginal' should be a character string.")
      rCDF[[l]] <- get(paste0("r",marginal), mode = "function")
      CDF_fit[[l]] <- get(paste0(marginal,"_fit"), mode = "function")
      if (GoFtest=="AD")	  pCDF <- get(paste0("p",marginal), mode = "function")
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
  #set.seed(10)
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
  par_day_list <- list()
  for(l in 1:length(data)){
    ### daily fitting of Kappa distribution
    ### fit the parameters of the Kappa distribution for each day separately.
    ### To enable a sufficient sample size by using daily values in moving window around day i (i.e., reduce uncertainty due to fitting)
    
    ### data[[l]]$index is somehow overwritten
    if(marginal_list[[l]]=='empirical'){
      marginal_list[[l]]<-'empirical'
    }
    if(marginal_list[[l]]=="kappa"){
      marginal_list[[l]] <- 'kappa'
      p_vals <- numeric(365) 
      par_day <- matrix(0, nrow=365, ncol=4)
      # density_kap <- list()
      ### define window length  
      win_length <- c(1:win_h_length)
      for(d in c(1:365)){
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
          
          data_kap <- rand.kappa(length(data_window), xi=kap_par$xi,alfa=kap_par$alfa, k=kap_par$k, h=kap_par$h)
          
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
        marginal_list[[l]]<-'empirical'
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
      par_day_list[[l]] <- par_day
    }      
    
    ### use either a predefined distribution in R or define own function
    if(marginal_list[[l]]!="kappa" & marginal_list[[l]]!="empirical"){
      marginal_list[[l]] <- marginal_list[[l]]
      p_vals <- numeric(365) 
      par_day <- matrix(0, nrow=365, ncol=n_par_list[[l]])
      for(d in c(1:365)){
        ### define window length
        win_length <- seq(1:15)
        ### define start and end of window
        before <- data[[l]]$index[d+365-win_length]
        after <- data[[l]]$index[d+365+win_length-1]
        ### define days within window
        ids <- c(before,after)
        ### determine values in window around day i
        data_window <- data[[l]]$Qobs[which(data[[l]]$index%in%ids)]
        theta <-  CDF_fit[[l]](xdat=data_window)
        
        ### goodness of fit test
        data_random <- rCDF[[l]](n=length(data_window), theta)
        
        # density_gengam[[d]] <- density(data_gengam)
        # hist(data_window)
        # hist(data_random,add=T,col="red")
        if (tolower(GoFtest)=="ks"){
          p_vals[d] <- ks_test(data_window,data_random)
          #            p_vals[d] <- ks.test(data_window,data_random)$p.value 
        }
        if (tolower(GoFtest)=="ad"){
          p_vals[d] <-  ad.test(data_window,pCDF,theta)$p.value
        }
        ### store parameters
        par_day[d,] <- theta
      } 
      par_day_list[[l]] <- par_day
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
        if(marginal_list[[l]]=="kappa"){
          colnames(par_day_list[[l]]) <- names(kap_par)
          
          ### use monthly Kappa distribution for backtransformation
          ### simulate random sample of size n from Kappa disribution
          data_day$kappa <- rand.kappa(length(data_day$Qobs),
                                       xi=par_day_list[[l]][d,"xi"],alfa=par_day_list[[l]][d,"alfa"],
                                       k=par_day_list[[l]][d,"k"],h=par_day_list[[l]][d,"h"])
          
          
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
        if(marginal_list[[l]]!="kappa" & marginal_list[[l]]!="empirical"){
          ### use monthly distribution for backtransformation
          ### simulate random sample of size n from disribution
          data_day$cdf <-   rCDF[[l]](n=length(data_day$Qobs), par_day_list[[l]][d,])
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
#filename<-paste("/media/tito/TIIGE/PRSIM/obs_outaouais_harm_bassins_sup.Rdata")
filename<-paste(dir_main, "obs_outaouais_harm.Rdata", sep='')
load(filename)
runoff_multi_sites<-tests


# start_sim_number<-as.numeric(args[1L])
# print(start_sim_number)
# 
# liste_des_distributions_choisies<-list('LNORM','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa')
# nbre_de_parametres_par_distribution<-list(2,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4)

# [1] "Bark Lake"                    "Baskatong"                    "Cabonga"                     
# [4] "Carillon et Hull"             "Chelsea"                      "Chenaux et Noire"            
# [7] "Chute-des-Chats"              "Des Joachims"                 "Dozois"                      
# [10] "High Falls"                   "Kamaniskeg"                   "Kiamika"                     
# [13] "Kipawa"                       "Lac des Quinze"               "Lac du poisson blanc"        
# [16] "Lac Temiscamingue a Angliers" "Lac Victoria et lac Granet"   "Lady Evelyn"                 
# [19] "Lower Notch et Indian Chute"  "Madawaska-Arnprior"           "Maniwaki"                    
# [22] "Masson"                       "Mistinikon"                   "Mitchinamecus"               
# [25] "Mont-Laurier"                 "Mountain Chute"               "Otto Holden"                 
# [28] "Paugan"                       "Rabbit Lake"                  "Rapide-2"                    
# [31] "Rapide-7"                     "Riviere Blanche"              "Riviere Bonnechere"          
# [34] "Riviere Coulonge"             "Riviere Dumoine"              "Riviere Kinojevis"           
# [37] "Riviere Mattawa"              "Riviere Mississippi"          "Riviere Petawawa"            
# [40] "Riviere Petite Nation"        "Riviere Rideau"               "Riviere Rouge"               
# [43] "Riviere South Nation"        

#Loi Gamma: Lac des Quinzes, Riv Mattawa, ne marchent pas, les bassins madawaska, Rapide-2 et Rapide-7 sont gardes LOGNORMAL 

liste_des_distributions_choisies<-list('kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa','kappa')
nbre_de_parametres_par_distribution<-list(2,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,2)

liste_des_distributions_choisies[1:43]<- 'GUMBEL'
nbre_de_parametres_par_distribution[1:43]<- 2
# #Baskatong
# liste_des_distributions_choisies[2]<- 'GUMBEL'
# nbre_de_parametres_par_distribution[2]<- 2
# 
# #
# liste_des_distributions_choisies[3]<- 'GUMBEL'
# nbre_de_parametres_par_distribution[3]<- 2
# 
# #
# liste_des_distributions_choisies[4]<- 'GUMBEL'
# nbre_de_parametres_par_distribution[4]<- 2
# 
# #Mattawa
# liste_des_distributions_choisies[37]<- 'GUMBEL'
# nbre_de_parametres_par_distribution[37]<- 2
# 
# #Madawaska
# liste_des_distributions_choisies[20]<- 'LNORM'
# nbre_de_parametres_par_distribution[20]<- 2
# #Rapide-2
# liste_des_distributions_choisies[30]<- 'LNORM'
# nbre_de_parametres_par_distribution[30]<- 2
# # Rapide-7
# liste_des_distributions_choisies[31]<- 'LNORM'
# nbre_de_parametres_par_distribution[31]<- 2


# CHANGE NAN VALUES AT THE END OF QOBS
runoff_multi_sites<-tests
data<-runoff_multi_sites
for (l in 1:length(data)){
  ### replace NA values
  if(length(which(is.na(data[[l]]$Qobs)))>0){
    data[[l]][which(is.na(data[[l]]$Qobs)),]$Qobs <- mean(data[[l]]$Qobs,na.rm=T)
    index=length(which(is.na(data[[l]]$Qobs)))>0
    # print(index)
    # print(names(data[l]))
  }
  # remplace les valeurs nulles ou negaatives par la moyenne
  index = which(data[[l]]$Qobs<= 0)
  if (length(index)>0){
    data[[l]][index,]$Qobs <- mean(data[[l]]$Qobs,na.rm=T)
    print(index)
    print(names(data[l]))
  }
}
runoff_multi_sites<- data

# Change for the number of SIMS
start_sim_number=1
for(i in start_sim_number:(start_sim_number)) {
  out <- prsim.wave.marg(data=runoff_multi_sites, number_sim=1, 
                        marginal_list=liste_des_distributions_choisies,
                        n_par_list=nbre_de_parametres_par_distribution)#probleme avec goftest ks_test ks.test

  ### Save the simulations
  names(out)<-names(runoff_multi_sites)
  # dir.create(paste0(dir_main, 'sims_final/'), showWarnings = FALSE) #stops warnings if folder already exists
  # filename<-paste(dir_main, "sims_final/stoch_sim_10_outaouais_Kappa_",as.character(i),"_9997_MB.Rdata",sep='')
  # 
  # save(out, file = filename)
  
  ### visualize simulations for all test stations
  ### store stochastically simulated time series
  setwd(dir_analysis)
  pdf('test_diff_marginals_test.pdf',width=8,height=8)
  marginal_list=liste_des_distributions_choisies
  par(mfrow = c(4,2),mar=c(2,3,2,1))
  ### station l
  for(l in 1:length(data)){
    sim <- out[[l]]$simulation
    title= paste(marginal_list[[l]], "_", names(data[l]))
    ### plot example of simulated time series
    par(mai=c(.9,.9,.1,.1))
    ### observed time series
    sz = length(sim$Qobs)
    plot(sim$timestamp[1:sz], (sim$Qobs[1:sz]), type="l", main=title,
        xlab="Time [d]", ylab=expression(paste("Discharge [m"^3,"/s]")))
    ### add simulations
    matlines(sim$timestamp[1:sz], (sim[1:sz, grep("r", names(sim))]),
            lty=1, col="gray")
    
    ### compare distributions without outliers
    ### without outliers
    # boxplot(sim$Qobs,sim$r1,outline=F,col=c('black','grey'),names=c('Obs','Sim'))
    ### with outliers
    boxplot((sim$Qobs),(sim$r1),col=c('black','grey'),names=c('Obs','Sim'))
  }
  dev.off()
}