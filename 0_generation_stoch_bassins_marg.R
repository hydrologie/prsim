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


source('prsim_wave_marg.R')


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
  
#ajouter visualisation
}