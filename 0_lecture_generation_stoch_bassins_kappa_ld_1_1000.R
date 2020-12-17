#!/usr/bin/env Rscript

library(PRSim)
library(tidyverse)
library(lubridate)
library(sparklyr)
library(foreach)
library(doParallel)
setwd("/media/tito/TIIGE/PRSIM/ApportPRSIM/ApportHarmonise_BassinSup/")
args = commandArgs(trailingOnly = TRUE)
tests<-list()
fichiers<-list.files()
for(fichier in fichiers){
  
  nom<-strsplit(fichier,'_')[[1]][4]
  
  date_debut<-strsplit(fichier,'_')[[1]][3]
  date_debut_sim<-"19700101" #adapte pour les bassins supplementaires
  date_fin<-substr(strsplit(fichier,'_')[[1]][5],1,8)
  
 
  Qobs<-read.delim(fichier,header = TRUE,sep = ',')#trouver une facon de rendre ceci robuste
  col_titles<-colnames(Qobs)
  Qobs<-Qobs[[col_titles[2]]]
  
  
  dates_complete<-seq(as.Date(date_debut_sim,'%Y%m%d'),as.Date("20191231",'%Y%m%d'),by='day')
  dates<-seq(as.Date(date_debut,'%Y%m%d'),as.Date(date_fin,'%Y%m%d'),by='day')
  
  
  Qobs<-Qobs[which(dates %in% dates_complete)]
  
  YYYY<-year(dates[which(dates %in% dates_complete)])
  MM<-month(dates[which(dates %in% dates_complete)])
  DD<-day(dates[which(dates %in% dates_complete)])
  
  df<-data.frame(YYYY,MM,DD,Qobs)
  

  
  
  tests[[nom]]<-df
  
}

bvs<-names(tests)
filename2<-paste("/media/tito/TIIGE/PRSIM/obs_outaouais_harm_bassins_sup.Rdata")
save(tests, file = filename2)
library(doParallel)
 cores<-detectCores()
 cl<-makeCluster(cores[1]-8)
 registerDoParallel(cl)
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

#trace(prsim.wave,edit=TRUE)#enlever les set seed.

source('/home/tito/Documents/R coding/PrSim_modifiedLD/PRSimLD.wave.R')

start_sim_number<-as.numeric(args[1L])
print(start_sim_number)
for(i in start_sim_number:(start_sim_number+100)) {
  #require(PRSim)
  stoch_sim<-prsimLD.wave(tests, "Qobs", 10, suppWarn=TRUE)
  names(stoch_sim)<-bvs
  dir.create(paste0('/media/tito/TIIGE/PRSIM/0.9997_harm/sims_final'), showWarnings = FALSE) #stops warnings if folder already exists
  filename<-paste("/media/tito/TIIGE/PRSIM/0.9997_harm/sims_final/stoch_sim_10_outaouais_Kappa_",as.character(i),"_9997_LD.Rdata",sep='')
  
  save(stoch_sim, file = filename)
  
}


