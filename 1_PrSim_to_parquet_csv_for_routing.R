#Fait par: Fabian Tito Arandia Martinez
#Date de creation: 5 juin 2020
#Objectif: creer les fichiers txt a partir des sorties de PRSIM afin de les rendre 
#ingerables par hecressim ou riverware

prsim_rdata_to_hecressim_csv<-function(path){
  
library(sparklyr)
library(reshape2)
library(readr)
library(stringr)
options(scipen = 999)
rm(list=ls())


path<-'/media/tito/TIIGE/PRSIM/0.9997_harm/'
fichiers<-list.files(paste0(path,'sims_final/'))
numeros_r<-c("r1","r2","r3","r4","r5","r6","r7","r8","r9","r10")

#Preallouer la matrice avec les annees
load(paste0(path,"sims_final/",fichiers[1]))
pre_test<-stoch_sim['Bark Lake'][[1]]$simulation

i=0
for(fichier in fichiers){
  try(rm(stoch_sim))
  
  load(paste(path,'sims_final/',fichier,sep=""))
  bvs<-names(stoch_sim)
  mainDir<-paste0(path,'bv_csv_hecressim')
  for(numero_r in numeros_r){
    conc_total<-pre_test$YYYY
    
    for(bv in bvs){
      test<-stoch_sim[bv][[1]]$simulation
      #conc_total<-test$YYYY
      
      conc<-test[numero_r]
      colnames(conc)<-bv
      conc_total<-cbind(conc_total,conc)
    }
    for(year in unique(conc_total$conc_total)){
      hydros_per_sim_per_year<-conc_total[conc_total$conc_total==year,]
      i=i+1
      filename<-paste(path,'bv_csv_hecressim/',str_pad(i, 7, pad = "0"),'.csv',sep='')
      hydros_per_sim_per_year <- hydros_per_sim_per_year[,-1] 
      write.csv(x = hydros_per_sim_per_year[complete.cases(hydros_per_sim_per_year),],filename,row.names = FALSE,quote = FALSE)
    }
  }
  }}
    
