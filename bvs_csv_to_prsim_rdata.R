#Fait par Fabian Tito Arandia Martinez


library(tidyverse)
library(lubridate)

dir_csv<-"/media/tito/TIIGE/PRSIM/ApportPRSIM/ApportHarmonise_BassinSup/"
setwd(dir_csv)

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