#Fait par: Fabian Tito Arandia Martinez
#Date de creation: 17 mars 2020
#Objectif: creer les fichiers txt a partir des sorties de PRSIM afin de les rendre 
#ingerables par Spark (calcul probabilites) et Modsim (laminage). Deux alternatives pour gagner en vitesse:sparkly r ou les librairies de parallelisation de R


library(reshape2)
library(readr)
library(parallel)
library(MASS)
library(foreach)
library(doParallel)

rm(list=ls())

numCores <- detectCores()
numCores

registerDoParallel(numCores) 

#sc <- spark_connect(master = "local")
path<-'/media/tito/TIIGE/PRSIM/0.9997/'
setwd(path)
fichiers<-list.files(paste0(path,'sims_final/'))
numeros_r<-c("r1","r2","r3","r4","r5","r6","r7","r8","r9","r10")

n_fichiers<-length(fichiers)
i=0
foreach (n=1:n_fichiers) %do% {
  fichier<-fichiers[n]
  try(rm(stoch_sim))
  i=i+1
  load(paste(path,'sims_final/',fichier,sep=""))
  numero_gen<-as.character(i)
  bvs<-names(stoch_sim)
  mainDir<-paste0(path,'bv_csv/')
  for(bv in bvs){
    test<-stoch_sim[bv][[1]]$simulation
    
    #length(test$r1)#24455/365 = 67 annees
    #unique(test$YYYY)
    
    
    
    
    #jours juliens
    tmp <- as.Date(test$timestamp, format = "%y-%m-%d")
    
    
    test['julian_day']<-format(tmp, "%j")
    test['gen_number']<-rep(i,nrow(test))
    
    test$season <- "winter"
    test$season[which(test$MM%in%c(3,4,5))] <- "spring"
    test$season[which(test$MM%in%c(6,7,8))] <- "summer"
    test$season[which(test$MM%in%c(9,10,11))] <- "fall"
    
    dir.create(file.path(mainDir, bv), showWarnings = FALSE)
    setwd(file.path(mainDir, bv))
    #melt matrix
    for(numero_r in numeros_r){
      
      test2<-melt(data=test,id.vars=c("julian_day","YYYY","season","gen_number"),measure.vars=numero_r)#c("r1","r2","r3","r4","r5","r6","r7","r8","r9","r10")
      
      #final<-dcast(test2,formula = julian_day+season+gen_number+sim_number~YYYY)
      filename<-paste(mainDir,'/',bv,'/',numero_gen,'-',bv,'-',numero_r,'.csv',sep='')
      write_csv(test2,filename)
      
      #filename<-paste(mainDir,'/',bv,'/',numero_gen,'-',bv,'-',numero_r,'-parquet',sep='')
      
      #final <- copy_to(sc, final,overwrite = TRUE)
      #Parquet = spark_write_parquet(final, filename, mode = "overwrite")
      
      
    }
    
  }
  
}
