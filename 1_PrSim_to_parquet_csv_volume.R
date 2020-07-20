#Fait par: Fabian Tito Arandia Martinez
#Date de creation: 17 mars 2020
#Objectif: creer les fichiers txt a partir des sorties de PRSIM afin de les rendre 
#ingerables par Spark (calcul probabilites) et Modsim (laminage)

library(sparklyr)
library(reshape2)
library(readr)
library(zoo)
library(dplyr)
library(roll)

rm(list=ls())
#sc <- spark_connect(master = "local")
setwd('/home/tito/Desktop')
fichiers<-list.files('/media/tito/TIIGE/PRSIM/0.9995/sims/')
numeros_r<-c("r1","r2","r3","r4","r5","r6","r7","r8","r9","r10")

i=0
for(fichier in fichiers){
  try(rm(stoch_sim))
  i=i+1
  load(paste('/media/tito/TIIGE/PRSIM/0.9995/sims/',fichier,sep=""))
  numero_gen<-as.character(i)
  bvs<-names(stoch_sim)
  mainDir<-'/media/tito/TIIGE/PRSIM/0.9995/bv_csv_volume/'
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
      
      res=test2%>%filter(season %in% 'spring')%>%group_by(YYYY,gen_number,variable)%>%mutate(vol =value*(1e-6*60*60*24) )%>%
        mutate(roll_sum = roll_sum(vol, 70))  %>%summarize(max=max(roll_sum,na.rm = TRUE)) %>%collect()
      #%>% mutate(rollapply_sum =rollapplyr(vol, 70, sum,fill=NA) )
      #filename<-paste('/home/tito/Desktop/',bv,'-PRSIM-Kappa.Rdata',sep='')
      #final<-dcast(test2,formula = julian_day+season+gen_number+sim_number~YYYY)
      filename<-paste(mainDir,'/',bv,'/',numero_gen,'-',bv,'-',numero_r,'-volume_prt_max_annuel.csv',sep='')
      write_csv(res,filename)
      
      #filename<-paste(mainDir,'/',bv,'/',numero_gen,'-',bv,'-',numero_r,'-parquet',sep='')
      
      #final <- copy_to(sc, final,overwrite = TRUE)
      #Parquet = spark_write_parquet(final, filename, mode = "overwrite")
      
      
    }
    
  }
  
}
