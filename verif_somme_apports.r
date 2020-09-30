#Fait par: Fabian Tito Arandia Martinez
#Date de création: 16 june 2020
#

library(tidyverse)
library(readxl)
library(reshape2)
rm(list=ls())

#historique
load('/media/tito/TIIGE/PRSIM/obs_outaouais.Rdata')#tests

qobs_total<-list()
for(bv in names(tests)){
  qobs_total[[bv]]<-tests[[bv]]$Qobs
  
}

qobs_df<-do.call(cbind,qobs_total)#55 annees de donnees
carillon_sum_obs<-rowSums(qobs_df)#cas pour carillon
q_obs_df_final<-as.data.frame(cbind(tests[[bv]]$YYYY,tests[[bv]]$MM,tests[[bv]]$DD,carillon_sum_obs))
colnames(q_obs_df_final)<-c('YYYY','MM','DD','qobs')



qobs_mean_per_day = q_obs_df_final %>% group_by(MM,DD) %>% summarise(AvgQ=mean(qobs))
qobs_max_per_day = q_obs_df_final %>% group_by(MM,DD) %>% summarise(MaxQ=max(qobs))
qobs_min_per_day = q_obs_df_final %>% group_by(MM,DD) %>% summarise(MinQ=min(qobs))

#graphiques des hydrogrammes observes chacun des bassins versants
bvs_names<-colnames(qobs_df)
qobs_df_outaouais<-as.data.frame(YYYY=cbind(tests[[bv]]$YYYY,MM=tests[[bv]]$MM,DD=tests[[bv]]$DD,qobs_df))
colnames(qobs_df_outaouais)[1]<-'YYYY'
colnames(qobs_df_outaouais)[2]<-'MM'
colnames(qobs_df_outaouais)[3]<-'DD'

#organiser les sorties par bassin versant avant de rentrer dans le calcul des statistiques sommaires
qobs_df_outaouais_melt<-melt(data = qobs_df_outaouais, id.vars = c("YYYY", "MM","DD"))


qobs_mean_per_day = qobs_df_outaouais_melt %>% group_by(MM,DD,variable) %>% summarise(AvgQ=mean(value))
qobs_min_per_day = qobs_df_outaouais_melt %>% group_by(MM,DD,variable) %>% summarise(MinQ=min(value))
qobs_max_per_day = qobs_df_outaouais_melt %>% group_by(MM,DD,variable) %>% summarise(MaxQ=max(value))



#exemple pour un bassin en particulier
bvs_names<-unique(qobs_mean_per_day$variable)
#bvs_names[] <- lapply(bvs_names, as.character)

for(bv_name in bvs_names){
  qobs_mean_per_bassin_per_day<-qobs_mean_per_day[qobs_mean_per_day$variable==bv_name,]
  qobs_min_per_bassin_per_day<-qobs_min_per_day[qobs_mean_per_day$variable==bv_name,]
  qobs_max_per_bassin_per_day<-qobs_max_per_day[qobs_mean_per_day$variable==bv_name,]
  
  
  #ggplot des statistiques sommaires des observations
  lim_max<-max(qobs_max_per_bassin_per_day$MaxQ)
  plot(qobs_mean_per_bassin_per_day$AvgQ,type='l',ylim=c(0,lim_max),main = bv_name,xlab='Jours juliens',ylab='Q')                
  points(qobs_max_per_bassin_per_day$MaxQ,type='l')
  points(qobs_min_per_bassin_per_day$MinQ,type='l')
  grid()

  
}


##########################################################
#Preparation des fichiers pour spark
#aller lire les fichiers dans le repertoire pour hecressim
mainDir<-'/media/tito/TIIGE/PRSIM/0.9995/bv_csv_hecressim'
setwd(mainDir)




mainDir2<-'/media/tito/TIIGE/PRSIM/0.9995/'
subDir<-'carillon_sum'
dir.create(file.path(mainDir2, subDir), showWarnings = FALSE)

fichiers<-list.files()

for(fichier in fichiers){ 
  df<-read_csv(fichier)
  carillon_row_sum<-rowSums(df)
  plot(carillon_row_sum,type='l')
  filename<-paste0(mainDir2,subDir,'/carillon_sum_',fichier)
  df<-as.data.frame(cbind(rep(fichier,length(carillon_row_sum)),seq(1,length(carillon_row_sum)),carillon_row_sum))
  colnames(df)<-c('sim_number','julian_day','carillon_sum')
  write_csv(as.data.frame(df),filename)
}

#lecture des csv des hydrogrammes a carillon
#strategie: 1) avoir les quantiles de l'hydrogramme final 2) calculer les quantiles de pointes

library(sparklyr)
library(dplyr)
library(tidyr)

config <- spark_config()

config$`sparklyr.shell.driver-memory` <- "2G"
config$`sparklyr.shell.executor-memory` <- "2G"
config$`spark.yarn.executor.memoryOverhead` <- "512"

# Connect to local cluster with custom configuration
sc <- spark_connect(master = "local", config = config)

spec_with_r <- sapply(read.csv('/media/tito/TIIGE/PRSIM/0.9995/carillon_sum/carillon_sum_0000001.csv', nrows = 1), class)

testo<-spark_read_csv(sc = sc,path = paste0(mainDir2,subDir),columns=spec_with_r,memory = FALSE)
src_tbls(sc)
testo$ops$vars
#testo %>% filter(carillon_sum > 2000)

df_mean_per_julian_day = testo %>% group_by(julian_day) %>% summarise(AvgQ=mean(carillon_sum))%>% collect()
df_max_per_julian_day = testo %>% group_by(julian_day) %>% summarise(MaxQ=max(carillon_sum))%>% collect()
df_min_per_julian_day = testo %>% group_by(julian_day) %>% summarise(MinQ=min(carillon_sum))%>% collect()

df_max_per_julian_day = testo %>% summarise(MaxQ=max(carillon_sum))%>% collect()


#mettre en ordre les statistiques sommaires de l'hydrogramme
df_mean_per_julian_day_ordered<-df_mean_per_julian_day[order(df_mean_per_julian_day$julian_day),]
df_max_per_julian_day_ordered<-df_max_per_julian_day[order(df_max_per_julian_day$julian_day),]
df_min_per_julian_day_ordered<-df_min_per_julian_day[order(df_min_per_julian_day$julian_day),]

#ggplot des statistiques sommaires des simulations PRSIM
plot(df_max_per_julian_day_ordered,type='l',ylim=c(0,20000))                
points(df_mean_per_julian_day_ordered,type='l')
points(df_min_per_julian_day_ordered,type='l')
#ajout des informations de debit observes
points(qobs_mean_per_day$AvgQ,type='l',col='red')                
points(qobs_max_per_day$MaxQ,type='l',col='red')
points(qobs_min_per_day$MinQ,type='l',col='red')
#sdf_pivot(testo, sim_number ~ carillon_sum)

#calcul de la pointe a carillon
#ajouter les saisons
#res=testo%>%filter(season %in% target)%>%group_by(sim_number)%>%summarize(max=max(carillon_sum))%>%collect()

#maximun par annee
res=testo%>%group_by(sim_number)%>%summarize(max=max(carillon_sum))%>%collect()
#ecdf cunnane
ecdf_cunnane<-function (x) 
{
  x <- sort(x)
  n <- length(x)
  if (n < 1) 
    stop("'x' must have 1 or more non-missing values")
  vals <- unique(x)
  rval <- approxfun(vals, cumsum(tabulate(match(x, vals))-0.4)/(n+0.2), 
                    method = "constant", yleft = 0, yright = 1, f = 0, ties = "ordered")
  class(rval) <- c("ecdf", "stepfun", class(rval))
  assign("nobs", n, envir = environment(rval))
  attr(rval, "call") <- sys.call()
  rval
}
ecdf_max_year<-ecdf_cunnane(res$max)
plot(ecdf_max_year)
Fn<- ecdf_cunnane(res$max)
#prendre les codes pour cunnane
quantiles_qinter<-data.frame(quantiles=quantile(Fn, prob = c((1-(1/10000)),(1-(1/2000)),(1-(1/1000)),(1-(1/200)),(1-(1/100)),(1-(1/50)),(1-(1/20)),(1-(1/10)),(1-(1/2))), names = FALSE),row.names=c(10000,2000,1000,200,100,50,20,10,2))

quantiles_qinter_2<-data.frame(quantiles=quantile(ecdf_max_year, prob = c((1-(1/10000)),(1-(1/2000)),(1-(1/1000)),(1-(1/200)),(1-(1/100)),(1-(1/50)),(1-(1/20)),(1-(1/10)),(1-(1/2))), names = FALSE),row.names=c(10000,2000,1000,200,100,50,20,10,2))

#

spark_disconnect(sc)

################################################################################################
#Preparation des fichiers pour spark
#aller lire les fichiers dans le repertoire pour hecressim

#lecture des csv des hydrogrammes aux bassins versants
#strategie: 1) avoir les quantiles de l'hydrogramme final 2) calculer les quantiles de pointes

library(sparklyr)
library(dplyr)
library(tidyr)

config <- spark_config()

config$`sparklyr.shell.driver-memory` <- "2G"
config$`sparklyr.shell.executor-memory` <- "2G"
config$`spark.yarn.executor.memoryOverhead` <- "512"

# Connect to local cluster with custom configuration
sc <- spark_connect(master = "local", config = config)

spec_with_r <- sapply(read.csv('/media/tito/TIIGE/PRSIM/0.9995/bv_csv_hecressim_summary//0000001.csv', nrows = 1), class)

mainDir2<-'/media/tito/TIIGE/PRSIM/0.9995/'

subDir<-'bv_csv_hecressim_summary'
testo<-spark_read_csv(sc = sc,path = paste0(mainDir2,subDir),columns=spec_with_r,memory = FALSE)#specifier par colonne
src_tbls(sc)
testo$ops$vars
#testo %>% filter(Cabonga > 400)%>%collect()




df_mean_per_julian_day = testo %>% group_by(julian_day,variable) %>% summarise(AvgQ=mean(value))%>% collect()
save(df_mean_per_julian_day, file = "df_mean_per_julian_day_outaouais.RData")

df_max_per_julian_day = testo %>% group_by(julian_day) %>% summarise(MaxQ=max(value))%>% collect()
save(df_max_per_julian_day, file = "df_max_per_julian_day_outaouais.RData")

df_min_per_julian_day = testo %>% group_by(julian_day) %>% summarise(MinQ=min(value))%>% collect()
save(df_min_per_julian_day, file = "df_min_per_julian_day_outaouais.RData")

#mettre en ordre les statistiques sommaires de l'hydrogramme
df_mean_per_julian_day_ordered<-df_mean_per_julian_day[order(df_mean_per_julian_day$julian_day),]

df_max_per_julian_day_ordered<-df_max_per_julian_day[order(df_max_per_julian_day$julian_day),]
df_min_per_julian_day_ordered<-df_min_per_julian_day[order(df_min_per_julian_day$julian_day),]

#ggplot des statistiques sommaires des simulations PRSIM
plot(df_max_per_julian_day_ordered,type='l',ylim=c(0,20000))                
points(df_mean_per_julian_day_ordered,type='l')
points(df_min_per_julian_day_ordered,type='l')
#ajout des informations de debit observes
points(qobs_mean_per_day$AvgQ,type='l',col='red')                
points(qobs_max_per_day$MaxQ,type='l',col='red')
points(qobs_min_per_day$MinQ,type='l',col='red')
#sdf_pivot(testo, sim_number ~ carillon_sum)

#calcul de la pointe a carillon
#ajouter les saisons
#res=testo%>%filter(season %in% target)%>%group_by(sim_number)%>%summarize(max=max(carillon_sum))%>%collect()

#maximun par annee
res=testo%>%group_by(sim_number)%>%summarize(max=max(carillon_sum))%>%collect()
#ecdf cunnane
ecdf_cunnane<-function (x) 
{
  x <- sort(x)
  n <- length(x)
  if (n < 1) 
    stop("'x' must have 1 or more non-missing values")
  vals <- unique(x)
  rval <- approxfun(vals, cumsum(tabulate(match(x, vals))-0.4)/(n+0.2), 
                    method = "constant", yleft = 0, yright = 1, f = 0, ties = "ordered")
  class(rval) <- c("ecdf", "stepfun", class(rval))
  assign("nobs", n, envir = environment(rval))
  attr(rval, "call") <- sys.call()
  rval
}
ecdf_max_year<-ecdf_cunnane(res$max)
plot(ecdf_max_year)
Fn<- ecdf_cunnane(res$max)
#prendre les codes pour cunnane
quantiles_qinter<-data.frame(quantiles=quantile(Fn, prob = c((1-(1/10000)),(1-(1/2000)),(1-(1/1000)),(1-(1/200)),(1-(1/100)),(1-(1/50)),(1-(1/20)),(1-(1/10)),(1-(1/2))), names = FALSE),row.names=c(10000,2000,1000,200,100,50,20,10,2))

quantiles_qinter_2<-data.frame(quantiles=quantile(ecdf_max_year, prob = c((1-(1/10000)),(1-(1/2000)),(1-(1/1000)),(1-(1/200)),(1-(1/100)),(1-(1/50)),(1-(1/20)),(1-(1/10)),(1-(1/2))), names = FALSE),row.names=c(10000,2000,1000,200,100,50,20,10,2))

#

spark_disconnect(sc)

