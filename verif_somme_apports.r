#Fait par: Fabian Tito Arandia Martinez
#Date de création: 16 june 2020
#

library(tidyverse)
library(readxl)
rm(list=ls())
#aller lire les fichiers dans le repertoire pour hecressim
mainDir<-'/media/tito/TIIGE/PRSIM/0.9995/bv_csv_hecressim'
setwd(mainDir)


total<-list()

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
testo %>% filter(carillon_sum > 2000)

df_mean_per_julian_day = testo %>% group_by(julian_day) %>% summarise(AvgQ=mean(carillon_sum))%>% collect()
df_max_per_julian_day = testo %>% group_by(julian_day) %>% summarise(MaxQ=max(carillon_sum))%>% collect()
df_min_per_julian_day = testo %>% group_by(julian_day) %>% summarise(MinQ=min(carillon_sum))%>% collect()

#mettre en ordre les statistiques sommaires de l'hydrogramme
df_mean_per_julian_day_ordered<-df_mean_per_julian_day[order(df_mean_per_julian_day$julian_day),]
df_max_per_julian_day_ordered<-df_max_per_julian_day[order(df_max_per_julian_day$julian_day),]
df_min_per_julian_day_ordered<-df_min_per_julian_day[order(df_min_per_julian_day$julian_day),]

#ggplot des statistiques sommaires des simulations PRSIM
plot(df_max_per_julian_day_ordered,type='l',ylim=c(0,20000))                
points(df_mean_per_julian_day_ordered,type='l')
points(df_min_per_julian_day_ordered,type='l')
#sdf_pivot(testo, sim_number ~ carillon_sum)



spark_disconnect(sc)

# 
# df <- do.call("cbind", total)
# 
# mean_hydro<-rowMeans(df)
# max_hydro<- apply(df, 1, max) 
# min_hydro<- apply(df, 1, min) 
# 
# #ecdf pointes printanières
# 
# #ecdf cunnane
# ecdf_cunnane<-function (x) 
# {
#   x <- sort(x)
#   n <- length(x)
#   if (n < 1) 
#     stop("'x' must have 1 or more non-missing values")
#   vals <- unique(x)
#   rval <- approxfun(vals, cumsum(tabulate(match(x, vals))-0.4)/(n+0.2), 
#                     method = "constant", yleft = 0, yright = 1, f = 0, ties = "ordered")
#   class(rval) <- c("ecdf", "stepfun", class(rval))
#   assign("nobs", n, envir = environment(rval))
#   attr(rval, "call") <- sys.call()
#   rval
# }
# colMax <- function(data) sapply(data, max, na.rm = TRUE)
# #res=df%>%summarize(max=max(df,na.rm = TRUE)) %>%collect()
# 
# df_max_col<-colMax(df)
# Fn<- ecdf_cunnane(df_max_col)
# 
# quantile(Fn, prob=(1-(1/10000)))
