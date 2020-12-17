#Fait par: Fabian Tito Arandia Martinez
#Date de creation: 18 mars 2020
#Objectif: lire les fichiers parquet ou csv dans spark et calculer la fonction de 
#densite cumulative empirique par bassin versant


library(sparklyr)
library(dplyr)
rm(list=ls())

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
#configuration pour les calculs effectues sur le spark dataframe
# Initialize configuration with defaults
path<-'/media/tito/TIIGE/PRSIM/0.9997/'
setwd(path)
config <- spark_config()

config$`sparklyr.shell.driver-memory` <- "2G"
config$`sparklyr.shell.executor-memory` <- "2G"
config$`spark.yarn.executor.memoryOverhead` <- "512"
config$`sparklyr.cores.local` <- 9

# Connect to local cluster with custom configuration
sc <- spark_connect(master = "local", config = config)

spec_with_r <- sapply(read.csv(paste0(path,"bv_csv_volume/Bark Lake/1-Bark Lake-r1-volume_prt_max_annuel.csv"), nrows = 10), class)

quantiles_qinter_bvs<-list()
#boucle a faire
bvs<-list.files('/media/tito/TIIGE/PRSIM/0.9997/bv_csv_volume/')

for(bv in bvs){
  testo<-spark_read_csv(sc = sc,path = paste(path,'bv_csv_volume/',bv,'/',sep=''),columns=spec_with_r,memory = FALSE)
  
  src_tbls(sc)
  #target <- c("summer", "winter")
  
  res=testo%>%collect()
  
  filename<-paste(path,'max_prt_volumes/',bv,'-PRSIM-Kappa.Rdata',sep='')
  save(res, file = filename)
  #testis<-testo%>%arrange(value)%>%collect()%trop de donnees pour lancer
  
  #il ny a que des 1950 pour 36500 entrees. Alors quil faut 36500 * 67
  #png(paste(bv,'_09993.png',sep=''))
  #plot(ecdf(testo$value),xlab = 'Q (m3/s)',main=bv)
  #dev.off()

  
  #calcul des quantiles
  
  Fn<- ecdf_cunnane(res$max)
  #prendre les codes pour cunnane
  quantiles_qinter<-data.frame(quantiles=quantile(Fn, prob = c((1-(1/10000)),(1-(1/2000)),(1-(1/1000)),(1-(1/200)),(1-(1/100)),(1-(1/50)),(1-(1/20)),(1-(1/10))), names = FALSE),row.names=c(10000,2000,1000,200,100,50,20,10))
  quantiles_qinter_bvs[[bv]]<-quantiles_qinter
}


spark_disconnect(sc)

#reorganiser la liste

df <- data.frame(matrix(unlist(quantiles_qinter_bvs), nrow=length(quantiles_qinter_bvs), byrow=T))
rownames(df)<-names(quantiles_qinter_bvs)
colnames(df)<-c(10000,2000,1000,200,100,50,20,10)
write.csv(df,'quantiles_prt_outaouais_prelim_prsim_kappaLD_printemps_volume_09997.csv')

