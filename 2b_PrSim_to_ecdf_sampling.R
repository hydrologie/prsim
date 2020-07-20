#Fait par: Fabian Tito Arandia Martinez
#Date de creation: 21 avril 2020
#Objectif: fonctions cumulatives suite à l'échantillonage des simulations stochastiques

#Données:
rm(list=ls())
library(tidyverse)

path<-'/media/tito/TIIGE/PRSIM/0.9995/max_prt_pointes/'
fichiers<-list.files(path=path)

#ecdf cunnane
ecdf_cunnane<-function (x) {
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

quantiles_qinter_bvs<-list()

for(fichier in fichiers){
  load(paste('/media/tito/TIIGE/PRSIM/0.9995/max_prt_pointes/',fichier,sep=''))
  
  #Premiere methode: 
  stoch_sim_size<-dim(res)[1]
  idx<-sample(1:stoch_sim_size,10000)#choisis 10000 indices parmi toutes les indices disponibles sans remise
  sampled_res<-res[idx,]
  
  bv<-strsplit(fichier,split  ='-PRSIM-')[[1]][1]
  Fn<- ecdf_cunnane(sampled_res$max)
  quantiles_qinter<-data.frame(quantiles=quantile(Fn, prob = c((1-(1/10000)),(1-(1/2000)),(1-(1/1000)),(1-(1/200)),(1-(1/100)),(1-(1/50)),(1-(1/20)),(1-(1/10))), names = FALSE),row.names=c(10000,2000,1000,200,100,50,20,10))
  quantiles_qinter_bvs[[bv]]<-quantiles_qinter$quantiles
}


df <- do.call("rbind",quantiles_qinter_bvs) #combine all vectors into a matrix
colnames(df)<-c('10000','2000','1000','200','100','50','20','10')


write.csv(df,'pointes_09995_sampled.csv')
