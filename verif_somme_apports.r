#Fait par: Fabian Tito Arandia Martinez
#Date de création: 16 june 2020

library(tidyverse)
library(readxl)

#aller lire les fichiers dans le repertoire pour hecressim
setwd('/media/tito/TIIGE/PRSIM/0.9995/bv_csv_hecressim')
fichiers<-list.files()

total<-list()
for(fichier in fichiers){ 
  df<-read_csv(fichier)
  carillon_row_sum<-rowSums(df)
  plot(carillon_row_sum,type='l')
  total[[fichier]]<-carillon_row_sum
  }

df <- do.call("cbind", total)

mean_hydro<-rowMeans(df)
max_hydro<- apply(df, 1, max) 
min_hydro<- apply(df, 1, min) 

#ecdf pointes printanières

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
colMax <- function(data) sapply(data, max, na.rm = TRUE)
#res=df%>%summarize(max=max(df,na.rm = TRUE)) %>%collect()

df_max_col<-colMax(df)
Fn<- ecdf_cunnane(df_max_col)

quantile(Fn, prob=(1-(1/10000)))
