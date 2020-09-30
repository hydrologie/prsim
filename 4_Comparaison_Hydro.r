
library(readxl)
library(tidyverse)



#' Comparaison de la somme des apports intermediaires avec les simulations de prsim
#'
#' @param path 
#'
#' @return
#' @export
#'
#' @examples
#' path<-'/home/tito/Documents/github/prsim/outaouais_sup_lynda/Frequentielles/Apports Latéraux Methode extrême PRT/Q10000/'
#' Comparaison_Temiscamingue(path)

Comparaison_Temiscamingue<-function(path){
  
  
  setwd(path)
  
  files<-list.files()
  
  
  total_temiscamingue<-list()
  
  for(file in files){
    
    df<-read_csv(file)
    
    #somme jusqu'a temiscamingue
    
    bvs_temiscamingue<-c( "Dozois","Lac Victoria et lac Granet","Rapide 7","Rapide 2","Riviere Kinojevis","Lac des Quinze","Mistinikon","Lady Evelyn","Lower Notch et Indian Chute","Rabbit Lake", "Kipawa",
                          "Lac Temiscamingue a Angliers","Riviere Blanche" )
    
    x<-df[bvs_temiscamingue]
    
    sum_temiscamingue<-rowSums (x, na.rm = FALSE, dims = 1)
    
    sum_carillon<-rowSums (df, na.rm = FALSE, dims = 1)
    
    total_temiscamingue[[file]]<-sum_temiscamingue
  }
  
  df_test<-do.call(cbind,total_temiscamingue)
  
  matplot(df_test,type='l',col='blue',lty=1)
  grid()
  return(sum_temiscamingue)
  
}

path<-'/home/tito/Documents/github/prsim/outaouais_sup_lynda/Frequentielles/Apports Latéraux Methode extrême PRT/Q10000/'
