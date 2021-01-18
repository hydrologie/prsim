#Fait par: Fabian Tito Arandia Martinez
#date de creation: 5 juin 2020
#modification 1: 23 juin 2020 :: formattage des tableaux excel

library(tidyverse)
library(readxl)
library(openxlsx)
options(scipen = 999)
rm(list=ls())
#repertoires
prsim_rep<-'/Users/caramelo/Documents/0000_Unu Engineering/github/prsim/'

#setwd('/home/tito/Documents/R coding/PrSim/plot_verifications/Pointes/')
quantiles_pointe_obs<-readxl::read_excel(paste0(prsim_rep,'QuantilesOutaouaisBeauharnois_prt_regional.xlsx'))

prsim_file<-paste(prsim_rep,'quantiles_prt_outaouais_prelim_prsim_kappaLD_printemps_new.csv',sep='')
quantiles_pointe_prsim<-read_csv(prsim_file)

#quantiles_pointe_prsim[,2:9]
#quantiles_pointe_obs[,4:11]

#Merge a partir des noms des bassins versants
res<-merge(x=quantiles_pointe_prsim,y=quantiles_pointe_obs,by.x='X1',by.y='Bassin PRSIM', fill=-9999,all.x=TRUE)

#rapide-2 et rapide-7 manquants dans quantiles de pointe prsim.
difference<-res[2:7]-res[8:13]
diff_percentage<-(difference/res[2:7])*100
test<-c(quantiles_pointe_prsim[,1],diff_percentage)
test2<-do.call(cbind.data.frame, test)

#test2$`2000.x` <- NULL
#test2$`200.x`<- NULL
names(test2)<-c('Bassins versants','10000 (%)','1000 (%)', '100 (%)', '50 (%)', '20 (%)', '10 (%)')
wb <- loadWorkbook('gabarit.xlsx')
writeData(wb,sheet='Ecart_Pointe_Printaniere',x=test2,startCol = 1, startRow = 2, rowNames = TRUE)
#showGridLines(wb,sheet='Ecart_Pointe_Printaniere',showGridLines = TRUE)
#options("openxlsx.borderColour" = "#4F80BD")
#options("openxlsx.borderStyle" = "thin")
saveWorkbook(wb, file = "revision_quantiles_kappaLD_09995.xlsx", overwrite = TRUE)#write.csv(test2,'revision_quantiles_kappaLD_09995_test_printemps_volume.csv')

########
#Volume#
########


setwd('/home/tito/Documents/R coding/PrSim/plot_verifications/Volume/')
quantiles_pointe_obs<-readxl::read_excel('QuantilesOutaouaisBeauharnois_volumes.xlsx')
quantiles_pointe_prsim<-read_csv('quantiles_prt_outaouais_prelim_prsim_kappaLD_printemps_09995.csv')

quantiles_pointe_prsim[,2:9]
quantiles_pointe_obs[,4:11]

#Merge a partir des noms des bassins versants
res<-merge(x=quantiles_pointe_prsim,y=quantiles_pointe_obs,by.x='X1',by.y='Bassin PRSIM', fill=-9999,all.x=TRUE)

#rapide-2 et rapide-7 manquants dans quantiles de pointe prsim.
difference<-res[2:9]-res[12:19]
diff_percentage<-(difference/res[2:9])*100
test<-c(quantiles_pointe_prsim[,1],diff_percentage)
test2<-do.call(cbind.data.frame, test)

test2$`2000.x` <- NULL
test2$`200.x`<- NULL
names(test2)<-c('Bassins versants','10000 (%)','1000 (%)', '100 (%)', '50 (%)', '20 (%)', '10 (%)')
wb <- loadWorkbook('revision_quantiles_kappaLD_09995_test_printemps_volume.xlsx')
writeData(wb,sheet='Ecart_Volume_Printanier_70j',x=test2,startCol = 1, startRow = 2, rowNames = TRUE)
showGridLines(wb,sheet='Ecart_Volume_Printanier_70j',showGridLines = TRUE)
options("openxlsx.borderColour" = "#4F80BD")
options("openxlsx.borderStyle" = "thin")
saveWorkbook(wb, file = "revision_quantiles_kappaLD_09995_test_printemps_volume2.xlsx", overwrite = TRUE)#write.csv(test2,'revision_quantiles_kappaLD_09995_test_printemps_volume.csv')


##############################################
#Revision des quantiles d'ete-automne a faire#
##############################################

library(tidyverse)
library(readxl)
library(openxlsx)
options(scipen = 999)
setwd('/home/tito/Documents/R coding/PrSim/plot_verifications/Pointes/')
quantiles_pointe_obs<-readxl::read_excel('QuantilesOutaouaisBeauharnois_ete.xlsx')
quantiles_pointe_prsim<-read_csv('quantiles_prt_outaouais_prelim_prsim_kappaLD_ete_09995.csv')

quantiles_pointe_prsim[,2:9]
quantiles_pointe_obs[,4:11]

#Merge a partir des noms des bassins versants
res<-merge(x=quantiles_pointe_prsim,y=quantiles_pointe_obs,by.x='X1',by.y='Bassin PRSIM', fill=-9999,all.x=TRUE)

#rapide-2 et rapide-7 manquants dans quantiles de pointe prsim.
difference<-res[2:9]-res[12:19]
diff_percentage<-(difference/res[2:9])*100
test<-c(quantiles_pointe_prsim[,1],diff_percentage)
test2<-do.call(cbind.data.frame, test)

test2$`2000.x` <- NULL
test2$`200.x`<- NULL
names(test2)<-c('Bassins versants','10000 (%)','1000 (%)', '100 (%)', '50 (%)', '20 (%)', '10 (%)')
#wb <- loadWorkbook('revision_quantiles_kappaLD_09995.xlsx')
writeData(wb,sheet='Ecart_Pointe_Ete_Automne',x=test2,startCol = 1, startRow = 2, rowNames = TRUE)
#showGridLines(wb,sheet='Ecart_Pointe_Ete_Automne',showGridLines = TRUE)
#options("openxlsx.borderColour" = "#4F80BD")
#options("openxlsx.borderStyle" = "thin")
saveWorkbook(wb, file = "revision_quantiles_kappaLD_09995.xlsx", overwrite = TRUE)#write.csv(test2,'revision_quantiles_kappaLD_09995_test_printemps_volume.csv')

