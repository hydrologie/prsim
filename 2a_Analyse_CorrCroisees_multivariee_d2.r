#Fait par: Fabian Tito Arandia Martinez
#Date de cr?ation: 2020-02-21
#Objectif: Analyse des corr?lations.

##################################
#Analyse des s?ries stochastiques#            
##################################

rm(list=ls())
filename<-paste("H:\\Projets_communs\\2020\\Outaouais PRsim\\01_Intrants\\Example Rdata pour Leslie\\stoch_sim_10_outaouais_Kappa_1_9997_LD.Rdata")

load(filename)
#sim <- stoch_sim[[1]]$simulation #Qobs, r1 

sim<-list(stoch_sim[[1]]$simulation,stoch_sim[[2]]$simulation,stoch_sim[[3]]$simulation,stoch_sim[[4]]$simulation)



######################################################################
#Tests sur la Gatineau                                               #
#Bassins ? tenir en compte: Cabonga,Baskatong,Maniwaki,Paugan,Chelsea#
######################################################################

filename2<-paste("H:\\Projets_communs\\2020\\Outaouais PRsim\\02_Calculs\\Resultats\\obs_outaouais_harm_complet-11-2021.Rdata")
load(filename2)
bvs<-names(stoch_sim)
names(stoch_sim)<-bvs

#index<-match(c('Cabonga','Baskatong','Maniwaki','Paugan','Chelsea'),bvs)


#filename<-paste('/home/tito/Desktop/sims_kappaLD/stoch_sim_10_outaouais_Kappa_1000_LD.Rdata')
#load(filename)

sim<-list()

choix_graphe<- readline(prompt="Quel système? (1: Gatineau,2:Outaouais Supérieur,) : ")

if(choix_graphe==1){
  #attention en construisant la liste par riviere
  for(i in c('Cabonga','Baskatong','Maniwaki','Paugan','Chelsea')){  
    sim[[i]]<-stoch_sim[[i]]$simulation
  }
} else {
  
  for(i in c('Dozois','Lac Victoria et lac Granet','Rapide-7','Kipawa','Lac des Quinze','Lac Temiscamingue a Angliers')){
    sim[[i]]<-stoch_sim[[i]]$simulation
  }
  
}


###


# #define plotting colors
 col_sim<-adjustcolor("#fd8d3c",alpha=0.8)
 col_sim_tran <- adjustcolor("#fd8d3c",alpha=0.2)
 col_obs <- adjustcolor( "black", alpha.f = 0.2)

# 
# ### plot cross-correlation function
par(mfrow=c(length(sim),length(sim)),mar=c(1,1,2,2))#

### run through each station comtination
#meilleure façon d'avoir un seul xlabel y ylabel c'est de passer par ggplot2.
#https://stackoverflow.com/questions/20163877/how-to-have-a-single-xlabel-and-ylabel-for-multiple-plots-on-the-same-page
bvs<-names(sim)
for(j in 1:length(sim)){
  for(i in 1:length(sim)){
    ### ccf of observations
    data_mat <- matrix(unlist(lapply(sim, "[", , "Qobs")),ncol=length(sim))
    ccf_obs <- ccf(data_mat[,i],data_mat[,j],plot=FALSE)
    ### plot ccfs of observations
    bv_bv_ccf<-paste(bvs[i],'-',bvs[j],sep='')
    plot(ccf_obs$lag,ccf_obs$acf,col=col_obs,type="l",ylim=c(0,1),main=bv_bv_ccf)#,xaxt='n',yaxt='n'
    
    ### simulated ccf
    ### run through each simulation run
    for(r in 1:10){
      data_mat_sim <- matrix(unlist(lapply(sim, "[", , paste("r",r,sep=""))),ncol=length(sim))
      ccf_sim <- ccf(na.omit(cbind(data_mat_sim[,i],data_mat_sim[,j]))[,1],na.omit(cbind(data_mat_sim[,i],data_mat_sim[,j]))[,2],plot=FALSE)
      ### add one ccf plot per simulation run
      lines(ccf_obs$lag,ccf_sim$acf,col=col_sim)
      grid (NULL,NULL, lty = 1, col = "grey") 
    }
    ### overplot observations again
    lines(ccf_obs$lag,ccf_obs$acf,col="black",lwd=2)
  }
}

