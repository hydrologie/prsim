source('prsim_wave_marg.R')
source('visualisation_boxplot_by_prsim_rdata.R')

generateur_prsim <- function(start_number,final_number,dir_main,dir_analysis){
#start_sim_number=1
#generateur de simulations prsim
for(i in start_sim_number:final_number) {
  out <- prsim.wave.marg(data=runoff_multi_sites, number_sim=10, 
                         marginal_list=liste_des_distributions_choisies,
                         n_par_list=nbre_de_parametres_par_distribution)#probleme avec goftest ks_test ks.test
  
  ### Save the simulations
  names(out)<-names(runoff_multi_sites)
  dir.create(paste0(dir_main, 'sims_final/'), showWarnings = FALSE) #stops warnings if folder already exists
  filename<-paste(dir_main, "sims_final/stoch_sim_10_outaouais_",as.character(i),"_MB.Rdata",sep='')
  # 
  save(out, file = filename)
  
  #ajouter visualisation
  source(visualisation_boxplot_by_prsim_rdata)
  visualisation_boxplot(dir_analysis) 
}}