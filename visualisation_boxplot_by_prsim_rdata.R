visualisation_boxplot <- function(dir_analysis){
### visualize simulations for all test stations
### store stochastically simulated time series
setwd(dir_analysis)
pdf('test_diff_marginals_test.pdf',width=8,height=8)
marginal_list=liste_des_distributions_choisies
par(mfrow = c(4,2),mar=c(2,3,2,1))
### station l
for(l in 1:length(data)){
  sim <- out[[l]]$simulation
  title= paste(marginal_list[[l]], "_", names(data[l]))
  ### plot example of simulated time series
  par(mai=c(.9,.9,.1,.1))
  ### observed time series
  sz = length(sim$Qobs)
  plot(sim$timestamp[1:sz], (sim$Qobs[1:sz]), type="l", main=title,
       xlab="Time [d]", ylab=expression(paste("Discharge [m"^3,"/s]")))
  ### add simulations
  matlines(sim$timestamp[1:sz], (sim[1:sz, grep("r", names(sim))]),
           lty=1, col="gray")
  
  ### compare distributions without outliers
  ### without outliers
  # boxplot(sim$Qobs,sim$r1,outline=F,col=c('black','grey'),names=c('Obs','Sim'))
  ### with outliers
  boxplot((sim$Qobs),(sim$r1),col=c('black','grey'),names=c('Obs','Sim'))
}
dev.off()

}