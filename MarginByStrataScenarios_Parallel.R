###########################################################################
###########################################################################
################################## START ##################################
###########################################################################
###########################################################################

###### 1: Load Libraries and Funtions
rm(list = ls())
set.seed(03252019)
#library(foreach)
library(doParallel)
source("Functions_ByStrata_Parallel.R")

N <- 50000;  n_strat <- c(0.70,0.30)*N
n_star <- c(1500,3500)
#### Scennario one in thesis
#alpha <- list(alpha0=c(0.5),alpha1=c(0,-0.5,-1)) #X1|Y1
#gamma <- list(gamma0=c(-0.25),gamma1=c(0,0.1,0.3),gamma2=c(-1.1)) #Nonresponse R|Y1,X1
#### Scennario two in thesis
#alpha <- list(alpha0=c(0.5),alpha1=c(0,-0.5,-1)) #X1|Y1
#gamma <- list(gamma0=c(-1),gamma1=c(0,-0.6,1.4),gamma2=c(-0.2)) #Nonresponse R|Y1,X1
#### Scenario three in thesis
#alpha <- list(alpha0=c(0.15),alpha1=c(0,-0.45,-0.15)) #X1|Y1
#gamma <- list(gamma0=c(-0.25),gamma1=c(0,0.1,0.3),gamma2=c(-1.1)) #Nonresponse R|Y1,X1
#### Scenario four in thesis
alpha <- list(alpha0=c(0.15),alpha1=c(0,-0.45,-0.15)) #X1|Y1
gamma <- list(gamma0=c(-1),gamma1=c(0,-0.6,1.4),gamma2=c(-0.2)) #Nonresponse R|Y1,X1

n_iter <- 10000; burn_in <- n_iter/2; nImputations <- 50; nRuns <- 10
no_cores <- detectCores()# - 1
cl <- makeCluster(no_cores)
registerDoParallel(cl)
RESULTS <- foreach(1:nRuns,.combine = 'CombineResults',.inorder=FALSE,
                   .packages = c(library(sampling),library(survey),library(mvtnorm),
                                 library(DirichletReg),library(nnet),library(coda)))  %dopar%  
  SingleRun(N,n_strat,n_star,alpha,gamma,n_iter,burn_in,nImputations)
stopCluster(cl)
#getDoParWorkers()
registerDoSEQ()



###### Scenario 0: Population and no missing data
lapply(RESULTS$Results0[c("TX_True","TX_hat_OrigData")], function(x) round(colMeans(x),2))
lapply(RESULTS$Results0[c("SE_TX_hat_OrigData")], function(x) round(sqrt(colMeans(x^2)),2))
###### Scenario 6: MAR+W
lapply(RESULTS$Results6[c("TX_hat_Imp","Alpha_hat_Imp","Gamma_hat_Imp")], function(x) round(colMeans(x),2))
lapply(RESULTS$Results6[c("SE_TX_hat_Imp","SE_Alpha_hat_Imp","SE_Gamma_hat_Imp")], function(x) round(sqrt(colMeans(x^2)),2))
###### Scenario 5: AN+W
lapply(RESULTS$Results5[c("TX_hat_Imp","Alpha_hat_Imp","Gamma_hat_Imp")], function(x) round(colMeans(x),2))
lapply(RESULTS$Results5[c("SE_TX_hat_Imp","SE_Alpha_hat_Imp","SE_Gamma_hat_Imp")], function(x) round(sqrt(colMeans(x^2)),2))
###### Scenario 2: AN+C
lapply(RESULTS$Results2[c("TX_hat_Imp","Alpha_hat_Imp","Gamma_hat_Imp")], function(x) round(colMeans(x),2))
lapply(RESULTS$Results2[c("SE_TX_hat_Imp","SE_Alpha_hat_Imp","SE_Gamma_hat_Imp")], function(x) round(sqrt(colMeans(x^2)),2))
round(mean(RESULTS$Results2$AccRatio),2)
range(RESULTS$Results2$AccRatio)
###### Scenario 3: AN+C+W
lapply(RESULTS$Results3[c("TX_hat_Imp","Alpha_hat_Imp","Gamma_hat_Imp")], function(x) round(colMeans(x),2))
lapply(RESULTS$Results3[c("SE_TX_hat_Imp","SE_Alpha_hat_Imp","SE_Gamma_hat_Imp")], function(x) round(sqrt(colMeans(x^2)),2))
round(mean(RESULTS$Results2$AccRatio),2)
range(RESULTS$Results3$AccRatio)
###### Scenario 1: AN+Margin
#lapply(RESULTS$Results1[c("TX_hat_Imp","Alpha_hat_Imp","Gamma_hat_Imp")], function(x) round(colMeans(x),2))
#lapply(RESULTS$Results1[c("SE_TX_hat_Imp","SE_Alpha_hat_Imp","SE_Gamma_hat_Imp")], function(x) round(sqrt(colMeans(x^2)),2))

sd(RESULTS$Results2$TX_hat_Imp)
sd(RESULTS$Results3$TX_hat_Imp)
apply(RESULTS$Results2$Alpha_hat_Imp,2,sd)
apply(RESULTS$Results3$Alpha_hat_Imp,2,sd)
apply(RESULTS$Results2$Gamma_hat_Imp,2,sd)
apply(RESULTS$Results3$Gamma_hat_Imp,2,sd)


