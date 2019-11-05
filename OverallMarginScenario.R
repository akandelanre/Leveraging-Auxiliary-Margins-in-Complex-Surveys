###########################################################################
###########################################################################
################################## START ##################################
###########################################################################
###########################################################################

###### 1: Load Libraries and Funtions
rm(list = ls())
#set.seed(100)

library(xtable)
source("Functions.R")

N <- 20000;
alpha <- list(alpha0=c(0.15),alpha1=c(0,-0.45,-0.15)) #X1|Y1
gamma <- list(gamma0=c(-0.25),gamma1=c(0,0.1,0.3),gamma2=c(-1.1)) #Nonresponse R|Y1,X1
n_star <- c(600,2400)
n_iter <- 100; burn_in <- n_iter/2; nImputations <- 50

RESULTS <- MultipleRuns(N,n_star,alpha,gamma,n_iter,burn_in,nImputations,nRuns=5)
###### Scenario 0: True results and results from original sample.
###### Scenario 1: No constraint on $\hat{t}_x$ but I augmented the data with the margin for $X_2$.
###### Scenario 2: Constraint on $\hat{t}_x$ but no augmented data.
###### Scenario 3: Constraint on $\hat{t}_x$ with $w_i$ (as a factor variable). No augmented data.
###### Scenario 4: Constraint on $\hat{t}_x$ with $w_i$ (as a factor variable) with interaction. No augmented data.
###### Scenario 5: No constraint on $\hat{t}_x$ and no augmented data -- AN with $w_i$ (as a factor variable).
###### Scenario 6: No constraint on $\hat{t}_x$ and no augmented data -- MAR with $w_i$ (as a factor variable).
###### Scenario 7: Constraint on $\hat{t}_x$ and augmented data.