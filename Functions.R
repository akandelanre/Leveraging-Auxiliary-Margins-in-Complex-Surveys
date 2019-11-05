###########################################################################
###########################################################################
################################## START ##################################
###########################################################################
###########################################################################


library(sampling)
library(survey)
library(mvtnorm)
library(DirichletReg)
library(nnet)
library(coda)


###### Functions ######
GenerateSample <- function(Population,GetData,n){
  S <- sampling:::strata(Population,stratanames="S",size=n,"srswor")
  Data <- getdata(Population,S)
  Sample <- Data[,c("Y1","X1")]; Sample$W <- 1/Data$Prob
  Sample <- Sample[,c("W","Y1","X1")]
  
  
  if(GetData){
    return(Sample)
  } else {
    SampleDesign <- svydesign(ids=~1,variables=~Y1+X1,weights=~W,data=Sample)
    SampleGLM <- svyglm(X1~Y1,design=SampleDesign,family=quasibinomial(probit))
    return(c(round(coef(SampleGLM),4),round(Sample$X1%*%Sample$W,4)))
  }
}


GibbsSampler <- function(Data,Model,n_iter,burn_in,UseContraint,mean_prop,sd_prop,nImputations,Augment,nAugment,Margin){
  n <- nrow(Data)
  
  
  ## augment the data
  AugData <- Data
  if(Augment){
    n_0 <- nAugment*n
    AugData <- rbind(AugData,data.frame(W=NA,Y1=NA,X1=rbinom(n_0,1,Margin$X1),R=NA))
  }
  
  
  ## initialize missing entries
  MissDataInd <- is.na(AugData)
  # Y1
  if(sum(MissDataInd[,"Y1"])>0){
    pr_Y1_miss <- predict(multinom(Y1~1, data = AugData[1:n,],trace=FALSE),AugData[MissDataInd[,"Y1"],],type="probs")
    Ran_unif_Y1_miss <- runif(nrow(pr_Y1_miss))
    cumul_Y1_miss <- pr_Y1_miss%*%upper.tri(diag(ncol(pr_Y1_miss)),diag=TRUE)
    AugData[MissDataInd[,"Y1"],"Y1"] <- rowSums(Ran_unif_Y1_miss>cumul_Y1_miss)+1L
  }
  # X1
  if(sum(MissDataInd[,"X1"])>0){
    AugData[MissDataInd[,"X1"],"X1"] <- rbinom(sum(MissDataInd[,"X1"]),1,
                                               predict(glm(X1~Y1,data=AugData,family=binomial(probit)),
                                                       AugData[MissDataInd[,"X1"],],type = "response"))
  }
  # R
  if(sum(MissDataInd[,"R"])>0){
    AugData[MissDataInd[,"R"],"R"] <- rbinom(sum(MissDataInd[,"R"]),1,
                                             predict(glm(R~Y1+X1,data=AugData,family=binomial(probit)),
                                                     AugData[MissDataInd[,"R"],],type = "response"))
  }

  
  ## initialize parameters
  pii <- rdirichlet(1,alpha = 1 + c(table(AugData[,"Y1"])))
  alpha <- rnorm(ncol(model.matrix(Model$X1,data=AugData)))
  gamma <- rnorm(ncol(model.matrix(Model$R,data=AugData)))
  b0_alpha <- alpha; b0_alpha[] <- 0
  sigma0_alpha <- diag(alpha); diag(sigma0_alpha) <- 1
  b0_gamma <- gamma; b0_gamma[] <- 0
  sigma0_gamma <- diag(gamma); diag(sigma0_gamma) <- 1
  Z <- data.frame(X1=rnorm(nrow(AugData)))
  Z$R <- rnorm(nrow(AugData))
  
  
  ## set mcmc parameters
  acc_ratio <- 0

  M_to_use_mc <- round(seq((burn_in +1),n_iter,length.out=nImputations))
  ImputedData <- vector("list", nImputations)
  PARA <- list(pii=NULL,alpha=NULL,gamma=NULL)
  counter <- 0
  
  
  for(mc in 1:n_iter){
    counter <- counter + 1
    
    ## sample pi, the parameter for Y1
    pii <- rdirichlet(1,alpha = 1 + c(table(AugData[,"Y1"])) )
    
    ## sample alpha, the parameters for X1|Y1
    X1_cond <- model.matrix(Model$X1,data=AugData)
    mu_alpha <- solve((t(X1_cond)%*%X1_cond)+solve(sigma0_alpha))%*%
      ((t(X1_cond)%*%Z$X1)+(solve(sigma0_alpha)%*%b0_alpha))
    sigma_alpha <- solve((t(X1_cond)%*%X1_cond)+solve(sigma0_alpha))
    alpha <- rmvnorm(1,mean=mu_alpha,sigma=sigma_alpha)
    
    ## sample Z_X1, the augmented variables for X1
    mu_Z_X1 <- X1_cond%*%t(alpha)
    U_X1 <- mu_Z_X1; U_X1[] <- 0
    U_X1[AugData$X1==0] <- runif(sum(AugData$X1==0),pnorm(-Inf-mu_Z_X1[AugData$X1==0]),pnorm(0-mu_Z_X1[AugData$X1==0]))
    U_X1[AugData$X1==1] <- runif(sum(AugData$X1==1),pnorm(0-mu_Z_X1[AugData$X1==1]),pnorm(Inf-mu_Z_X1[AugData$X1==1]))
    Z$X1 <- mu_Z_X1 + qnorm(U_X1)
    
    ## sample gamma, the parameters for R|Y1,X1
    R_cond <- model.matrix(Model$R,data=AugData)
    mu_gamma <- solve((t(R_cond)%*%R_cond)+solve(sigma0_gamma))%*%
      ((t(R_cond)%*%Z$R)+(solve(sigma0_gamma)%*%b0_gamma))
    sigma_gamma <- solve((t(R_cond)%*%R_cond)+solve(sigma0_gamma))
    gamma <- rmvnorm(1,mean=mu_gamma,sigma=sigma_gamma)
    
    ## sample Z_R, the augmented variables for R
    mu_Z_R <- R_cond%*%t(gamma)
    U_Z <- mu_Z_R; U_Z[] <- 0
    U_Z[AugData$R==0] <- runif(sum(AugData$R==0),pnorm(-Inf-mu_Z_R[AugData$R==0]),pnorm(0-mu_Z_R[AugData$R==0]))
    U_Z[AugData$R==1] <- runif(sum(AugData$R==1),pnorm(0-mu_Z_R[AugData$R==1]),pnorm(Inf-mu_Z_R[AugData$R==1]))
    Z$R <- mu_Z_R + qnorm(U_Z)
    
    
    ## sample missing X1|Y1,R
    if(sum(MissDataInd[,"X1"])>0){
      pr_X1_miss <- matrix(0,ncol=2,nrow=sum(MissDataInd[,"X1"]))
      colnames(pr_X1_miss) <- c("0","1")
      Miss_cond_X1 <- AugData[MissDataInd[,"X1"],]
      Miss_cond_X1$X1[] <- 0
      pr_X1_miss[,"0"] <- dnorm(Z$R[MissDataInd[,"X1"]],model.matrix(Model$R,data=Miss_cond_X1)%*%t(gamma),1)*
        pnorm(0,model.matrix(Model$X1,data=Miss_cond_X1)%*%t(alpha),1)
      Miss_cond_X1$X1[] <- 1
      pr_X1_miss[,"1"] <- dnorm(Z$R[MissDataInd[,"X1"]],model.matrix(Model$R,data=Miss_cond_X1)%*%t(gamma),1)*
        (1-pnorm(0,model.matrix(Model$X1,data=Miss_cond_X1)%*%t(alpha),1))
      pr_X1_miss <- pr_X1_miss/matrix(rowSums(pr_X1_miss),ncol=2,nrow=sum(MissDataInd[,"X1"]))
      Ran_unif_X1_miss <- runif(nrow(pr_X1_miss))
      cumul_X1_miss <- pr_X1_miss%*%upper.tri(diag(ncol(pr_X1_miss)),diag=TRUE)
      X1_prop <- AugData$X1
      X1_prop[MissDataInd[,"X1"]] <- rowSums(Ran_unif_X1_miss>cumul_X1_miss)
      
      if(UseContraint){
        log_MH_ratio <- dnorm((X1_prop[1:n]%*%AugData$W[1:n]),mean_prop,sd=sd_prop,log=TRUE) - 
          dnorm((AugData$X1[1:n]%*%AugData$W[1:n]),mean_prop,sd=sd_prop,log=TRUE)  
        if(log(runif(1)) < log_MH_ratio){
          AugData$X1[MissDataInd[,"X1"]] <- X1_prop[MissDataInd[,"X1"]]
          acc_ratio <- acc_ratio + 1
        }
      } else {
        AugData$X1[MissDataInd[,"X1"]] <- X1_prop[MissDataInd[,"X1"]]
      }
    }
    
    
    ## sample missing Y1|X1,R
    if(sum(MissDataInd[,"Y1"])>0){
      pr_Y1_miss <- matrix(0,ncol=length(pii),nrow=sum(MissDataInd[,"Y1"]))
      Miss_cond_Y1 <- AugData[MissDataInd[,"Y1"],]
      for(l in 1:length(pii)){
        Miss_cond_Y1$Y1[] <- l
        pr_Y1_miss[,l] <- dnorm(Z$R[MissDataInd[,"Y1"]],model.matrix(Model$R,data=Miss_cond_Y1)%*%t(gamma),1)*
          dnorm(Z$X1[MissDataInd[,"Y1"]],model.matrix(Model$X1,data=Miss_cond_Y1)%*%t(alpha),1)*pii[l]
      }
      pr_Y1_miss <- pr_Y1_miss/matrix(rowSums(pr_Y1_miss),ncol=length(pii),nrow=sum(MissDataInd[,"Y1"]))
      Ran_unif_Y1_miss <- runif(nrow(pr_Y1_miss))
      cumul_Y1_miss <- pr_Y1_miss%*%upper.tri(diag(ncol(pr_Y1_miss)),diag=TRUE)
      AugData[MissDataInd[,"Y1"],"Y1"] <- rowSums(Ran_unif_Y1_miss>cumul_Y1_miss) + 1L
    }
    
    
    ## sample missing R|Y1,X1
    if(sum(MissDataInd[,"R"])>0){
      Z_R_miss <- rnorm(sum(MissDataInd[,"R"]),model.matrix(Model$R,data=AugData[MissDataInd[,"R"],])%*%t(gamma),1)
      AugData[MissDataInd[,"R"],"R"] <- ifelse(Z_R_miss>0,1,0)
    }
    
    
    ## save
    if(mc > burn_in){
      PARA$pii <- rbind(PARA$pii,pii)
      PARA$alpha <- rbind(PARA$alpha,alpha)
      PARA$gamma <- rbind(PARA$gamma,gamma)
      if(is.element(mc,M_to_use_mc)){
        ImputedData[[which(M_to_use_mc==mc)]] <- AugData[1:n,]
      }
    }
    
  }
  
  #cat("\n")
  if(UseContraint){
    return(list(PARA=PARA,ImputedData=ImputedData,acc_ratio=acc_ratio/counter))
  } else {
    return(list(PARA=PARA,ImputedData=ImputedData,acc_ratio=NA))
  }
  
}


AnalyzeResults <- function(Results,PopulationGLM,nImputations){
  Alpha_Imp <- matrix(0,nrow=length(coef(PopulationGLM)),ncol=nImputations)
  Alpha_SE_Imp <- matrix(0,nrow=length(coef(PopulationGLM)),ncol=nImputations)
  Tx_Imp <- matrix(0,nrow=1,ncol=nImputations)
  Tx_SE_Imp <- matrix(0,nrow=1,ncol=nImputations)
  for(mm in 1:nImputations){
    Imp_mm <- Results$ImputedData[[mm]]
    ImpDesign <- svydesign(ids=~1,variables=~Y1+X1,weights=~W,data=Imp_mm)
    ImpGLM <- svyglm(X1~Y1,design=ImpDesign,family=quasibinomial(probit))
    
    Alpha_Imp[,mm] <- c(round(coef(ImpGLM),4))
    Alpha_SE_Imp[,mm] <- c(round(summary(ImpGLM)$coefficient[,2],4))
    Tx_Imp[mm] <- round(coef(svytotal(~X1,design=ImpDesign)),4)
    Tx_SE_Imp[mm] <- round(SE(svytotal(~X1,design=ImpDesign)),4)
  }
  
  #cat(
  #paste("No tuning. Acceptance ratio based on the MH proposal is: "),
  #paste(round(Results$acc_ratio,2)),"\n",
  #paste("The HT estimator of t_x from the imputed datasets is: "),
  #paste(round(mean(c(Tx_Imp)),2)),"\n",
  #paste("with standard error: "),
  #paste(round(sqrt(((1+1/nImputations)*var(c(Tx_Imp))) + mean(c(Tx_SE_Imp^2))),2)),"\n",
  #paste("For comparison, the true value of t_x is equal to "),
  #paste(round(sum(Population$X1),2)),"\n",
  #paste("and the HT estimator of t_x from the original complete data is: "),
  #paste(round(coef(svytotal(~X1,design=CompleteDataDesign)),2)),"\n",
  #paste("with standard error: "),
  #paste(round(SE(svytotal(~X1,design=CompleteDataDesign)),2)),"\n",
  #paste("The weighted estimate of alpha from the imputed datasets is (using svyglm) is: "),
  #paste(round(rowMeans(Alpha_Imp),2)),"\n",
  #paste("with standard error: "),
  #paste(round(sqrt(((1+1/nImputations)*apply(Alpha_Imp,1,var)) + apply(Alpha_SE_Imp^2,1,mean)),2)),"\n",
  #paste("For comparison, the true value of alpha is equal to: "),
  #paste(round(coef(PopulationGLM),2)),"\n",
  #paste("and the weighted estimate of alpha from the original complete data (also using svyglm) is: "),
  #paste(round(summary(CompleteDataGLM)$coefficient[,1],2)),"\n",
  #paste("with standard error: "),
  #paste(round(summary(CompleteDataGLM)$coefficient[,2],2)),"\n",
  #paste("To check identifiability, the posterior estimate of gamma is equal to: "),
  #paste(round(colMeans(Results$PARA$gamma),2)),"\n",
  #paste("with posterior standard error: "),
  #paste(round(apply(Results$PARA$gamma,2,sd),2)),"\n",
  #paste("and the true of gamma is equal to: "),
  #paste(c(gamma$gamma0,gamma$gamma1[-1],gamma$gamma2)),"\n"
  #)
  
  list(
    AccRatio=round(Results$acc_ratio,2),
    TX_hat_Imp=round(mean(c(Tx_Imp)),2),
    SE_TX_hat_Imp=round(sqrt(((1+1/nImputations)*var(c(Tx_Imp))) + 
                               mean(c(Tx_SE_Imp^2))),2),
    Alpha_hat_Imp=round(rowMeans(Alpha_Imp),2),
    SE_Alpha_hat_Imp=round(sqrt(((1+1/nImputations)*apply(Alpha_Imp,1,var)) + 
                                  apply(Alpha_SE_Imp^2,1,mean)),2),
    Gamma_hat_Imp=round(colMeans(Results$PARA$gamma),2),
    SE_Gamma_hat_Imp=round(apply(Results$PARA$gamma,2,sd),2)
       )
}


MultipleRuns <- function(N,n_star,alpha,gamma,n_iter,burn_in,nImputations,nRuns){
  RESULTS <- vector("list", length = 8)
  names(RESULTS) <-  paste("Results",c(0:7),sep="")
  RESULTS <- lapply(RESULTS,function(x) list(AccRatio=NULL,TX_hat_Imp=NULL,SE_TX_hat_Imp=NULL,
                                             Alpha_hat_Imp=NULL,SE_Alpha_hat_Imp=NULL,
                                             Gamma_hat_Imp=NULL,SE_Gamma_hat_Imp=NULL) )
  RESULTS$Results0 <- list(TX_True=NULL,TX_hat_OrigData=NULL,SE_TX_hat_OrigData=NULL,
                           Alpha_True=NULL,Alpha_hat_OrigData=NULL,SE_Alpha_hat_OrigData=NULL,
                           Gamma_True=c(gamma$gamma0,gamma$gamma1[-1],gamma$gamma2))
  
  ProgressBar <- txtProgressBar(min = 0, max = nRuns, style = 3)
  
  for(i in 1:nRuns){
    setTxtProgressBar(ProgressBar, i)
    
    ###### 2: Create Strata
    n_strat <- c(0.75,0.25)*N
    Population <- data.frame(S=rep(c(1,2),n_strat))
    
    
    ###### 3: Create Y1
    pii <- data.frame(c(0.5,0.15,0.35),c(0.1,0.45,0.45))
    colnames(pii) <- NULL; pii <- t(pii)
    pi_Y1 <- pii[Population$S,]
    Ran_unif_Y1 <- runif(nrow(pi_Y1))
    cumul_Y1 <- pi_Y1%*%upper.tri(diag(ncol(pi_Y1)),diag=TRUE)
    Population$Y1 <- factor(rowSums(Ran_unif_Y1>cumul_Y1)+1L)
    
    remove(pii); remove(pi_Y1); remove(Ran_unif_Y1); remove(cumul_Y1)
    
    
    ###### 4: Create X1
    names(alpha$alpha1) <- as.character(sort(unique(Population$Y1)))
    mu_Z_X1 <- cbind(rep(1,N)*alpha$alpha0 + alpha$alpha1[as.character(Population$Y1)])
    Z_X1 <- rnorm(N,mu_Z_X1,1)
    Population$X1 <- ifelse(Z_X1>0,1,0)
    PopulationGLM <- glm(X1~Y1,family=binomial(probit),data=Population)
    
    remove(mu_Z_X1); remove(Z_X1)
    
    
    ###### 5: Create Sample
    Data <- GenerateSample(Population,GetData=TRUE,n=n_star)
    #mean(Population$X1)
    #mean(Data$X1)
    #table(Data[,c("Y1","X1")])/sum(table(Data[,c("Y1","X1")]))
    #table(Population[,c("Y1","X1")])/sum(table(Population[,c("Y1","X1")]))
    HT_Estimate <- list(Replications=t(replicate(200,GenerateSample(Population,GetData=FALSE,n=n_star))))
    HT_Estimate$Replications <- as.data.frame(HT_Estimate$Replications)
    HT_Estimate$Mean <- colMeans(HT_Estimate$Replications)
    HT_Estimate$SE <- apply(HT_Estimate$Replications,2,sd)
    
    
    
    ###### 6: Item nonresponse
    names(gamma$gamma1) <- as.character(sort(unique(Population$Y1)))
    mu_Z_R <- cbind(rep(1,nrow(Data))*gamma$gamma0 + gamma$gamma1[as.character(Data$Y1)] + Data$X1*gamma$gamma2)
    Z_R <- rnorm(nrow(Data),mu_Z_R,1)
    Data$R <- ifelse(Z_R>0,1,0)
    #table(Data$R)/sum(table(Data$R))
    CompleteData <- Data
    Data$X1[Data$R==1] <- NA
    
    remove(mu_Z_R); remove(Z_R)
    
    
    ###### 7: Check observed table vs actual true table
    CompleteDataDesign <- svydesign(ids=~1,variables=~Y1+X1,weights=~W,data=CompleteData)
    CompleteDataGLM <- svyglm(X1~Y1,design=CompleteDataDesign,family=quasibinomial(probit))
    #CCDesign <- svydesign(ids=~1,variables = ~Y1+X1,weights=~W,data=na.omit(Data))
    #CCGLM <- svyglm(X1~Y1,design=CCDesign,family=quasibinomial(probit))
    
    
    ###### 8: Fit model and examine results
    Margin <- list(X1=sum(Population$X1)/nrow(Population))
    
    
    ###### 9:
    ###### Scenario 0: True results and results from original sample.
    Results0 <- list(TX_True=round(sum(Population$X1),2),
                     TX_hat_OrigData=round(coef(svytotal(~X1,design=CompleteDataDesign)),2),
                     SE_TX_hat_OrigData=round(SE(svytotal(~X1,design=CompleteDataDesign)),2),
                     Alpha_True=round(coef(PopulationGLM),2),
                     Alpha_hat_OrigData=round(summary(CompleteDataGLM)$coefficient[,1],2),
                     SE_Alpha_hat_OrigData=round(summary(CompleteDataGLM)$coefficient[,2],2) )
    for(k in 1:length(RESULTS$Results0)){
      RESULTS$Results0[[k]] <- rbind(RESULTS$Results0[[k]],Results0[[k]])
    }
    
    
    ###### Scenario 1: No constraint on $\hat{t}_x$ but I augmented the data with the margin for $X_2$.
    Model <- list(X1=~Y1,R=~Y1+X1)
    Results1 <- AnalyzeResults(GibbsSampler(Data,Model,n_iter,burn_in,UseContraint=FALSE,
                                            mean_prop=sum(Population$X1),
                                            sd_prop = HT_Estimate$SE[length(HT_Estimate$SE)],
                                            nImputations,Augment=TRUE,nAugment=2,Margin),
                               PopulationGLM,nImputations)
    for(j in 2:length(RESULTS$Results1)){
      RESULTS$Results1[[j]] <- rbind(RESULTS$Results1[[j]],Results1[[j]])
    }
    
    
    ###### Scenario 2: Constraint on $\hat{t}_x$ but no augmented data.
    Model <- list(X1=~Y1,R=~Y1+X1)
    Results2 <- AnalyzeResults(GibbsSampler(Data,Model,n_iter,burn_in,UseContraint=TRUE,
                                            mean_prop=sum(Population$X1),
                                            sd_prop = HT_Estimate$SE[length(HT_Estimate$SE)],
                                            nImputations,Augment=FALSE,nAugment=2,Margin),
                               PopulationGLM,nImputations)
    for(j in 1:length(RESULTS$Results2)){
      RESULTS$Results2[[j]] <- rbind(RESULTS$Results2[[j]],Results2[[j]])
    }
    
    
    ###### Scenario 3: Constraint on $\hat{t}_x$ with $w_i$ (as a factor variable). No augmented data.
    Model <- list(X1=~Y1+as.factor(W),R=~Y1+X1)
    Results3 <- AnalyzeResults(GibbsSampler(Data,Model,n_iter,burn_in,UseContraint=TRUE,
                                            mean_prop=sum(Population$X1),
                                            sd_prop = HT_Estimate$SE[length(HT_Estimate$SE)],
                                            nImputations,Augment=FALSE,nAugment=2,Margin),
                               PopulationGLM,nImputations)
    for(j in 1:length(RESULTS$Results3)){
      RESULTS$Results3[[j]] <- rbind(RESULTS$Results3[[j]],Results3[[j]])
    }
    
    
    ###### Scenario 4: Constraint on $\hat{t}_x$ with $w_i$ (as a factor variable) with interaction. No augmented data.
    Model <- list(X1=~Y1+as.factor(W)+Y1*as.factor(W),R=~Y1+X1)
    Results4 <- AnalyzeResults(GibbsSampler(Data,Model,n_iter,burn_in,UseContraint=TRUE,
                                            mean_prop=sum(Population$X1),
                                            sd_prop = HT_Estimate$SE[length(HT_Estimate$SE)],
                                            nImputations,Augment=FALSE,nAugment=2,Margin),
                               PopulationGLM,nImputations)
    for(j in 1:length(RESULTS$Results4)){
      RESULTS$Results4[[j]] <- rbind(RESULTS$Results4[[j]],Results4[[j]])
    }
    
    
    ###### Scenario 5: No constraint on $\hat{t}_x$ and no augmented data -- AN with $w_i$ (as a factor variable).
    Model <- list(X1=~Y1+as.factor(W),R=~Y1+X1)
    Results5 <- AnalyzeResults(GibbsSampler(Data,Model,n_iter,burn_in,UseContraint=FALSE,
                                            mean_prop=sum(Population$X1),
                                            sd_prop = HT_Estimate$SE[length(HT_Estimate$SE)],
                                            nImputations,Augment=FALSE,nAugment=2,Margin),
                               PopulationGLM,nImputations)
    for(j in 2:length(RESULTS$Results5)){
      RESULTS$Results5[[j]] <- rbind(RESULTS$Results5[[j]],Results5[[j]])
    }
    
    
    ###### Scenario 6: No constraint on $\hat{t}_x$ and no augmented data -- MAR with $w_i$ (as a factor variable).
    Model <- list(X1=~Y1+as.factor(W)+Y1,R=~Y1)
    Results6 <- AnalyzeResults(GibbsSampler(Data,Model,n_iter,burn_in,UseContraint=FALSE,
                                            mean_prop=sum(Population$X1),
                                            sd_prop = HT_Estimate$SE[length(HT_Estimate$SE)],
                                            nImputations,Augment=FALSE,nAugment=2,Margin),
                               PopulationGLM,nImputations)
    for(j in 2:length(RESULTS$Results6)){
      RESULTS$Results6[[j]] <- rbind(RESULTS$Results6[[j]],Results6[[j]])
    }
    
    
    ###### Scenario 7: Constraint on $\hat{t}_x$ and augmented data.
    Model <- list(X1=~Y1,R=~Y1+X1)
    Results7 <- AnalyzeResults(GibbsSampler(Data,Model,n_iter,burn_in,UseContraint=TRUE,
                                            mean_prop=sum(Population$X1),
                                            sd_prop = HT_Estimate$SE[length(HT_Estimate$SE)],
                                            nImputations,Augment=TRUE,nAugment=2,Margin),
                               PopulationGLM,nImputations)
    for(j in 1:length(RESULTS$Results7)){
      RESULTS$Results7[[j]] <- rbind(RESULTS$Results7[[j]],Results7[[j]])
    }
  }
  
  close(ProgressBar)
  cat("\n")
  return(RESULTS)
}


###### End of Functions ######





