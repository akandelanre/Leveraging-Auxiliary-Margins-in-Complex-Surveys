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
  Sample <- Data[,c("S","Y1","X1")]; Sample$W <- 1/Data$Prob
  Sample <- Sample[,c("S","W","Y1","X1")]
  
  
  if(GetData){
    return(Sample)
  } else {
    SampleDesign <- svydesign(ids=~1,variables=~Y1+X1,weights=~W,data=Sample)
    SampleGLM <- svyglm(X1~Y1,design=SampleDesign,family=quasibinomial(probit))
    return(c(round(coef(SampleGLM),4),
             by(Sample[,c("W","X1")],Sample$S,function(x) round(sum(x$X1%*%x$W),4))))
    
  }
}


GibbsSampler <- function(Data,Model,n_iter,burn_in,UseContraint,mean_prop,sd_prop,nImputations,Augment,nAugment,Margin){
  n <- nrow(Data)
  
  
  ## augment the data
  AugData <- Data
  if(Augment){
    n_0 <- nAugment*n
    n_0_strata <- n_0 * table(AugData$S)/n; strata_values <- unique(AugData$S)
    weight_values <- unique(AugData$W)
    for(St in 1:length(strata_values)){
      AugData <- rbind(AugData,data.frame(S=strata_values[St],W=weight_values[St],Y1=NA,
                                          X1=rbinom(n_0_strata[as.character(strata_values[St])],1,
                                                    Margin$X1[as.character(strata_values[St])]),R=NA))
    }
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
  acc_ratio <- rep(0,length(mean_prop)); names(acc_ratio) <- names(mean_prop)
  
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
        for(St in 1:length(mean_prop)){
          St_name <- names(mean_prop)[St]
          log_MH_ratio <- dnorm(((X1_prop[1:n])[Data$S==as.numeric(St_name)]%*%(AugData$W[1:n])[Data$S==as.numeric(St_name)]),
                                mean_prop[St_name],sd=sd_prop[St_name],log=TRUE) - 
            dnorm(((AugData$X1[1:n])[Data$S==as.numeric(St_name)]%*%(AugData$W[1:n])[Data$S==as.numeric(St_name)]),
                  mean_prop[St_name],sd=sd_prop[St_name],log=TRUE)  
          if(log(runif(1)) < log_MH_ratio){
            AugData$X1[AugData$S==as.numeric(St_name) & MissDataInd[,"X1"]] <- 
              X1_prop[AugData$S==as.numeric(St_name) & MissDataInd[,"X1"]]
            acc_ratio[St_name] <- acc_ratio[St_name] + 1
          }
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


AnalyzeResults <- function(Results,PopulationGLM,nImputations,R_Model){
  Alpha_Imp <- matrix(0,nrow=length(coef(PopulationGLM)),ncol=nImputations)
  Alpha_SE_Imp <- matrix(0,nrow=length(coef(PopulationGLM)),ncol=nImputations)
  Gamma_Imp <- matrix(0,nrow=ncol(Results$PARA$gamma),ncol=nImputations) #delete when using post samples directly
  Gamma_SE_Imp <- matrix(0,nrow=ncol(Results$PARA$gamma),ncol=nImputations) #delete when using post samples directly
  Tx_Imp <- matrix(0,nrow=1,ncol=nImputations)
  Tx_SE_Imp <- matrix(0,nrow=1,ncol=nImputations)
  for(mm in 1:nImputations){
    Imp_mm <- Results$ImputedData[[mm]]
    ImpDesign <- svydesign(ids=~1,variables=~Y1+X1+R,weights=~W,data=Imp_mm)
    ImpGLM <- svyglm(X1~Y1,design=ImpDesign,family=quasibinomial(probit))
    #ImpGLM_R <- svyglm(R_Model,design=ImpDesign,family=quasibinomial(probit))
    ImpGLM_R <- glm(R_Model,family=binomial(probit),data=Imp_mm) #delete when using post samples directly
    
    Alpha_Imp[,mm] <- c(round(coef(ImpGLM),4))
    Alpha_SE_Imp[,mm] <- c(round(summary(ImpGLM)$coefficient[,2],4))
    Gamma_Imp[,mm] <- c(round(coef(ImpGLM_R),4)) #delete when using post samples directly
    Gamma_SE_Imp[,mm] <- c(round(summary(ImpGLM_R)$coefficient[,2],4)) #delete when using post samples directly
    Tx_Imp[mm] <- round(coef(svytotal(~X1,design=ImpDesign)),4)
    Tx_SE_Imp[mm] <- round(SE(svytotal(~X1,design=ImpDesign)),4)
  }
  
  list(
    AccRatio=round(Results$acc_ratio,2),
    TX_hat_Imp=round(mean(c(Tx_Imp)),2),
    SE_TX_hat_Imp=round(sqrt(((1+1/nImputations)*var(c(Tx_Imp))) + 
                               mean(c(Tx_SE_Imp^2))),2),
    Alpha_hat_Imp=round(rowMeans(Alpha_Imp),2),
    SE_Alpha_hat_Imp=round(sqrt(((1+1/nImputations)*apply(Alpha_Imp,1,var)) + 
                                  apply(Alpha_SE_Imp^2,1,mean)),2),
    #Gamma_hat_Imp=round(colMeans(Results$PARA$gamma),2),
    #SE_Gamma_hat_Imp=round(apply(Results$PARA$gamma,2,sd),2)
    Gamma_hat_Imp=round(rowMeans(Gamma_Imp),2),
    SE_Gamma_hat_Imp=round(sqrt(((1+1/nImputations)*apply(Gamma_Imp,1,var)) + 
                                  apply(Gamma_SE_Imp^2,1,mean)),2)
  )
}


SingleRun <- function(N,n_strat,n_star,alpha,gamma,n_iter,burn_in,nImputations){
  RESULTS <- NULL
  
  ###### 2: Create Strata
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
  Margin <- list(X1=by(Population$X1,Population$S,function(x) sum(x)/length(x)))
  
  
  ###### 9:
  ###### Scenario 0: True results and results from original sample.
  RESULTS$Results0 <- list(TX_True=round(sum(Population$X1),2),
                           TX_hat_OrigData=round(coef(svytotal(~X1,design=CompleteDataDesign)),2),
                           SE_TX_hat_OrigData=round(SE(svytotal(~X1,design=CompleteDataDesign)),2),
                           Alpha_True=round(coef(PopulationGLM),2),
                           Alpha_hat_OrigData=round(summary(CompleteDataGLM)$coefficient[,1],2),
                           SE_Alpha_hat_OrigData=round(summary(CompleteDataGLM)$coefficient[,2],2),
                           Gamma_True=c(gamma$gamma0,gamma$gamma1[-1],gamma$gamma2))
  
  
  ###### Scenario 1: AN+Margin
  #Model <- list(X1=~Y1,R=~Y1+X1)
  #RESULTS$Results1 <- AnalyzeResults(GibbsSampler(Data,Model,n_iter,burn_in,UseContraint=FALSE,
  #                                                mean_prop=by(Population$X1,Population$S,function(x) sum(x)),
  #                                                sd_prop = HT_Estimate$SE[as.character(unique(Data$S))],
  #                                                nImputations,Augment=TRUE,nAugment=2,Margin),
  #                                   PopulationGLM,nImputations,R_Model=R~Y1+X1)
  
  
  ###### Scenario 2: AN+Constraint
  Model <- list(X1=~Y1,R=~Y1+X1)
  RESULTS$Results2 <- AnalyzeResults(GibbsSampler(Data,Model,n_iter,burn_in,UseContraint=TRUE,
                                                  mean_prop=by(Population$X1,Population$S,function(x) sum(x)),
                                                  sd_prop = HT_Estimate$SE[as.character(unique(Data$S))],
                                                  nImputations,Augment=FALSE,nAugment=2,Margin),
                                     PopulationGLM,nImputations,R_Model=R~Y1+X1)
  
  
  ###### Scenario 3: AN+Constraint+Weight
  Model <- list(X1=~Y1+as.factor(W),R=~Y1+X1)
  RESULTS$Results3 <- AnalyzeResults(GibbsSampler(Data,Model,n_iter,burn_in,UseContraint=TRUE,
                                                  mean_prop=by(Population$X1,Population$S,function(x) sum(x)),
                                                  sd_prop = HT_Estimate$SE[as.character(unique(Data$S))],
                                                  nImputations,Augment=FALSE,nAugment=2,Margin),
                                     PopulationGLM,nImputations,R_Model=R~Y1+X1)
  
  
  ###### Scenario 4: Augmented data with $w_i$ (as a factor variable). No constraint
  #Model <- list(X1=~Y1+as.factor(W),R=~Y1+X1)
  #RESULTS$Results4 <- AnalyzeResults(GibbsSampler(Data,Model,n_iter,burn_in,UseContraint=FALSE,
  #                                                mean_prop=by(Population$X1,Population$S,function(x) sum(x)),
  #                                                sd_prop = HT_Estimate$SE[as.character(unique(Data$S))],
  #                                                nImputations,Augment=TRUE,nAugment=2,Margin),
  #                                   PopulationGLM,nImputations,R_Model=R~Y1+X1)
  
  
  ###### Scenario 5: AN+Weight
  Model <- list(X1=~Y1+as.factor(W),R=~Y1+X1)
  RESULTS$Results5 <- AnalyzeResults(GibbsSampler(Data,Model,n_iter,burn_in,UseContraint=FALSE,
                                                  mean_prop=by(Population$X1,Population$S,function(x) sum(x)),
                                                  sd_prop = HT_Estimate$SE[as.character(unique(Data$S))],
                                                  nImputations,Augment=FALSE,nAugment=2,Margin),
                                     PopulationGLM,nImputations,R_Model=R~Y1+X1)
  
  
  ###### ###### Scenario 6: MAR+Weight
  Model <- list(X1=~Y1+as.factor(W),R=~Y1)
  RESULTS$Results6 <- AnalyzeResults(GibbsSampler(Data,Model,n_iter,burn_in,UseContraint=FALSE,
                                                  mean_prop=by(Population$X1,Population$S,function(x) sum(x)),
                                                  sd_prop = HT_Estimate$SE[as.character(unique(Data$S))],
                                                  nImputations,Augment=FALSE,nAugment=2,Margin),
                                     PopulationGLM,nImputations,R_Model=R~Y1)
  
  
  ###### Scenario 7: AN+Margin+Constraint
  #Model <- list(X1=~Y1,R=~Y1+X1)
  #RESULTS$Results7 <- AnalyzeResults(GibbsSampler(Data,Model,n_iter,burn_in,UseContraint=TRUE,
  #                                                mean_prop=by(Population$X1,Population$S,function(x) sum(x)),
  #                                                sd_prop = HT_Estimate$SE[as.character(unique(Data$S))],
  #                                                nImputations,Augment=TRUE,nAugment=2,Margin),
  #                                   PopulationGLM,nImputations,R_Model=R~Y1+X1)
  return(RESULTS)
}


CombineResults <- function(RESULTS1,RESULTS2){
  for(k in 1:(length(RESULTS1$Results0)-1)){
    RESULTS1$Results0[[k]] <- rbind(RESULTS1$Results0[[k]],RESULTS2$Results0[[k]])
  }
  for(j in 2:length(RESULTS1$Results5)){
    #RESULTS1$Results1[[j]] <- rbind(RESULTS1$Results1[[j]],RESULTS2$Results1[[j]])
    RESULTS1$Results5[[j]] <- rbind(RESULTS1$Results5[[j]],RESULTS2$Results5[[j]])
    RESULTS1$Results6[[j]] <- rbind(RESULTS1$Results6[[j]],RESULTS2$Results6[[j]])
  }
  for(j in 1:length(RESULTS1$Results2)){
    RESULTS1$Results2[[j]] <- rbind(RESULTS1$Results2[[j]],RESULTS2$Results2[[j]])
    RESULTS1$Results3[[j]] <- rbind(RESULTS1$Results3[[j]],RESULTS2$Results3[[j]])
    #RESULTS1$Results4[[j]] <- rbind(RESULTS1$Results4[[j]],RESULTS2$Results4[[j]])
    #RESULTS1$Results7[[j]] <- rbind(RESULTS1$Results7[[j]],RESULTS2$Results7[[j]])
  }
  
  return(RESULTS1)
}
###### End of Functions ######





