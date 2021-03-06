library("lme4")
library("mvtnorm")
library("MCMCpack")
library("nloptr")
library("ggplot2")
library("dplyr")
library("stats4")
############################################################################
# Empirical Bayes Estimator of beta, sigR2 and sigS2; one layer of laplace
#### Gaussian prior
# logpost11G calculate the posterior probability of model 1
# logpost21G calculate the posterior probability of model 2
# logpost31G calculate the posterior probability of model 3
# logpost41G calculate the posterior probability of model 4
###########################################################################

logpost11G<-function(D,par){
    aR<-par$aRs
    bR<-par$bRs
    aS<-par$aSs
    bS<-par$bSs
    K<-nlevels(D$Rep)
    # Design matrix
    ZR=as.matrix(model.matrix(~0+(factor(D[,4]))))
    # Fit GLMM
    mod1=glmer(formula = cbind(YI, NI-YI)~ 0+(1 |Rep),
               family = binomial, data = D) 
    bhati<-unlist(lapply(ranef(mod1),'[[',1)) 
    sigR2hat<-lapply(VarCorr(mod1),'[[',1)$Rep
    sigR2tulta<-((K-1)*sigR2hat+2*bR)/(K-1+2*aR)
    # integral out random effect b
    # bhat H1
    nfb<-function(x){
        sigmax=diag(sigR2tulta,K)
        P<-as.vector(exp(ZR%*%x)/(1+exp(ZR%*%x)))
        sum<-sum(log(dbinom(D[,3],D[,2],P)))
        return(-sum-dmvnorm(x,mean=rep(0,K),sigma=sigmax,log=TRUE))
    }
    bopt<-optim(bhati,nfb,hessian=TRUE)
    bhat<-bopt$par
    H1<-bopt$hessian
    ##plug in everything
    logbin<-0
    sigmax=diag(sigR2tulta,K)
    P<-as.vector(exp(ZR%*%bhat)/(1+exp(ZR%*%bhat)))
    logbin<-sum(log(dbinom(D[,3],D[,2],P)))
    #Final expression
    return(0.5*K*log(2*pi)-0.5*log(det(H1))+
               logbin+dmvnorm(bhat,mean=rep(0,K),sigma=sigmax,log=TRUE))
}


logpost21G<-function(D,par){
    aR<-par$aRs
    bR<-par$bRs
    aS<-par$aSs
    bS<-par$bSs
    mu<-par$mlemus
    s<-paras$mlesigmas
    K2<-nlevels(D$Rep)
    ZR2=as.matrix(model.matrix(~0+(factor(D[,4]))))
    # fit GLMM
    # betahati2 and bhati2: starting value in optim() 
    mod2=glmer(formula = cbind(YI, NI-YI)~ 1+(1 |Rep),
               family = binomial, data = D)
    betahati2<-as.vector(fixef(mod2))
    varbeta<-(summary(mod2)$coefficients[2])^2
    betatulta<-(mu*varbeta+s^2*betahati2)/(varbeta+s^2)
    bhati2<-as.vector(unlist(lapply(ranef(mod2),'[[',1)))
    sigR2hat<-lapply(VarCorr(mod2),'[[',1)$Rep
    sigR2tulta<-((K2-1)*sigR2hat+2*bR)/(K2-1+2*aR)
    sigma2f<-diag(sigR2tulta,K2)
    #integral out b 
    nffunction<-function(b) {
        P<-as.vector(exp(betatulta+ZR2%*%b)/(1+exp(betatulta+ZR2%*%b)))
        sum21<-sum(log(dbinom(D[,3],D[,2],P)))
        logpaib2<-dmvnorm(b,mean=rep(0,K2),sigma=sigma2f,log=TRUE)
        return(-sum21-logpaib2)
    }
    optAll<-optim(bhati2,nffunction,hessian=TRUE)
    H1<-optAll$hessian
    bhatf<-optAll$par
    ## plug in everything
    P<-as.vector(exp(betatulta+ZR2%*%bhatf)/(1+exp(betatulta+ZR2%*%bhatf)))
    logbin2<-sum(log(dbinom(D[,3],D[,2],P)))
    prior2f<-dnorm(betatulta,mu,s,log=TRUE)
    mutinom<-dmvnorm(bhatf,mean=rep(0,K2),sigma=sigma2f,log=TRUE)
    #final expression 
    return((K2)/2*log(2*pi)-0.5*(log(det(H1)))+logbin2+prior2f+mutinom)
} 

logpost31G<-function(D,par){
    K3<-nlevels(D$Rep)
    J3<-nlevels(D$SNP)
    aR<-par$aRs
    bR<-par$bRs
    aS<-par$aSs
    bS<-par$bSs
    # Design matrix
    ZR3=cbind(as.matrix(model.matrix(~0+(factor(D[,1])))),as.matrix(model.matrix(~0+(factor(D[,4])))))
    # Fit GLMM
    mod3=glmer(formula = cbind(YI, NI-YI)~ 0+(1 |Rep)+(1|SNP),
               family = binomial, data = D)
    bhati3<-c(lapply(ranef(mod3),'[[',1)$SNP,lapply(ranef(mod3),'[[',1)$Rep)
    sigR2hat<-lapply(VarCorr(mod3),'[[',1)$Rep
    sigS2hat<-lapply(VarCorr(mod3),'[[',1)$SNP
    sigR2tulta<-((K3-1)*sigR2hat+2*bR)/(K3-1+2*aR)
    sigS2tulta<-((J3-1)*sigS2hat+2*bS)/(J3-1+2*aS)
    sigma3f<-as.matrix(bdiag(diag(sigS2tulta,J3),diag(sigR2tulta,K3)))
    # integral out random effect b
    nfb3<-function(x){
        P<-as.vector(exp(ZR3%*%x)/(1+exp(ZR3%*%x)))
        sum31<-sum(log(dbinom(D[,3],D[,2],P)))
        logpaib31<-dmvnorm(x,mean=rep(0,(K3+J3)),sigma=sigma3f,log=TRUE)
        return(-sum31-logpaib31)
    }
    optb3<-optim(bhati3,nfb3,hessian=TRUE)
    bhat3f<-optb3$par
    H1<-optb3$hessian
    detH1<-det(H1)
    if (detH1<=0) {detH1<-0.0001}
    # Plug in everything
    const3<-0.5*(K3+J3)*log(2*pi)-0.5*log(detH1)
    P<-as.vector(exp(ZR3%*%bhat3f)/(1+exp(ZR3%*%bhat3f)))
    logbin3<-sum(log(dbinom(D[,3],D[,2],P)))
    logmutinorm<-dmvnorm(bhat3f,mean=rep(0,(J3+K3)),sigma=sigma3f,log=TRUE)
    # final expression of function logpost3
    return(const3+logbin3+logmutinorm)
}

logpost41G<-function(D,par){
    aR<-par$aRs
    bR<-par$bRs
    aS<-par$aSs
    bS<-par$bSs
    mu<-par$mlemus
    s<-par$mlesigmas
    K4<-nlevels(D$Rep)
    J4<-nlevels(D$SNP)
    # Fit GLMM
    mod4=glmer(formula = cbind(YI, NI-YI)~ 1+(1 |Rep)+(1|SNP),
               family = binomial, data =D)
    bhati4<-c(lapply(ranef(mod4),'[[',1)$SNP,lapply(ranef(mod4),'[[',1)$Rep)
    betahati4<-as.vector(fixef(mod4))
    varbeta<-(summary(mod4)$coefficients[2])^2
    betatulta<-(mu*varbeta+s^2*betahati4)/(varbeta+s^2)
    sigR2hat<-lapply(VarCorr(mod4),'[[',1)$Rep
    sigS2hat<-lapply(VarCorr(mod4),'[[',1)$SNP
    sigR2tulta<-((K4-1)*sigR2hat+2*bR)/(K4-1+2*aR)
    sigS2tulta<-((J4-1)*sigS2hat+2*bS)/(J4-1+2*aS)
    sigma4f<-as.matrix(bdiag(diag(sigS2tulta,J4),diag(sigR2tulta,K4)))
    #Design matrix
    ZR4=cbind(as.matrix(model.matrix(~0+(factor(D[,1])))),as.matrix(model.matrix(~0+(factor(D[,4])))))
    #### integral out b 
    nffunction<-function(b) {
        P<-as.vector(exp(betatulta+ZR4%*%b)/(1+exp(betatulta+ZR4%*%b)))
        sum41<-sum(log(dbinom(D[,3],D[,2],P)))
        logpaib4<-dmvnorm(b,mean=rep(0,(K4+J4)),sigma=sigma4f,log=TRUE)
        return(-sum41-logpaib4)
    }
    optAll<-optim(bhati4,nffunction,hessian=TRUE)
    H1<-optAll$hessian
    detH1<-det(H1)
    if (detH1<=0) {detH1<-0.0001}
    bhat4f<-optAll$par
    ### plug in everyting
    const4f<-0.5*(K4+J4)*log(2*pi)-0.5*log(detH1)
    P<-as.vector(exp(betatulta+ZR4%*%bhat4f)/(1+exp(betatulta+ZR4%*%bhat4f)))
    logbin4<-sum(log(dbinom(D[,3],D[,2],P)))
    prior4f<-dnorm(betatulta,mu,s,log=TRUE)
    otherpartf<-dmvnorm(bhat4f,mean=rep(0,(K4+J4)),sigma=sigma4f,log=TRUE)
    return(const4f+logbin4+prior4f+otherpartf)
}

PPsFUN<-function(par,index,data){
    print(index)
    D<-GDD(index,data)
    return(tryCatch(matrix(c(index,
                             logpost11G(D,par=par),
                             logpost21G(D,par=par),
                             logpost31G(D,par=par),
                             logpost41G(D,par=par)
    ),nrow=1,ncol=5), 
    error=function(e){cat("error: ",conditionMessage(e), "\n")}))
}
