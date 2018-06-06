###############################################################################

# GDD(i,data) is a function to re-structurize the data to a structure that    # 
# contains necessary information for analysis, i.e., SNPs, Replicates, counts #
# from maternal allel (YI) and total counts (NI). This structure is friendly  #
# to GLMM fitting in glmer() function in "lme4" package. This function can be #                      #
# easily modified to apply to other similar real data structures.             #

#i is the gene number of each gene. Take the sample dataset for instance, the #
# gene number of the first gene in this dataset is 277.                       #

# > unique(sample$gene)[1]
# [1] 277
# > GDD(277,data=sample)
#    SNP NI YI Rep
#     1  22  6   1
#     1  23 11   3
#     1  16 11   4
#     2  22  6   1
#     2  23 11   3
#     2  18 10   4
#     3  47 30   1
#     3  23 11   2
#     3  39 17   3
#     3  51 23   4

##############################################################################


GDD=function(i,data){
    GeneDat=subset(data,gene==i)
    S=length(unique(GeneDat[,3]))##number of SNP
    GD=matrix(0, ncol=4, nrow=4*S)
    GD[,1]=kronecker(c(1:S),rep(1,4))  
    GD[,4]=kronecker(rep(1,S),c(1:4))  
    for(ss in 1:S){
        for(l in 1:4){
            GD[(ss-1)*4+l,2:3]=as.matrix(GeneDat[ss,((4+(l-1)*2):(5+(l-1)*2))])
        }
    }
    ##### Remove NA from Dataset
    GDNA= matrix(GD[complete.cases(GD),],ncol=4)
    GD=matrix(GDNA[(GDNA[,2]>0),],ncol=4)
    colnames(GD)=c("SNP","NI","YI","Rep")
    GD=data.frame(GD)
    GDV=data.frame(GD)
    GDV[,1]=as.factor(GDV[,1])
    GDV[,4]=as.factor(GDV[,4])
    return(GDV)
}


###############################################################################
#para.est(data) is a function estimate hyper parameters
##############################################################################
filter<-function(data){
    geneNum<-unique(data$gene)
    clean<- rep(NA,length(geneNum))
    for (i in geneNum){
        D<-GDD(i,data)
        tryCatch({
            print(i)
            mod1=glmer(formula = cbind(YI, NI-YI)~ 0+(1 |Rep),
                       family = binomial, data = D)
            mod2=glmer(formula = cbind(YI, NI-YI)~ 1+(1 |Rep),
                       family = binomial, data = D)
            mod3=glmer(formula = cbind(YI, NI-YI)~ 0+(1 |Rep)+(1|SNP),
                       family = binomial, data = D)
            mod4=glmer(formula = cbind(YI, NI-YI)~ 1+(1 |Rep)+(1|SNP),
                       family = binomial, data =D)
            clean[i]<-i
            
        },error=function(e){cat("error: ",conditionMessage(e), "\n")})
    }
    cleanf<-na.omit(clean)
    return(cleanf) 
} 


para.est<-function(data,index){
    ###### add two obs and then fit GLMM 
    geneNum<-index
    sigR<- sigS<- beta<- rep(NA,length(geneNum))
    for (i in geneNum){
        DATA<-GDD(i,data)
        D<-data.frame(SNP=DATA$SNP,NI=DATA$NI+2,YI=DATA$YI+1,Rep=DATA$Rep)
        tryCatch({
            print(i)
            mod1=glmer(formula = cbind(YI, NI-YI)~ 0+(1 |Rep),
                       family = binomial, data = D)
            mod2=glmer(formula = cbind(YI, NI-YI)~ 1+(1 |Rep),
                       family = binomial, data = D)
            mod3=glmer(formula = cbind(YI, NI-YI)~ 0+(1 |Rep)+(1|SNP),
                       family = binomial, data = D)
            mod4=glmer(formula = cbind(YI, NI-YI)~ 1+(1 |Rep)+(1|SNP),
                       family = binomial, data =D)
            A1<-BIC(mod1)
            A2<-BIC(mod2)
            A3<-BIC(mod3)
            A4<-BIC(mod4)
            if (A1==min(A1,A2,A3,A4)) {
                sigR[i]<-lapply(VarCorr(mod1),'[[',1)$Rep
            }
            else if (A2==min(A1,A2,A3,A4)) {
                sigR[i]<-lapply(VarCorr(mod2),'[[',1)$Rep
                beta[i]<- as.vector(fixef(mod2))  
            }
            else if (A3==min(A1,A2,A3,A4)) {
                sigS[i]<-lapply(VarCorr(mod3),'[[',1)$SNP
                sigR[i]<-lapply(VarCorr(mod3),'[[',1)$Rep
            }  
            else {
                sigS[i]<-lapply(VarCorr(mod4),'[[',1)$SNP
                sigR[i]<-lapply(VarCorr(mod4),'[[',1)$Rep
                beta[i]<- as.vector(fixef(mod4)) 
            }
        },error=function(e){cat("error: ",conditionMessage(e), "\n")})
    }
    # extract GLMM estiamtes
    sigRf<-sigR[!is.na(sigR)]
    sigSf<-sigS[!is.na(sigS)]
    betasf<-beta[!is.na(beta)]
    ## MOM to estimate paras of Inverse Gamma dist. for sigmaRsquare
    uRs<-mean(sigRf)
    vRs<-var(sigRf)
    aRs<-uRs^2/vRs+2
    bRs<-uRs*(uRs^2/vRs+1)
    ## MOM to estimate paras of Inverse Gamma dist. for sigmaSsquareS
    uSs<-mean(sigSf)
    vSs<-var(sigSf)
    aSs<-uSs^2/vSs+2
    bSs<-uSs*(uSs^2/vSs+1)
    ## MLE to estimate paras of Gaussian prior for fixed gene effect
    mlemus<- mean(betasf)
    mlesigmas<-sum((betasf-mlemus)^2)/length(betasf)
    calibration<-data.frame(aRs=aRs,bRs=bRs,aSs=aSs, bSs=bSs, 
                            mlemus=mlemus,mlesigmas=mlesigmas)
    return(calibration) 
} 


