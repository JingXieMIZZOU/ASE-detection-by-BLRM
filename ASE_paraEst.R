###############################################################################

# GDD(i,data) is a function to re-structurize the data to a structure that    # 
# contains necessary information for analysis, i.e., SNPs, Replicates, counts #
# from maternal allel (YI) and total counts (NI). This structure is friendly  #
# to GLMM fitting in glmer() function in "lme4" package. This function can be #                      
# easily modified to apply to other similar real data structures.             #

#i is the gene number of each gene. Take the sample dataset for instance, the #
# gene number of the first gene in this dataset is 28.                       #

# > unique(mysample$gene)[1]
# [1] 28
#  > GDD(28,data=mysample)
#   SNP NI YI Rep
#    1  6  3   1
#    1  7  2   2
#    1  3  1   3
#    2  6  5   2
#    3 14  4   1
#    3  4  1   2
#    3 18  8   3
#    3 13  7   4
#    4  8  2   1
#    4  7  3   2
#    4  5  4   3
#    4 10  3   4
#    5  8  2   1
#    5  5  2   2
#    5 10  6   3
#    5 13  7   4
#    6  8  2   1
#    6  9  3   2
#    6 10  4   3
#    6 12  5   4
#    7  5  2   1
#    7  4  1   2
#    7  9  0   3
#    7  9  2   4

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
#filter(data) is a function to filter the genes with computational problems
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
    # index is the result returned by above function filter(data)
    ###### add two obs and then fit GLMM 
    geneNum<-index
    P_SNP<- P_gene<- B4<- R4<- S4<- rep(NA,length(geneNum))
    for (i in geneNum){
        DATA<-GDD(i,data)
        D<-data.frame(SNP=DATA$SNP,NI=DATA$NI+2,YI=DATA$YI+1,Rep=DATA$Rep)
        tryCatch({
            print(i)
            mod2=glmer(formula = cbind(YI, NI-YI)~ 1+(1 |Rep),
                       family = binomial, data = D)
            mod3=glmer(formula = cbind(YI, NI-YI)~ 0+(1 |Rep)+(1|SNP),
                       family = binomial, data = D)
            mod4=glmer(formula = cbind(YI, NI-YI)~ 1+(1 |Rep)+(1|SNP),
                       family = binomial, data =D)
            B4[i]<- as.vector(fixef(mod4))
            R4[i]<-lapply(VarCorr(mod4),'[[',1)$Rep
            S4[i]<-lapply(VarCorr(mod4),'[[',1)$SNP
            P_SNP[i]<-anova(mod2,mod4)$`Pr(>Chisq)`[2]
            P_gene[i]<- anova(mod3,mod4)$`Pr(>Chisq)`[2]
            }
        ,error=function(e){cat("error: ",conditionMessage(e), "\n")})
    }
    # organize resutls
    res<- data.frame(geneNum=geneNum, P_gene=P_gene, P_SNP=P_SNP, beta=B4, Rep=R4, SNP=S4)
    res$q_gene<-qvalue(res$P_gene)$qvalues
    res$q_SNP<-qvalue(res$P_SNP)$qvalues
    set.seed(23457)
    tiebreaker <- sample(1:nrow(res), replace=TRUE)
    ## MOM to estimate paras of Inverse Gamma dist. for sigmaRsquare
    sigRf<-na.omit(res$Rep)
    uRs<-mean(sigRf)
    vRs<-var(sigRf)
    aRs<-uRs^2/vRs+2
    bRs<-uRs*(uRs^2/vRs+1)
    ## MOM to estimate paras of Inverse Gamma dist. for sigmaSsquareS
    res_SNP<-res[order(res$P_SNP,tiebreaker,decreasing = F),]
    sigSf<-na.omit(res_SNP[res_SNP$q_SNP<=0.05,]$SNP)
    uSs<-mean(sigSf)
    vSs<-var(sigSf)
    aSs<-uSs^2/vSs+2
    bSs<-uSs*(uSs^2/vSs+1)
    ## MLE to estimate paras of Gaussian prior for fixed gene effect
    res_gene<-res[order(res$P_gene,tiebreaker,decreasing = F),]
    betas<-res_gene[res_gene$q_gene<=0.05,]$beta
    betasf<-betas[!is.na(betas)]
    mlemus<- mean(betasf)
    mlesigmas<-sqrt(sum((betasf-mlemus)^2)/length(betasf))
    calibration<-data.frame(aRs=aRs,bRs=bRs,aSs=aSs, bSs=bSs, 
                            mlemus=mlemus,mlesigmas=mlesigmas)
    return(calibration) 
} 


