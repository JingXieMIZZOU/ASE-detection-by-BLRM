###############################################################

#detection(data) returns a list contains  genes which has gene
# effect, SNP effect, and gene & SNP effect.

##############################################################



detection<-function(data){
    
    # note function para.est() accounts for data augmentation process but PPSFUN doesn't.
    # so, need to do data augmentation before calculate posterior probability 
    
    # hyper parameter estimation
    clean_index<-filter(data)
    paras<-para.est(data,index=clean_index)
    # data augmentation before calculate posterior probability
    data[,c(4,6,8,10)]<-data[,c(4,6,8,10)]+2
    data[,c(5,7,9,11)]<-data[,c(5,7,9,11)]+1
    ## realPPS is the log of PP before standardize
    temp<-lapply(clean_index,PPsFUN,par=paras,data=data)
    realPP<-do.call(rbind, lapply(temp, head, 1))
    ## normalized PPs
    realPPS<-exp(realPP[,2:5])
    realPPSS<-matrix(nrow=nrow(realPP),ncol=6)
    realPPSS<-as.data.frame(realPPSS)
    realPPSS[,1]<-realPP[,1]
    realPPSS[,2]<-unlist(lapply(as.vector(realPP[,1]),function(x) unique(data[data$gene==x,]$ensemble_id)))  
    realPPSS[,3]<-realPPS[,1]/(apply(realPPS,1,sum))
    realPPSS[,4]<-realPPS[,2]/(apply(realPPS,1,sum))
    realPPSS[,5]<-realPPS[,3]/(apply(realPPS,1,sum))
    realPPSS[,6]<-realPPS[,4]/(apply(realPPS,1,sum))
    colnames(realPPSS)<-c("GeneNum","GeneID","PP1","PP2","PP3","PP4")
    
    # re-identification
    # gene effect
    real_gene<-data.frame(geneNum=realPPSS$GeneNum,
                          geneID=realPPSS$GeneID,
                          PP=realPPSS$PP2+realPPSS$PP4,
                          FPP=1-(realPPSS$PP2+realPPSS$PP4))
    real_gene<-real_gene[order(real_gene$PP,decreasing = TRUE),]
    real_gene$FDR<-cummean(real_gene$FPP)
    real_gene_res<-real_gene[real_gene$FDR<=0.05,]
    
    # SNP effect
    # re-identification
    real_SNP<-data.frame(geneNum=realPPSS$GeneNum,
                         geneID=realPPSS$GeneID,
                         PP=realPPSS$PP3+realPPSS$PP4,
                         FPP=1-(realPPSS$PP3+realPPSS$PP4))
    real_SNP<-real_SNP[order(real_SNP$PP,decreasing = TRUE),]
    real_SNP$FDR<-cummean(real_SNP$FPP)
    real_SNP_res<-real_SNP[real_SNP$FDR<=0.05,]
    
    # Gene & SNP effect
    real_GS<-data.frame(geneNum=realPPSS$GeneNum,
                        geneID=realPPSS$GeneID,
                        PP=realPPSS$PP4,
                        FPP=1-(realPPSS$PP4))
    real_GS<-real_GS[order(real_GS$PP,decreasing=TRUE),]
    real_GS$FDR<-cummean(real_GS$FPP)
    real_GS_res<-real_GS[real_GS$FDR<=0.05,]

    return(list(GeneEffect=na.omit(real_gene_res),
                SNPEffect=na.omit(real_SNP_res),
                GSEffect=na.omit(real_GS_res)))
}


