library(sva) ### untidied but more completed code sources for reference: BatchCorrection.R; BatchTemp.R 
###functions to call

### k batch for KNN testing: input expression data and batch matrix, return rejection rate matrix
k=2
KNNtest <- function(exprDat_norm,MetaRaw){
  ### k nearest neighborhood batch effect test
  for (roundCount in 1:1){
    pc_norm <- prcomp(exprDat_norm,scale = TRUE)
    topPCn <- which(get_eigenvalue(pc_norm)$cumulative.variance.percent>80)[1]
    distPairs = pc_norm$x[,1:topPCn]
    DisM = as.matrix(dist(distPairs, method = "euclidean", diag = TRUE, upper = TRUE))
    
    sampSizeS = seq(10,500,by=100)  ### different percentage of 1%~25% samples for batch effect test (based on the paper)
    kNearS = seq(10,500,by=100)          ### different neiboughood sizes, based on testing, >250 no rejection at all
    positiveRate = matrix(NA,nrow=length(kNearS),ncol=length(sampSizeS))
    batchTest = vector(mode="list",length=k)
    
    
    for (batchCf in 1:ncol(MetaRaw)){ 
      tblWh = table(MetaRaw[,batchCf])
      EpectedRio = tblWh/sum(tblWh)
      
      for (sampSct in 1:length(sampSizeS)){
        sampSize = sampSizeS[sampSct]
        randSampId = sample(1:nrow(exprDat_norm),sampSize) ### random select which samples to be tested
        
        for (neiSize in 1:length(kNearS)){
          
          kNear = kNearS[neiSize]   ###set neighborhood 
          
          idInMetaR = matrix(NA,nrow=kNear-1,ncol=1)
          testSampYN = matrix(NA,nrow=sampSize,ncol=1)
          
          
          for (testP in  1:sampSize){
            
            testSamp = DisM[randSampId[testP],]  ###define test sample
            neighborD = sort(testSamp)[2:kNear]  ### find its nearest kNear neighbors
            neiborOriID = names(neighborD)       ### find its neighbor name, for batch label matching
            
            neiNameList = rownames(MetaRaw)
            
            for (neiCounter in 1:(kNear-1)){ 
              idInMetaR[neiCounter] = which(neiNameList %in% neiborOriID[neiCounter]) ### meta for tested sample
            }
            
            neiDat = MetaRaw[idInMetaR,batchCf]
            tblPart = table(neiDat)
            
            missingCtg = apply(as.matrix(rownames(as.matrix(tblPart))),2,function(x){which(!(rownames(as.matrix(tblWh)) %in% x))})
            if (length(missingCtg) != 0) {### for categories with 0 count
              missingVl = matrix(0,ncol=length(missingCtg))
              colnames(missingVl) <- rownames(as.matrix(tblWh))[missingCtg]
              tblPart = cbind(t(as.matrix(tblPart)),missingVl)
              tblTest = rbind(tblPart[,rownames(as.matrix(tblWh))],tblWh) ### make col orders consistent
            }
            
            else{tblTest = rbind(tblPart,tblWh)}
            
            ### chi square test or exact test, check expected value for each element >5?
            if(min(rowSums(tblTest) %*% t(colSums(tblTest)) / sum(tblTest)) >5){
              reportYN = chisq.test(tblTest)$p.value}
            else{reportYN = fisher.test(tblTest,simulate.p.value = TRUE)$p.value}
            
            ### chi square based multinomial test, power is less than above one
            # chi.sq.value <- sum((tblTest[1,]/sampSize - EpectedRio)^2/EpectedRio)
            #  
            # reportYN <- 1 - pchisq(chi.sq.value, df=length(unique(MetaRaw[,batchCf]))-1) #p-value for the result
            
            ### hypothesis test, local batch label distribution against the global label distribution
            if(reportYN < 0.05/sampSize) {testSampYN [testP] = 1} ###using Bonferroni adjusted p values
            else {testSampYN[testP] = 0}
            
          } 
          
          positiveRate[neiSize,sampSct] = length(which(testSampYN == 1)) / length(testSampYN)
          colnames(positiveRate) = sampSizeS
          rownames(positiveRate) = kNearS
          batchTest[[batchCf]]=positiveRate ### columns: samples for test, rows: neighborhood
          
        }
        
      }
      
    }
    
    # boxplot(batchTest[[1]],xlab = "tested sample amount",ylab="Percentage of rejection H0 (Neighborhood size 25 ~ 250)",main = "One Round of Plate Batch Effect Test on Raw RFUs",cex.lab=0.7,cex.main=0.8,cex.axis=0.8)
    
    BatchEffectM = sapply(batchTest,mean)
    # names(BatchEffectM) = c("PlateID","Disease Group","Corhort","BloodStain","SampleAge")
    BatchEffectM <- as.table(t(as.matrix(BatchEffectM)))
    rownames(BatchEffectM) = "Rejection Rate"
  }
  return(BatchEffectM)
}

### filter to exclude buffer|calibrator, select only clinical samples and human proteins
### input RFU matrix after normalisation steps, output filtered RFU matrix only with human data 
filterHM <- function(MySoma){
  HMpro <- which(!grepl("HybControlElution|NonBiotin|None",colnames(MySoma)))
  HMsam <- which(grepl("Sample",MySoma[,"SampleType"]))
  MySomaDone <- MySoma[HMsam,HMpro]
  return(MySomaDone)
}

### get plate batch per tranche
GetPlateBatch <- function(MySomaPlate){
  plates <- levels(as.factor(MySomaPlate))
  
  PlateBatch = matrix(0,ncol=1,nrow=length(MySomaPlate))
  for(palteCounter in 1:length(MySomaPlate)){
    PlateBatch[palteCounter] = which(plates == MySomaPlate[palteCounter])
  }
  return(PlateBatch)
}


### Main body begin
### combine tranches after normalisation separately
MySomaAll <- rbind(MySoma1[,25:ncol(MySoma1)],MySoma2[,26:ncol(MySoma2)])

### construct the plate batch and tranche batch
tranche <- c(rep(1,nrow(MySoma1)),rep(2,nrow(MySoma2)))
MySomaPlate <- as.matrix(c(MySoma1[,"PlateId"],MySoma2[,"PlateId"]),ncol=1)
PlateBatch <- GetPlateBatch(MySomaPlate)
MetaRaw = cbind(tranche,PlateBatch)
rownames(MetaRaw)=rownames(MySomaAll)
colnames(MetaRaw) = c("Tranche Batch","Plate Batch")
PlateBatch <- as.vector(PlateBatch) ### Combat argument requirement

### Combat normalized RFU matrix
###1. source: combined(tranche1,2); combat tranche only; parametric
combat_all_batch <- sva::ComBat(t(MySomaAll), tranche, mod=NULL, par.prior = TRUE, prior.plots = FALSE)

###2. source: log(combined(tranche1,2)); combat tranche only; parametric
combat_all_batch <- sva::ComBat(t(log(MySomaAll)), tranche, mod=NULL, par.prior = TRUE, prior.plots = FALSE)

###3. source: combined(tranche1,2)); combat tranche only; nonparametric
combat_all_batch <- sva::ComBat(t(MySomaAll), tranche, mod=NULL, par.prior = FALSE, prior.plots = FALSE)

###4. source: log(combined(tranche1,2)); combat tranche only; nonparametric
combat_all_batch <- sva::ComBat(t(log(MySomaAll)), tranche, mod=NULL, par.prior = FALSE, prior.plots = FALSE)

###5. source: combined(tranche1,2); combat plate only; parametric
combat_all_batch <- sva::ComBat(t(MySomaAll), PlateBatch, mod=NULL, par.prior = TRUE, prior.plots = FALSE)

###6. source: log(combined(tranche1,2)); combat plate only; parametric
combat_all_batch <- sva::ComBat(t(log(MySomaAll)), PlateBatch, mod=NULL, par.prior = TRUE, prior.plots = FALSE)

###7. source: combined(tranche1,2)); combat plate only; nonparametric
combat_all_batch <- sva::ComBat(t(MySomaAll), PlateBatch, mod=NULL, par.prior = FALSE, prior.plots = FALSE)

###8. source: log(combined(tranche1,2)); combat plate only; nonparametric
combat_all_batch <- sva::ComBat(t(log(MySomaAll)), PlateBatch, mod=NULL, par.prior = FALSE, prior.plots = FALSE)

###9. source: combined(tranche1,2); combat tranche + plate as covariate; parametric
mod1=model.matrix(~as.factor(PlateBatch))
combat_all_batch <- sva::ComBat(t(MySomaAll), tranche, mod=mod1, par.prior = TRUE, prior.plots = FALSE)

###10. source: log(combined(tranche1,2)); combat tranche + plate as covariate; parametric
combat_all_batch <- sva::ComBat(t(log(MySomaAll)), tranche, mod=mod1, par.prior = TRUE, prior.plots = FALSE)

###11. source: combined(tranche1,2)); combat tranche + plate as covariate; nonparametric
combat_all_batch <- sva::ComBat(t(MySomaAll), tranche, mod=mod1, par.prior = FALSE, prior.plots = FALSE)

###12. source: log(combined(tranche1,2)); combat tranche + plate as covariate; nonparametric
combat_all_batch <- sva::ComBat(t(log(MySomaAll)), tranche, mod=mod1, par.prior = FALSE, prior.plots = FALSE)

###13. source: combined(tranche1,2); combat tranche + plate as covariate; parametric
mod2=model.matrix(~as.factor(tranche))
combat_all_batch <- sva::ComBat(t(MySomaAll), PlateBatch, mod=mod2, par.prior = TRUE, prior.plots = FALSE)

###14. source: log(combined(tranche1,2)); combat tranche + plate as covariate; parametric
combat_all_batch <- sva::ComBat(t(log(MySomaAll)), PlateBatch, mod=mod2, par.prior = TRUE, prior.plots = FALSE)

###15. source: combined(tranche1,2)); combat tranche + plate as covariate; nonparametric
combat_all_batch <- sva::ComBat(t(MySomaAll), PlateBatch, mod=mod2, par.prior = FALSE, prior.plots = FALSE)

###16. source: log(combined(tranche1,2)); combat tranche + plate as covariate; nonparametric
combat_all_batch <- sva::ComBat(t(log(MySomaAll)), PlateBatch, mod=mod2, par.prior = FALSE, prior.plots = FALSE)

### Check Combat effectiveness: we recommend log(combined tranche1,2) + plate combat + parametric
###PCA display(add title manually) 
combat_all_PC <- prcomp(as.matrix(t(combat_all_batch)),scale = TRUE)
temp <- combat_all_PC$x[,1:3]
par(mfrow=c(1,2))
pairs(temp,col=tranche)
title("Tranche effect")
pairs(temp5,col=PlateBatch)
title("Plate effect")
### KNN test
TrancheEffect = KNNtest(t(combat_all_batch),MetaRaw) 


###Include only the human samples and huma proteins, perform the same analysis from main body
MySoma1Done <- filterHM(MySoma1)
MySoma2Done <- filterHM(MySoma2)
MySomaAll <- rbind(MySoma1Done[,25:ncol(MySoma1Done)],MySoma1Done[,26:ncol(MySoma2Done)])
