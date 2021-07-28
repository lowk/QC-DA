library(sva) ### untidied but more completed code sources for reference: BatchCorrection.R; BatchTemp.R 

###functions to call
### k batch for KNN testing: input expression data and batch matrix, return rejection rate matrix
k=2
BatchEffectM = list()
KNNtest <- function(exprDat_norm,MetaRaw){
  ### k nearest neighborhood batch effect test
  for (roundCount in 1:10){
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
    
    BatchEffectM[[roundCount]] = sapply(batchTest,mean)
    # names(BatchEffectM) = c("PlateID","Disease Group","Corhort","BloodStain","SampleAge")
    # BatchEffectM <- as.table(t(as.matrix(BatchEffectM)))
    # rownames(BatchEffectM) = "Rejection Rate"
  }
  return(BatchEffectM)
}

### filter to exclude buffer|calibrator, select only clinical samples and human proteins
### input RFU matrix after normalisation steps, output filtered RFU matrix only with human data 
filterHM <- function(MySoma,BioMeta){
  HMpro <- which(!grepl("HybControlElution|NonBiotin|None",colnames(MySoma)))
  HMsam <- which(grepl("Sample",MySoma[,"SampleType"]))
  MySomaDone <- MySoma[HMsam,HMpro]
  BioMetaDone <- BioMeta[HMsam,]
  return(list(MySomaDone,BioMetaDone))
}

### get plate batch
GetPlateBatch <- function(MySomaPlate){
  plates <- levels(as.factor(MySomaPlate))
  
  PlateBatch = matrix(0,ncol=1,nrow=length(MySomaPlate))
  for(palteCounter in 1:length(MySomaPlate)){
    PlateBatch[palteCounter] = which(plates == MySomaPlate[palteCounter])
  }
  return(PlateBatch)
}


### Main body begin

### run QCdemoTranche2.R before line 59: calucate normalised tranche2 data

### load normalised tranche1 
load("MySoma1.RData")
load("MetaRawM1.RData")
BioMeta1 = MetaRawM1
BioMeta = rbind(BioMeta1[,c("diseaseGroup","Corhort")],BioMeta2[,c("diseaseGroup","Corhort")])
RawM1 <- read.adat(paste(myFilePath,"SS-200008.ADat",sep=""))

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
Sys.time()
combat_all_batch1 <- sva::ComBat(t(MySomaAll), tranche, mod=NULL, par.prior = TRUE, prior.plots = FALSE)
TrancheEffect1 = KNNtest(t(combat_all_batch1),MetaRaw) 

###2. source: log(combined(tranche1,2)); combat tranche only; parametric
Sys.time()
combat_all_batch2 <- sva::ComBat(t(log(MySomaAll)), tranche, mod=NULL, par.prior = TRUE, prior.plots = FALSE)
TrancheEffect2 = KNNtest(t(combat_all_batch2),MetaRaw) 

###3. source: combined(tranche1,2)); combat tranche only; nonparametric
Sys.time()
combat_all_batch3 <- sva::ComBat(t(MySomaAll), tranche, mod=NULL, par.prior = FALSE, prior.plots = FALSE)
TrancheEffect3 = KNNtest(t(combat_all_batch3),MetaRaw) 

###4. source: log(combined(tranche1,2)); combat tranche only; nonparametric
Sys.time()
combat_all_batch4 <- sva::ComBat(t(log(MySomaAll)), tranche, mod=NULL, par.prior = FALSE, prior.plots = FALSE)
TrancheEffect4 = KNNtest(t(combat_all_batch4),MetaRaw) 

###5. source: combined(tranche1,2); combat plate only; parametric
Sys.time()
combat_all_batch5 <- sva::ComBat(t(MySomaAll), PlateBatch, mod=NULL, par.prior = TRUE, prior.plots = FALSE)
TrancheEffect5 = KNNtest(t(combat_all_batch5),MetaRaw) 

###6. source: log(combined(tranche1,2)); combat plate only; parametric
Sys.time()
combat_all_batch6 <- sva::ComBat(t(log(MySomaAll)), PlateBatch, mod=NULL, par.prior = TRUE, prior.plots = FALSE)
TrancheEffect6 = KNNtest(t(combat_all_batch6),MetaRaw) 

###7. source: combined(tranche1,2)); combat plate only; nonparametric
Sys.time()
combat_all_batch7 <- sva::ComBat(t(MySomaAll), PlateBatch, mod=NULL, par.prior = FALSE, prior.plots = FALSE)
TrancheEffect7 = KNNtest(t(combat_all_batch7),MetaRaw) 

###8. source: log(combined(tranche1,2)); combat plate only; nonparametric
Sys.time()
combat_all_batch8 <- sva::ComBat(t(log(MySomaAll)), PlateBatch, mod=NULL, par.prior = FALSE, prior.plots = FALSE)
TrancheEffect8 = KNNtest(t(combat_all_batch8),MetaRaw) 
### combat_all_batch8 have constant protein feature, PCA need to exclude such vector
which(apply(combat_all_batch8,1,function(x) {length(unique(x))==1}))
combat_all_batch8 <- combat_all_batch8[-3255,]

###9. source: combined(tranche1,2); combat tranche + plate as covariate; parametric
Sys.time()
mod1=model.matrix(~as.factor(PlateBatch))
combat_all_batch9 <- sva::ComBat(t(MySomaAll), tranche, mod=mod1, par.prior = TRUE, prior.plots = FALSE)
TrancheEffect9 = KNNtest(t(combat_all_batch9),MetaRaw) 

###10. source: log(combined(tranche1,2)); combat tranche + plate as covariate; parametric
Sys.time()
combat_all_batch10 <- sva::ComBat(t(log(MySomaAll)), tranche, mod=mod1, par.prior = TRUE, prior.plots = FALSE)
TrancheEffect210 = KNNtest(t(combat_all_batch10),MetaRaw) 

###11. source: combined(tranche1,2)); combat tranche + plate as covariate; nonparametric
Sys.time()
combat_all_batch11 <- sva::ComBat(t(MySomaAll), tranche, mod=mod1, par.prior = FALSE, prior.plots = FALSE)
TrancheEffect11 = KNNtest(t(combat_all_batch11),MetaRaw) 

###12. source: log(combined(tranche1,2)); combat tranche + plate as covariate; nonparametric
Sys.time()
combat_all_batch12 <- sva::ComBat(t(log(MySomaAll)), tranche, mod=mod1, par.prior = FALSE, prior.plots = FALSE)
TrancheEffect12 = KNNtest(t(combat_all_batch12),MetaRaw) 

###13. source: combined(tranche1,2); combat tranche + plate as covariate; parametric
Sys.time()
mod2=model.matrix(~as.factor(tranche))
combat_all_batch13 <- sva::ComBat(t(MySomaAll), PlateBatch, mod=mod2, par.prior = TRUE, prior.plots = FALSE)
TrancheEffect13 = KNNtest(t(combat_all_batch13),MetaRaw) 

###14. source: log(combined(tranche1,2)); combat tranche + plate as covariate; parametric
Sys.time()
combat_all_batch14 <- sva::ComBat(t(log(MySomaAll)), PlateBatch, mod=mod2, par.prior = TRUE, prior.plots = FALSE)
TrancheEffect14 = KNNtest(t(combat_all_batch14),MetaRaw) 

###15. source: combined(tranche1,2)); combat tranche + plate as covariate; nonparametric
Sys.time()
combat_all_batch15 <- sva::ComBat(t(MySomaAll), PlateBatch, mod=mod2, par.prior = FALSE, prior.plots = FALSE)
TrancheEffect15 = KNNtest(t(combat_all_batch15),MetaRaw) 

###16. source: log(combined(tranche1,2)); combat tranche + plate as covariate; nonparametric
Sys.time()
combat_all_batch16 <- sva::ComBat(t(log(MySomaAll)), PlateBatch, mod=mod2, par.prior = FALSE, prior.plots = FALSE)
TrancheEffect16 = KNNtest(t(combat_all_batch16),MetaRaw) 

### Check Combat effectiveness: we recommend log(combined tranche1,2) + plate combat + parametric
###PCA display(add title manually) 
combat_all_PC<-combat_all_PC1
combat_all_PC <- prcomp(as.matrix(t(combat_all_batch)),scale = TRUE)
temp <- combat_all_PC$x[,1:3]
pairs(temp,col=tranche)
title("Tranche effect")
pairs(temp,col=PlateBatch)
title("Plate effect")
pairs(temp,col=as.factor(BioMeta[,"diseaseGroup"]))
title("Disease Group effect after Combat")
pairs(temp,col=as.factor(BioMeta[,"Corhort"]))
title("Cirhort effect after Combat")


###Include only the human samples and huma proteins, and perform the same analysis from main body line 132
MySoma1Done <- filterHM(MySoma1,BioMeta1)[[1]]
MySoma2Done <- filterHM(MySoma2,BioMeta2)[[1]]
BioMeta1Done <- filterHM(MySoma1,BioMeta1)[[2]]
BioMeta2Done <- filterHM(MySoma2,BioMeta2)[[2]]

MySomaAll <- rbind(MySoma1Done[,25:ncol(MySoma1Done)],MySoma2Done[,26:ncol(MySoma2Done)])
tranche <- c(rep(1,nrow(MySoma1Done)),rep(2,nrow(MySoma2Done)))
MySomaPlate <- as.matrix(c(MySoma1Done[,"PlateId"],MySoma2Done[,"PlateId"]),ncol=1)
PlateBatch <- GetPlateBatch(MySomaPlate)
MetaRaw = cbind(tranche,PlateBatch)
rownames(MetaRaw)=rownames(MySomaAll)
colnames(MetaRaw) = c("Tranche Batch","Plate Batch")
PlateBatch <- as.vector(PlateBatch) ### Combat argument requirement
BioMeta <- rbind(BioMeta1Done[,2:3],BioMeta2Done[,1:2])

###BioMeta to match MySomaAll rownames
### Before Combat
MySoma_PC <- prcomp(as.matrix(log(MySomaAll)),scale = TRUE)
temptemp <- MySoma_PC$x[,1:3]
PlotDat <- data.frame(cbind(temptemp,MetaRaw))

ggpairs(PlotDat, title="Tranche Effect before Combat",
        columns=1:3, aes(color= as.factor(Tranche.Batch)),
        diag=list(continuous=wrap("densityDiag",alpha=0.4)),
        lower=list(continuous = wrap("points",alpha=0.9,size=0.1)),
        upper = list(continuous = "blank"),
        legend = c(1,1)) + labs(fill = colnames(MetaRaw)[1]) + theme(plot.title=element_text(size=15,hjust=0.5))

ggpairs(PlotDat, title="Plate Effect before Combat",
        columns=1:3, aes(color= as.factor(Plate.Batch)),
        diag=list(continuous=wrap("densityDiag",alpha=0.4)),
        lower=list(continuous = wrap("points",alpha=0.9,size=0.1)),
        upper = list(continuous = "blank"),
        legend = c(1,1)) + labs(fill = colnames(MetaRaw)[2]) + theme(plot.title=element_text(size=15,hjust=0.5))

ggpairs(PlotDat, title="Disease Effect before Combat",
        columns=1:3, aes(color= as.factor(BioMeta[,"diseaseGroup"])),
        diag=list(continuous=wrap("densityDiag",alpha=0.4)),
        lower=list(continuous = wrap("points",alpha=0.9,size=0.1)),
        upper = list(continuous = "blank"),
        legend = c(1,1)) + labs(fill = colnames(BioMeta)[1]) + theme(plot.title=element_text(size=15,hjust=0.5))

ggpairs(PlotDat, title="Cohort Effect before Combat",
        columns=1:3, aes(color= as.factor(BioMeta[,"Corhort"])),
        diag=list(continuous=wrap("densityDiag",alpha=0.4)),
        lower=list(continuous = wrap("points",alpha=0.9,size=0.1)),
        upper = list(continuous = "blank"),
        legend = c(1,1)) + labs(fill = colnames(BioMeta)[2]) + theme(plot.title=element_text(size=15,hjust=0.5))

###after combat
combat_all_PC <- prcomp(as.matrix(t(combat_all_batch6)),scale = TRUE)
temptemptemp <- combat_all_PC$x[,1:3]
PlotDat2 <- data.frame(cbind(temptemptemp,MetaRaw))

ggpairs(PlotDat2, title="Tranche Effect after Combat, approach6",
        columns=1:3, aes(color= as.factor(Tranche.Batch)),
        diag=list(continuous=wrap("densityDiag",alpha=0.4)),
        lower=list(continuous = wrap("points",alpha=0.9,size=0.1)),
        upper = list(continuous = "blank"),
        legend = c(1,1)) + labs(fill = colnames(MetaRaw)[1]) + theme(plot.title=element_text(size=15,hjust=0.5))

ggpairs(PlotDat2, title="Plate Effect after Combat, approach6",
        columns=1:3, aes(color= as.factor(Plate.Batch)),
        diag=list(continuous=wrap("densityDiag",alpha=0.4)),
        lower=list(continuous = wrap("points",alpha=0.9,size=0.1)),
        upper = list(continuous = "blank"),
        legend = c(1,1)) + labs(fill = colnames(MetaRaw)[2]) + theme(plot.title=element_text(size=15,hjust=0.5))

ggpairs(PlotDat2, title="Disease Effect after Combat, selected approach",
        columns=1:3, aes(color= as.factor(BioMeta[,"diseaseGroup"])),
        diag=list(continuous=wrap("densityDiag",alpha=0.4)),
        lower=list(continuous = wrap("points",alpha=0.9,size=0.1)),
        upper = list(continuous = "blank"),
        legend = c(1,1)) + labs(fill = colnames(BioMeta)[1]) + theme(plot.title=element_text(size=15,hjust=0.5))

ggpairs(PlotDat2, title="Cohort Effect after Combat, selected approach",
        columns=1:3, aes(color= as.factor(BioMeta[,"Corhort"])),
        diag=list(continuous=wrap("densityDiag",alpha=0.4)),
        lower=list(continuous = wrap("points",alpha=0.9,size=0.1)),
        upper = list(continuous = "blank"),
        legend = c(1,1)) + labs(fill = colnames(BioMeta)[2]) + theme(plot.title=element_text(size=15,hjust=0.5))


pairs(temptemptemp,col=tranche)
title("Tranche effect after combat",line=3)
pairs(temptemp,col=PlateBatch)
title("Plate effect")
pairs(temptemptemp,col=as.factor(BioMeta[,"diseaseGroup"]),line=1)
title("Disease Group effect after norm but before Combat, selected approach")
pairs(temptemptemp,col=as.factor(BioMeta[,"Corhort"]))
title("Cirhort effect after norm but before Combat, selected approch")

RawM1 <- read.adat(inputfile1)
RawAll <- rbind(RawM1[,25:ncol(RawM1)],RawM2[,26:ncol(RawM2)])
RawAll_PC <- prcomp(as.matrix(RawAll),scale = TRUE)
temptemptemp <- RawAll_PC$x[,1:3]
pairs(temptemptemp,col=tranche)
title("Tranche effect Raw")
pairs(temptemptemp,col=PlateBatch)
title("Plate effect Raw")
pairs(temptemptemp,col=as.factor(BioMeta[,"diseaseGroup"]))
title("Disease Group effect Raw")
pairs(temptemptemp,col=as.factor(BioMeta[,"Corhort"]))
title("Corhort effect Raw")

###Tranche2: levels(as.factor): plate order is this: "P0028435" "P0028436" "P0028443" "P0028444" "P0028445" "P0028446" "P0028453"
###But in the adat file, plate order is: "P0028435" "P0028436" "P0028443" "P0028444" "P0028453" "P0028445" "P0028446"
### always check the row order of meta data and raw and MySoma

###KNN test
batcheMean=matrix(0,ncol=2,nrow=8)
batchSd=matrix(0,ncol=2,nrow=8)
for (i in c(1,2,5,6,7,8)){
  BatchEffectM=list()
  
  testM <- paste("combat_all_batch",i,sep="")
  BatchEffectM = KNNtest(t(get(testM)),MetaRaw) 
  BatchEffectM2 = matrix(unlist(BatchEffectM),nrow=10,ncol=2,byrow=TRUE)
  batcheMean[i,] <- apply(BatchEffectM2,2,mean)
  batchSd[i,] <- apply(BatchEffectM2,2,sd)
}


