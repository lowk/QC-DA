library(sva)

### k batch for testing
k=2
KNNtest <- function(exprDat_norm,MetaRaw){
  ### k nearest neighborhood batch effect test
  for (roundCount in 1:1){
    pc_norm <- prcomp(exprDat_norm,scale = TRUE)
    distPairs = pc_norm$x[,1:10]
    DisM = as.matrix(dist(distPairs, method = "euclidean", diag = TRUE, upper = TRUE))
    
    sampSizeS = seq(10,500,by=100)  ### different percentage of 1%~25% samples for batch effect test (based on the paper)
    kNearS = seq(10,500,by=100)          ### different neiboughood sizes, based on testing, >250 no rejection at all
    positiveRate = matrix(NA,nrow=length(kNearS),ncol=length(sampSizeS))
    batchTest = vector(mode="list",length=k)
    
    
    for (batchCf in 1:ncol(MetaRaw)){ ###batchCf: "PlateID" "diseaseGroup" "Corhort" "bloodStain" "sampleAge"   
      
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
            
            
            ### chi square based multinomial test
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


### check the row names and col names of these two tranche data
### combine the two tranche data: MySomaDatAll
all(colnames(MySoma1[,25:5308])==colnames(MySoma2[,26:5309]))
MySomaDatAll <- rbind(MySoma1[,25:5308],MySoma2[,26:5309])

load("MySoma1.RData")

### construct the plate batch and tranche batch
batchV <- c(rep(1,nrow(MySoma1)),rep(2,nrow(MySoma2)))

PlateBatch1 <- matrix(0,ncol=1)
for (i in 1:nrow(MySoma1)){
  if(MySoma1[i,1]=="P0025226"){PlateBatch1[i]=1}
  if(MySoma1[i,1]=="P0025228"){PlateBatch1[i]=2}
  if(MySoma1[i,1]=="P0025230"){PlateBatch1[i]=3}
  if(MySoma1[i,1]=="P0025232"){PlateBatch1[i]=4}
  if(MySoma1[i,1]=="P0025234"){PlateBatch1[i]=5}
}

PlateBatch2 <- matrix(0,ncol=1)
for (i in 1:nrow(MySoma2)){
  if(MySoma2[i,1]=="P0028435"){PlateBatch2[i]=6}
  if(MySoma2[i,1]=="P0028436"){PlateBatch2[i]=7}
  if(MySoma2[i,1]=="P0028443"){PlateBatch2[i]=8}
  if(MySoma2[i,1]=="P0028444"){PlateBatch2[i]=9}
  if(MySoma2[i,1]=="P0028445"){PlateBatch2[i]=10}
  if(MySoma2[i,1]=="P0028446"){PlateBatch2[i]=11}
  if(MySoma2[i,1]=="P0028453"){PlateBatch2[i]=12}
}

PlateBatch <- c(PlateBatch1,PlateBatch2)

MetaRaw = cbind(batchV,PlateBatch)
rownames(MetaRaw)=rownames(MySomaDatAll)


### PCA check batch effect before combat
pc_norm <- prcomp(log10(as.matrix(MySomaDatAll)),scale = TRUE)
temp = pc_norm$x[,1:3]
names(temp) <- paste("PC",1:3)
pairs(temp,col=batchV)
title("Trache effect before combat")
pairs(temp,col=PlateBatch)
title("Plate effect before combat")
TrancheEffect1 = KNNtest(MySomaDatAll,MetaRaw) 

combat_all_batch <- sva::ComBat(t(MySomaDatAll), batchV, mod=NULL, par.prior = TRUE, prior.plots = FALSE)
combat_all_combat_all_batch_pc <- prcomp(as.matrix(t(combat_all_batch)),scale = TRUE)
tempAll = combat_all_combat_all_batch_pc$x[,1:3]
names(tempAll) <- paste("PC",1:3)
pairs(tempAll,col=batchV)
title("Trache effect after combat only on tranche batch")
pairs(tempAll,col=PlateBatch)
title("Plate effect after combat only on tranche batch")
TrancheEffect2 = KNNtest(t(combat_all_batch),MetaRaw) 

combat_all_plate_noTrache <- sva::ComBat(t(MySomaDatAll), PlateBatch, mod=NULL, par.prior = TRUE, prior.plots = FALSE)
combat_all_plate_noTrache_PC <- prcomp(as.matrix(t(combat_all_plate_noTrache)),scale = TRUE)
temp4 <- combat_all_plate_noTrache_PC$x[,1:3]
pairs(temp4,col=batchV)
title("Tranche effect after combat on Plate only")
pairs(temp4,col=PlateBatch)
title("Plate effect after combat on Plate only")
TrancheEffect3 = KNNtest(t(combat_all_plate_noTrache),MetaRaw) 

combat_all_plate2 <- sva::ComBat(combat_all_batch, PlateBatch, mod=NULL, par.prior = TRUE, prior.plots = FALSE)
combat_all_plate2_PC <- prcomp(as.matrix(t(combat_all_plate2)),scale = TRUE)
temp3 <- combat_all_plate2_PC$x[,1:3]
pairs(temp3,col=batchV)
title("Tranche effect after combat on Trache|Plate sequentially")
pairs(temp3,col=PlateBatch)
title("Plate effect after combat on Trache|Plate sequentially")
TrancheEffect4 = KNNtest(t(combat_all_plate2),MetaRaw) 

combat_all_plate_tranche <- sva::ComBat(t(MySomaDatAll), batchV, mod=PlateBatch, par.prior = TRUE, prior.plots = FALSE)
combat_all_plate_tranche_PC <- prcomp(as.matrix(t(combat_all_plate_tranche)),scale = TRUE)
temp5 <- combat_all_plate_tranche_PC$x[,1:3]
pairs(temp5,col=batchV)
title("Tranche effect after combat on Trache|Plate together")
pairs(temp5,col=PlateBatch)
title("Plate effect after combat on Trache|Plate together")
TrancheEffect5 = KNNtest(t(combat_all_plate_tranche),MetaRaw) 

combat_all_plate_tranche2 <- sva::ComBat(t(MySomaDatAll), PlateBatch, mod=batchV, par.prior = TRUE, prior.plots = FALSE)
combat_all_plate_tranche_PC2 <- prcomp(as.matrix(t(combat_all_plate_tranche2)),scale = TRUE)
temp5 <- combat_all_plate_tranche_PC2$x[,1:3]
pairs(temp5,col=batchV)
title("Tranche effect after combat on Trache|Plate together")
pairs(temp5,col=PlateBatch)
title("Plate effect after combat on Trache|Plate together")

### plate effect test based on normalised 2 tranches
Tranche1plate <-as.matrix(MySoma1[,1],ncol=1)
rownames(Tranche1plate) <- rownames(MySoma1)
Tranche2plate <-as.matrix(MySoma2[,1],ncol=1)
rownames(Tranche2plate) <- rownames(MySoma2)

TrancheEffect6 = KNNtest(MySoma1[,25:ncol(MySoma1)],Tranche1plate) 
TrancheEffect7 = KNNtest(MySoma2[,26:ncol(MySoma2)],Tranche2plate) 

### Combat plate on tranche 2 
express2 <- MySoma2[,26:ncol(MySoma2)]
combat_plate_tranche2 <- sva::ComBat(t(express2), PlateBatch2, mod=NULL, par.prior = TRUE, prior.plots = FALSE)
NewExp <- rbind(MySoma1[,25:ncol(MySoma1)],t(combat_plate_tranche2))
combat_final <- sva::ComBat(t(NewExp), batchV, mod=NULL, par.prior = TRUE, prior.plots = FALSE)
TrancheEffect8 = KNNtest(t(combat_final),MetaRaw) 
combat_final2 <- sva::ComBat(t(NewExp), batchV, mod=PlateBatch, par.prior = TRUE, prior.plots = FALSE)
TrancheEffect9 = KNNtest(t(combat_final2),MetaRaw) 


### check the plate effect of tranche 2 Raw data
Tranche2RawPlate <-as.matrix(RawM[,1],ncol=1)
rownames(Tranche2RawPlate) <- rownames(RawM)
TrancheEffect10 = KNNtest(RawM[,26:ncol(RawM)],Tranche2RawPlate) 
pc_tranche2_raw <- prcomp(RawM[,26:ncol(RawM)],scale = TRUE)
pc_tranche2_raw_top <- pc_tranche2_raw$x[,1:3]
pairs(pc_tranche2_raw_top,col=as.factor(Tranche2RawPlate))
title("tranche2 plate effect before Combat")
pc_tranche2_plateCombat <- prcomp(t(combat_plate_tranche2),scale = TRUE)
pc_tranche2_plateCombat_top <- pc_tranche2_plateCombat$x[,1:3]
pairs(pc_tranche2_plateCombat_top,col=as.factor(Tranche2RawPlate))
title("tranche2 plate effect after Combat")
