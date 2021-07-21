### too many duplicated codes, simplify as soon as possible

Tranche1plate <-as.matrix(MySoma1[,1],ncol=1)
rownames(Tranche1plate) <- rownames(MySoma1)
Tranche2plate <-as.matrix(MySoma2[,1],ncol=1)
rownames(Tranche2plate) <- rownames(MySoma2)

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
MySomaDatAll <- rbind(MySoma1[,25:5308],MySoma2[,26:5309])
rownames(MetaRaw)=rownames(MySomaDatAll)

### Normalize and Combat plate on tranche 1 

### Tranche1 Raw  
Tranche1RawPlate <-as.matrix(RawM1[,1],ncol=1)
rownames(Tranche1RawPlate) <- rownames(RawM1)
TrancheEffect11 = KNNtest(RawM1[,25:ncol(RawM1)],Tranche1RawPlate) 

### Tranche1 Normalisation
TrancheEffect6 = KNNtest(MySoma1[,25:ncol(MySoma1)],Tranche1plate) 

### Tranche1 Normalisation + plate Combat
express1 <- MySoma1[,25:ncol(MySoma1)]
combat_plate_tranche1 <- sva::ComBat(t(express1), PlateBatch1, mod=NULL, par.prior = TRUE, prior.plots = FALSE)
PlateBatchT1 <- as.matrix(MySoma1[,1],ncol=1)  ### Combat and KNNtest function required different class type of PlateBatch
rownames(PlateBatchT1) <- rownames(MySoma1)
TrancheEffect12 = KNNtest(t(combat_plate_tranche1),PlateBatchT1) 

### Tranche1 RawM + plate combat
expressRaw1 <- RawM1[,25:ncol(RawM1)]
combat_plate_trancheRaw1 <- sva::ComBat(t(expressRaw1), PlateBatch1, mod=NULL, par.prior = TRUE, prior.plots = FALSE)
PlateBatchRawT1 <- as.matrix(RawM1[,1],ncol=1)  ### Combat and KNNtest function required different class type of PlateBatch
rownames(PlateBatchRawT1) <- rownames(RawM1)
TrancheEffect13 = KNNtest(t(combat_plate_trancheRaw1),PlateBatchT1) 


### Normalize and Combat plate on tranche 2 

### Tranche2 Raw  
Tranche2RawPlate <-as.matrix(RawM2[,1],ncol=1)
rownames(Tranche2RawPlate) <- rownames(RawM2)
TrancheEffect22 = KNNtest(RawM2[,26:ncol(RawM2)],Tranche2RawPlate) 

### Tranche1 Normalisation
TrancheEffect7 = KNNtest(MySoma2[,26:ncol(MySoma2)],Tranche2plate) 

### Tranche1 Normalisation + plate Combat
express1 <- MySoma1[,25:ncol(MySoma1)]
combat_plate_tranche1 <- sva::ComBat(t(express1), PlateBatch1, mod=NULL, par.prior = TRUE, prior.plots = FALSE)
PlateBatchT1 <- as.matrix(MySoma1[,1],ncol=1)  ### Combat and KNNtest function required different class type of PlateBatch
rownames(PlateBatchT1) <- rownames(MySoma1)
TrancheEffect12 = KNNtest(t(combat_plate_tranche1),PlateBatchT1) 

### Tranche1 RawM + plate combat
expressRaw1 <- RawM1[,25:ncol(RawM1)]
combat_plate_trancheRaw1 <- sva::ComBat(t(expressRaw1), PlateBatch1, mod=NULL, par.prior = TRUE, prior.plots = FALSE)
PlateBatchRawT1 <- as.matrix(RawM1[,1],ncol=1)  ### Combat and KNNtest function required different class type of PlateBatch
rownames(PlateBatchRawT1) <- rownames(RawM1)
TrancheEffect13 = KNNtest(t(combat_plate_trancheRaw1),PlateBatchT1) 

### Normalize and Combat plate on tranche 1 

### Tranche1 Raw  
Tranche1RawPlate <-as.matrix(RawM1[,1],ncol=1)
rownames(Tranche1RawPlate) <- rownames(RawM1)
TrancheEffect11 = KNNtest(RawM1[,25:ncol(RawM1)],Tranche1RawPlate) 

### Tranche1 Normalisation
TrancheEffect6 = KNNtest(MySoma1[,25:ncol(MySoma1)],Tranche1plate) 

### Tranche1 Normalisation + plate Combat
express1 <- MySoma1[,25:ncol(MySoma1)]
combat_plate_tranche1 <- sva::ComBat(t(express1), PlateBatch1, mod=NULL, par.prior = TRUE, prior.plots = FALSE)
PlateBatchT1 <- as.matrix(MySoma1[,1],ncol=1)  ### Combat and KNNtest function required different class type of PlateBatch
rownames(PlateBatchT1) <- rownames(MySoma1)
TrancheEffect12 = KNNtest(t(combat_plate_tranche1),PlateBatchT1) 

### Tranche1 RawM + plate combat
expressRaw1 <- RawM1[,25:ncol(RawM1)]
combat_plate_trancheRaw1 <- sva::ComBat(t(expressRaw1), PlateBatch1, mod=NULL, par.prior = TRUE, prior.plots = FALSE)
PlateBatchRawT1 <- as.matrix(RawM1[,1],ncol=1)  ### Combat and KNNtest function required different class type of PlateBatch
rownames(PlateBatchRawT1) <- rownames(RawM1)
TrancheEffect13 = KNNtest(t(combat_plate_trancheRaw1),PlateBatchT1) 

###plate effect on PCA plot: RawM VS norm VS norm+PlateCombat VS RawM+PlateCombat
plot1 <- prcomp(as.matrix(RawM1[,25:ncol(RawM1)]),scale = TRUE)
plot2 <- prcomp(as.matrix(MySoma1[,25:ncol(MySoma1)]),scale = TRUE)
plot3 <- prcomp(as.matrix(t(combat_plate_tranche1)),scale = TRUE)
plot4 <- prcomp(as.matrix(t(combat_plate_trancheRaw1)),scale = TRUE)
pairs(plot1$x[,1:3],col=PlateBatch1)
title("Plate effect of Tranche1 Raw")
pairs(plot2$x[,1:3],col=PlateBatch1)
title("Plate effect of Tranche1 Normalised")
pairs(plot3$x[,1:3],col=PlateBatch1)
title("Plate effect of Tranche1 Normalising+PlateCombat")
pairs(plot4$x[,1:3],col=PlateBatch1)
title("Plate effect of Tranche1 Raw+PlateCombat")


### Tranche2 Raw  
Tranche2RawPlate <-as.matrix(RawM2[,1],ncol=1)
rownames(Tranche2RawPlate) <- rownames(RawM2)
TrancheEffect22 = KNNtest(RawM2[,26:ncol(RawM2)],Tranche2RawPlate) 

### Tranche2 Normalisation
TrancheEffect7 = KNNtest(MySoma2[,26:ncol(MySoma2)],Tranche2plate) 

### Tranche2 Normalisation + plate Combat
express2 <- MySoma2[,26:ncol(MySoma2)]
combat_plate_tranche2 <- sva::ComBat(t(express2), PlateBatch2, mod=NULL, par.prior = TRUE, prior.plots = FALSE)
PlateBatchT2 <- as.matrix(MySoma2[,1],ncol=1)  ### Combat and KNNtest function required different class type of PlateBatch
rownames(PlateBatchT2) <- rownames(MySoma2)
TrancheEffect12 = KNNtest(t(combat_plate_tranche2),PlateBatchT2) 

### Tranche2 RawM + plate combat
expressRaw2 <- RawM2[,26:ncol(RawM2)]
combat_plate_trancheRaw2 <- sva::ComBat(t(expressRaw2), PlateBatch2, mod=NULL, par.prior = TRUE, prior.plots = FALSE)
PlateBatchRawT2 <- as.matrix(RawM2[,1],ncol=1)  ### Combat and KNNtest function required different class type of PlateBatch
rownames(PlateBatchRawT2) <- rownames(RawM2)
TrancheEffect13 = KNNtest(t(combat_plate_trancheRaw2),PlateBatchT2) 

### Tranche 2 Plate Combat + both tranches tranche Combat
ewExp <- rbind(MySoma1[,25:ncol(MySoma1)],t(combat_plate_tranche2))
combat_final <- sva::ComBat(t(NewExp), batchV, mod=NULL, par.prior = TRUE, prior.plots = FALSE)
TrancheEffect8 = KNNtest(t(combat_final),MetaRaw) 
### Tranche 2 Plate Combat + both tranches tranche+plate Combat(can not work)
# mod=model.matrix(~as.factor(PlateBatch))
# combat_final2 <- sva::ComBat(t(NewExp), batchV, mod=mod, par.prior = TRUE, prior.plots = FALSE)
# TrancheEffect9 = KNNtest(t(combat_final2),MetaRaw) 
combat_final2 <- sva::ComBat(t(NewExp), batchV, mod=NULL, par.prior = FALSE, prior.plots = FALSE)
TrancheEffect88 = KNNtest(t(combat_final2),MetaRaw)

### Tranche 2 Plate Combat + both plate tranche Combat
combat_final2 <- sva::ComBat(t(NewExp),PlateBatch, mod=NULL, par.prior = TRUE, prior.plots = FALSE)
TrancheEffect99 = KNNtest(t(combat_final2),MetaRaw) 
### Tranche 2 Plate Combat + both tranches tranche+plate Combat(can not work)
# mod=model.matrix(~as.factor(PlateBatch))
# combat_final2 <- sva::ComBat(t(NewExp), batchV, mod=mod, par.prior = TRUE, prior.plots = FALSE)
# TrancheEffect9 = KNNtest(t(combat_final2),MetaRaw) 
combat_final2 <- sva::ComBat(t(NewExp), batchV, mod=NULL, par.prior = FALSE, prior.plots = FALSE)
TrancheEffect88 = KNNtest(t(combat_final2),MetaRaw)


###plate effect on PCA plot on tranche2: RawM VS norm VS norm+PlateCombat VS RawM+PlateCombat

plot12 <- prcomp(as.matrix(RawM2[,26:ncol(RawM2)]),scale = TRUE)
plot22 <- prcomp(as.matrix(MySoma2[,26:ncol(MySoma2)]),scale = TRUE)
plot32 <- prcomp(as.matrix(t(combat_plate_tranche2)),scale = TRUE)
plot42 <- prcomp(as.matrix(t(combat_plate_trancheRaw2)),scale = TRUE)
pairs(plot12$x[,1:3],col=PlateBatch2)
title("Plate effect of Tranche2 Raw")
pairs(plot22$x[,1:3],col=PlateBatch2)
title("Plate effect of Tranche2 Normalised")
pairs(plot32$x[,1:3],col=PlateBatch2)
title("Plate effect of Tranche2 Normalising+PlateCombat")
pairs(plot42$x[,1:3],col=PlateBatch2)
title("Plate effect of Tranche2 Raw+PlateCombat")

### plate and batch on PCA plot of both tranches
plotF <- prcomp(as.matrix(t(combat_final)),scale = TRUE)
pairs(plotF$x[,1:3],col=batchV)
title("Trache effect after norm tranche2 + tranche combat both")
pairs(plotF$x[,1:3],col=PlateBatch)
title("Plate effect after norm tranche2 + tranche combat both")
TrancheEffect2 = KNNtest(t(combat_all_batch),MetaRaw) 

plotFF <- prcomp(as.matrix(t(combat_final2)),scale = TRUE)
pairs(plotFF$x[,1:3],col=batchV)
title("Trache effect after norm tranche2 + plate combat both")
pairs(plotFF$x[,1:3],col=PlateBatch)
title("Plate effect after norm tranche2 + plate combat both")

RFUdone = t(combat_all_batch)
RFUpc <- prcomp(as.matrix(RFUdone),scale = TRUE)$x[,1:10]

fks = vector(mode="numeric",length=10)
for (myK in 1:10)
{fks[myK] <-FkStatistic(RFUpc,myK)[[1]]}
idK = which(fks<0.85)

pairs(RFUpc[,1:3],col=batchV)
plot(1:10,fks,type="p",main="define cluster number of two tranches data",xlab="cluster number",ylab="f(K) statistic")
lines(1:10,fks,type="l",col="blue")
abline(h=0.85,lty=2)     
grid()
points(idK,fks[idK],pch=23,col="red")

EndoLabel <- NbClust(data = pcDat ,distance = "euclidean", min.nc = 2, max.nc = 10, 
                            method = "kmeans", index = "all", alphaBeale = 0.1)$Best.partition

### try log for combat
express22 <- log(MySoma2[,26:ncol(MySoma2)])
combat_plate_tranche22 <- sva::ComBat(t(express22), PlateBatch2, mod=NULL, par.prior = TRUE, prior.plots = FALSE)
TrancheEffect22 = KNNtest(t(combat_plate_tranche22),PlateBatchT2) 

MySomaAllLog <- rbind(log(MySoma1[,25:5308]),t(combat_plate_tranche22))
combat_log_logComBatTranche2 <- sva::ComBat(t(MySomaAllLog), batchV, mod=NULL, par.prior = TRUE, prior.plots = FALSE)
combat_log_logComBatTranche2_PC <- prcomp(as.matrix(t(combat_log_logComBatTranche2)),scale = TRUE)
temp5 <- combat_log_logComBatTranche2_PC$x[,1:3]
pairs(temp5,col=batchV)
title("Tranche 2 plate combat + Tranche effect after log combat on tranche only")
pairs(temp5,col=PlateBatch)
title("Tranche 2 plate combat + Plate effect after combat on tranche only")
TrancheEffect3 = KNNtest(t(combat_log_logComBatTranche2),MetaRaw) 

combat_log_logComBatTranche22 <- sva::ComBat(t(MySomaAllLog), PlateBatch, mod=NULL, par.prior = FALSE, prior.plots = FALSE)
combat_log_logComBatTranche22_PC <- prcomp(as.matrix(t(combat_log_logComBatTranche22)),scale = TRUE)
temp6 <- combat_log_logComBatTranche22_PC$x[,1:3]
pairs(temp6,col=batchV)
title("Tranche 2 plate combat + Tranche effect after log combat on plate only")
pairs(temp6,col=PlateBatch)
title("Tranche 2 plate combat + Plate effect after combat on Plate only")
TrancheEffect33 = KNNtest(t(combat_log_logComBatTranche22),MetaRaw) 

start_time1 <- Sys.time()
combat_log_logComBatTranche2 <- sva::ComBat(t(MySomaAllLog), batchV, mod=NULL, par.prior = FALSE, prior.plots = FALSE) ###8 min
TrancheEffect44 = KNNtest(t(combat_log_logComBatTranche2),MetaRaw) 

end_time1 <- Sys.time()
start_time2 <- Sys.time()
combat_log_logComBatTranche22 <- sva::ComBat(t(MySomaAllLog), PlateBatch, mod=NULL, par.prior = FALSE, prior.plots = FALSE) ###40min

TrancheEffect44 = KNNtest(t(combat_log_logComBatTranche22),MetaRaw) 
end_time2 <- Sys.time() 

Sys.time()
combat_all_batch <- sva::ComBat(t(MySomaDatAll), batchV, mod=NULL, par.prior = FALSE, prior.plots = FALSE)
combat_all_batch_PC <- prcomp(as.matrix(t(combat_all_batch)),scale = TRUE)
TrancheEffectcombat_all_batch = KNNtest(t(combat_all_batch),MetaRaw) 

Sys.time()
combat_all_batch_log <- sva::ComBat(t(log(MySomaDatAll)), batchV, mod=NULL, par.prior = FALSE, prior.plots = FALSE)
combat_all_batch_log_PC <- prcomp(as.matrix(t(combat_all_batch_log)),scale = TRUE)
TrancheEffectcombat_all_batch_log = KNNtest(t(combat_all_batch_log),MetaRaw) 


Sys.time()
combat_all_plate <- sva::ComBat(t(MySomaDatAll), PlateBatch, mod=NULL, par.prior = FALSE, prior.plots = FALSE)
combat_all_plate_PC <- prcomp(as.matrix(t(combat_all_plate)),scale = TRUE)
TrancheEffectcombat_all_plate = KNNtest(t(combat_all_plate),MetaRaw) 

Sys.time()
combat_all_plate_log <- sva::ComBat(t(log(MySomaDatAll)), PlateBatch, mod=NULL, par.prior = FALSE, prior.plots = FALSE)
combat_all_plate_PC <- prcomp(as.matrix(t(combat_all_plate_log)),scale = TRUE)
TrancheEffectcombat_all_plate_log = KNNtest(t(combat_all_plate_log),MetaRaw) 

Sys.time()
