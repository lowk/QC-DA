### load in the stored data for convenience
# normList = readRDS("normList.rds")
# RawM1 = normList$RawM1
# RawM2 = normList$RawM2
# TestRaw = normList$TestRaw
# MySoma1Hyb = normList$MySoma1Hyb
# MySoma2Hyb = normList$MySoma2Hyb
# TestHyb = normList$TestHyb
# TestHybCombat = normList$TestHybCombat
# MySoma1Full = normList$MySoma1Full
# MySoma2Full = normList$MySoma2Full
# TestMyFullOnly = normList$TestMyFullOnly
# TestMyFullCombat = normList$TestMyFullCombat
# Soma1 = normList$Soma1
# Soma2 = normList$Soma2
# TestSomaOnly = normList$TestSomaOnly
# TestSomaCombat = normList$TestSomaCombat

library(readxl) 
library(SomaDataIO)
source("QCnormTranche2.R")
library(sva)
library(umap)
library(GGally)
library(factoextra)

load("Luke2021.10.28.RData")

### retrieve pure data block and batch meta data from stored data frame as rds 
CombinedRaw = TestRaw[,which(colnames(TestRaw)=="CRYBB2.10000.28"):ncol(TestRaw)] ###combined expression data
batchMeta_Raw = TestRaw[,c("Tranche Batch","Plate Batch")]

CombinedHyb = TestHyb[,which(colnames(TestHyb)=="CRYBB2.10000.28"):ncol(TestHyb)] 
batchMeta_HybOnly = TestHyb[,c("Tranche Batch","Plate Batch")]

combat_MySomaHyb <- TestHybCombat[,which(colnames(TestHybCombat)=="CRYBB2.10000.28"):ncol(TestHybCombat)]
batchMeta_MySomaHyb = TestHybCombat[,c("Tranche Batch","Plate Batch")]

MyFullOnly <- TestMyFullOnly[,which(colnames(TestMyFullOnly)=="CRYBB2.10000.28"):ncol(TestMyFullOnly)] 
batchMeta_MyFullOnly = TestMyFullOnly[,c("Tranche Batch","Plate Batch")]

combat_MySomaFull <- TestMyFullCombat[,which(colnames(TestMyFullCombat)=="CRYBB2.10000.28"):ncol(TestMyFullCombat)] 
batchMeta_MySomaFul = TestMyFullCombat[,c("Tranche Batch","Plate Batch")]

SomaOnly = TestSomaOnly[,which(colnames(TestSomaOnly)=="CRYBB2.10000.28"):ncol(TestSomaOnly)] 
batchMeta_SomaOnly = TestSomaOnly[,c("Tranche Batch","Plate Batch")]

combat_Soma = TestSomaCombat[,which(colnames(TestSomaCombat)=="CRYBB2.10000.28"):ncol(TestSomaCombat)]  
batchMeta_Soma = TestSomaCombat[,c("Tranche Batch","Plate Batch")]

### set variables Dd and Dk, for the convenience of analysis of different combinations of RFU resources and Dilution group.
Dd=1 ### different RFU resources: 1:7, which are 
# Raw data: RawM1,RawM2, CombinedRaw is the combined raw data.
# Hybnorm only: HubNorm only on seperate tranche individually: MySoma1Hyb, MySoma2Hyb; CombinedHyb is the expression profile of combined hyb normalised data.
# Hybnorm only + ComBat: combat_MySomaHyb is the expression profile by further combatted  on CombinedHyb.
# “Full norm”: full norm on seprate tranche individually: MySoma1Full,MySoma2Full; MyFullOnly is the expression profile of combined full normalised data.
# Full norm + ComBat: combat_MySomaFull is the expression profile by further combatted on MyFullOnly.
# Darryl norm: Daryyl normliased: Soma1,Soma2; SomaOnly is the expression profile of combined Soma1 and Soma2.
# Darryl norm + ComBat: combat_Soma is the expression profile by further combatted on SomaOnly.

# Redo the PCA and UMAP stratifying proteins by:
#   a) dilution group

startDatId = which(colnames(RawM1)=="CRYBB2.10000.28")
Dilution = as.matrix(attributes(RawM1)$Col.Meta[,"Dilution"])
keepCol = which(!grepl("HybControlElution|Non|Spuriomer",colnames(RawM1)[startDatId:ncol(RawM1)]))
DilutionKeep = Dilution[keepCol] ### Dilution is from Col_Meta, need to exclude non human protein
Dilx = unique(DilutionKeep) ### pay attention, unique, levels(as.factor), change the col order
### up to now, DilutionKeep has the same col order as CombinedRaw and all of our stored processed RFUs

### for convenient PCA and UMAP plot, define the following RFUresource, batchResource,titleMessage
FrameResource = list(TestRaw,TestHyb,TestHybCombat,TestMyFullOnly,TestMyFullCombat,TestSomaOnly,TestSomaCombat)
RFUresource = list(CombinedRaw,CombinedHyb,TestHybCombat,MyFullOnly,combat_MySomaFull,SomaOnly,combat_Soma)
batchResource = list(batchMeta_Raw,batchMeta_HybOnly,combat_MySomaHyb,batchMeta_MyFullOnly,batchMeta_MySomaFul,batchMeta_Soma)
titleMessage = c("Raw combined RFU", "Hyb combined RFU", "Hyb combined + Combat corrected RFU", "full normalised combined RFU",
                 "full normalised combined + Combat corrected", "Soma normliased combined RFU", "Soma normliased + Combat corrected RFU")

### I have a question about this variable - LJD
Dk=1 ### different dilution concentration group 1:3

### PCA/UMAP plot based on the dilution group "Dk" from the RFUresource of "Dd"

PlotPCA(RFUresource[[Dd]][,which(DilutionKeep==Dilx[Dk])],batchResource[[Dd]],3,1,paste("PCA on ", titleMessage[[Dd]], "from dilution ",Dilx[Dk])) 
PlotPCA(RFUresource[[Dd]][,which(DilutionKeep==Dilx[Dk])],batchResource[[Dd]],3,2,paste("PCA on ", titleMessage[[Dd]], "from dilution ",Dilx[Dk])) 

PlotUmap(RFUresource[[Dd]][,which(DilutionKeep==Dilx[Dk])],batchResource[[Dd]],1,paste("UMAP on ", titleMessage[[Dd]], "from dilution ",Dilx[Dk])) 
PlotUmap(RFUresource[[Dd]][,which(DilutionKeep==Dilx[Dk])],batchResource[[Dd]],1,paste("UMAP on ", titleMessage[[Dd]], "from dilution ",Dilx[Dk])) 

###b) cutting off highly expressed proteins
totalPro = colSums(RFUresource[[Dd]])
exProId = totalPro>(mean(totalPro)+2*sd(totalPro)) ### define protein to be exluded by mean+2*sd of total expression level for each protein
#exProId = totalPro>(mean(totalPro)+2*sd(totalPro)) | totalPro<(mean(totalPro)-2*sd(totalPro))  ### define protein to be exluded by mean+2*sd of total expression level for each protein
ProKeep = RFUresource[[Dd]][,!exProId] ### proteins to keep for plot

### PCA/UMAP 
PlotPCA(ProKeep,batchResource[[Dd]],3,1,paste("Excluding highly expressed proteins, PCA on ", titleMessage[[Dd]])) 
PlotPCA(ProKeep,batchResource[[Dd]],3,2,paste("Excluding highly expressed proteins, PCA on ", titleMessage[[Dd]])) 

PlotUmap(ProKeep,batchResource[[Dd]],1,paste("Excluding highly expressed proteins, UMAP on ", titleMessage[[Dd]])) 
PlotUmap(ProKeep,batchResource[[Dd]],2,paste("Excluding highly expressed proteins, UMAP on ", titleMessage[[Dd]])) 

###c) using the "strict" QC proteins
# strictList<- read.csv("/Users/ydeng/Documents/QCstepOA/normComp/protein_filters.txt",sep="\t") ### directly read in proteins pass strict filters
# strictPro <- strictList[which(strictList$inj_filter_status=="Pass"&strictList$oa_filter_status=="Pass"),1]
strictRFU <- RFUresource[[Dd]][,strictPro]### RFU frames from each resource with strictly filtered proteins.

### PCA/UMAP 
PlotPCA(strictRFU,batchResource[[Dd]],3,1,paste("Strictly filtered proteins, PCA on ", titleMessage[[Dd]])) 
PlotPCA(strictRFU,batchResource[[Dd]],3,2,paste("Strictly filtered proteins, PCA on ", titleMessage[[Dd]])) 
PlotUmap(strictRFU,batchResource[[Dd]],1,paste("Strictly filtered proteins, UMAP on ", titleMessage[[Dd]])) 
PlotUmap(strictRFU,batchResource[[Dd]],2,paste("Strictly filtered proteins, UMAP on ", titleMessage[[Dd]])) 

### comprehensive strictly filter code coming soon.

# Make plots of shift (i.e. ie plot changes in mean of each protein between tranche 1 and tranche 2) for:
#   a) plasma calibrators
### RawM1 vs RawM2
Raw1mean <- meanShift(RawM1,RawM2,"Calibrator")[[1]]
Raw2mean <- meanShift(RawM1,RawM2,"Calibrator")[[2]]    
hist((Raw2mean-Raw1mean)/Raw1mean)
### hybnorm only: MyHyb1 vs MyHub2
Hyb1mean <- meanShift(MySoma1Hyb,MySoma2Hyb,"Calibrator")[[1]]
Hyb2mean <- meanShift(MySoma1Hyb,MySoma2Hyb,"Calibrator")[[2]]    
hist((Hyb2mean-Hyb1mean)/Hyb1mean)
### full norm: MyFull1 vs MyFull2
Full1mean <- meanShift(MySoma1Full,MySoma2Full,"Calibrator")[[1]]
Full2mean <- meanShift(MySoma1Full,MySoma2Full,"Calibrator")[[2]]    
hist((Full2mean-Full1mean)/Full1mean)
### Soma Nomr: Soma1 vs Soma2
Soma1mean <- meanShift(Soma1,Soma2,"Calibrator")[[1]]
Soma2mean <- meanShift(Soma1,Soma2,"Calibrator")[[2]]    
hist((Soma2mean-Soma1mean)/Soma1mean)

#   c) buffers, for simplicity, use the same variables as "Calibrator" mean shift
### RawM1 vs RawM2
Raw1mean <- meanShift(RawM1,RawM2,"Buffer")[[1]]
Raw2mean <- meanShift(RawM1,RawM2,"Buffer")[[2]]    
hist((Raw2mean-Raw1mean)/Raw1mean)
### hybnorm only: MyHyb1 vs MyHub2
Hyb1mean <- meanShift(MySoma1Hyb,MySoma2Hyb,"Buffer")[[1]]
Hyb2mean <- meanShift(MySoma1Hyb,MySoma2Hyb,"Buffer")[[2]]    
hist((Hyb2mean-Hyb1mean)/Hyb1mean)
### full norm: MyFull1 vs MyFull2
Full1mean <- meanShift(MySoma1Full,MySoma2Full,"Buffer")[[1]]
Full2mean <- meanShift(MySoma1Full,MySoma2Full,"Buffer")[[2]]    
hist((Full2mean-Full1mean)/Full1mean)
### Soma Nomr: Soma1 vs Soma2
Soma1mean <- meanShift(Soma1,Soma2,"Buffer")[[1]]
Soma2mean <- meanShift(Soma1,Soma2,"Buffer")[[2]]    
hist((Soma2mean-Soma1mean)/Soma1mean)

# b) injury and oa pools, for simplicity, use the same variables as "Calibrator" mean shift
### RawM1 vs RawM2
Raw1mean <- meanShiftDis(RawM1,RawM2,"OA POOL")[[1]]
Raw2mean <- meanShiftDis(RawM1,RawM2,"OA POOL")[[2]]    
hist((Raw2mean-Raw1mean)/Raw1mean)
Raw1mean <- meanShiftDis(RawM1,RawM2,"INJ POOL")[[1]]
Raw2mean <- meanShiftDis(RawM1,RawM2,"INJ POOL")[[2]]    
hist((Raw2mean-Raw1mean)/Raw1mean)
### hybnorm only: MyHyb1 vs MyHub2
Hyb1mean <- meanShiftDis(MySoma1Hyb,MySoma2Hyb,"OA POOL")[[1]]
Hyb2mean <- meanShiftDis(MySoma1Hyb,MySoma2Hyb,"OA POOL")[[2]]    
hist((Hyb2mean-Hyb1mean)/Hyb1mean)
Hyb1mean <- meanShiftDis(MySoma1Hyb,MySoma2Hyb,"INJ POOL")[[1]]
Hyb2mean <- meanShiftDis(MySoma1Hyb,MySoma2Hyb,"INJ POOL")[[2]]    
hist((Hyb2mean-Hyb1mean)/Hyb1mean)
### full norm: MyFull1 vs MyFull2
Full1mean <- meanShiftDis(MySoma1Full,MySoma2Full,"OA POOL")[[1]]
Full2mean <- meanShiftDis(MySoma1Full,MySoma2Full,"OA POOL")[[2]]    
hist((Full2mean-Full1mean)/Full1mean)
Full1mean <- meanShiftDis(MySoma1Full,MySoma2Full,"INJ POOL")[[1]]
Full2mean <- meanShiftDis(MySoma1Full,MySoma2Full,"INJ POOL")[[2]]    
hist((Full2mean-Full1mean)/Full1mean)
### Soma Nomr: Soma1 vs Soma2
Soma1mean <- meanShiftDis(Soma1,Soma2,"OA POOL")[[1]]
Soma2mean <- meanShiftDis(Soma1,Soma2,"OA POOL")[[2]]
hist((Soma2mean-Soma1mean)/Soma1mean)
Soma1mean <- meanShiftDis(Soma1,Soma2,"INJ POOL")[[1]]
Soma2mean <- meanShiftDis(Soma1,Soma2,"INJ POOL")[[2]]
hist((Soma2mean-Soma1mean)/Soma1mean)

# f) freshly treated OA pool 
### freshly treated OA pool vs pre treated OA pool within tranche2
meanOA1 = colMeans(RawM2[grepl("OA POOL",RawM2$SampleId) & RawM2$SampleId!="OA POOL-HT-26/29",which(colnames(RawM2)=="CRYBB2.10000.28"):ncol(RawM2)])
meanOA2 = as.matrix(RawM2[which(RawM2$SampleId =="OA POOL-HT-26/29"),which(colnames(RawM2)=="CRYBB2.10000.28"):ncol(RawM2)])
hist((meanOA1-meanOA2)/meanOA2)
### freshly treated OA pool vs pre treated OA pool within tranche1
meanOA1 = colMeans(RawM1[grepl("OA POOL",RawM1$SampleId),which(colnames(RawM1)=="CRYBB2.10000.28"):ncol(RawM1)])
meanOA2 = as.matrix(RawM2[which(RawM2$SampleId =="OA POOL-HT-26/29"),which(colnames(RawM2)=="CRYBB2.10000.28"):ncol(RawM2)])
hist((meanOA1-meanOA2)/meanOA2)

###d) oa samples; e) injury samples
# clinicFile1 <- "/Users/ydeng/Documents/QCstepOA/normComp/STEpUP_QCData_Tranche1.xlsx"
# clinicFile2 <- "/Users/ydeng/Documents/QCstepOA/normComp/STEpUP_QCData_Tranche2_09Mar2021.xlsx"  
# ClinicMeta1 <- ExtractClinicG(RawM1,clinicFile1,1)
# ClinicMeta2 <- ExtractClinicG(RawM2,clinicFile2,2)
ClinicMeta <- rbind(ClinicMeta1,ClinicMeta2)
# MetaWhole <- cbind(ClinicMeta,rbind(RawM1[,-(23:24)],RawM2[,-(23:25)]))
### ClinicMeta has the same row order as RawM1 and RawM2, but not definitely the same as MySoma and Soma
### MetaWhole: most comprehensive meta data frame up to date 26 OCT 2021. Using self written function: ExtractClinicG + addCLinicMeta, make sure the order of combined meta have consistency.

### combine CLinicMeta to combined and combat corrected RFU frames. Pay attention to the row order for each clinicMeta info combination.
ClinicRFU <- addCLinicMeta(ClinicMeta,FrameResource[[Dd]])

### mean shift for OA patients between two tranches
meanRFU1 <- meanShiftClinic(ClinicRFU,"OA")[[1]]
meanRFU2 <- meanShiftClinic(ClinicRFU,"OA")[[2]]
hist((meanRFU2-meanRFU1)/meanRFU1)
### mean shift for injury patients between two tranches
meanRFU1 <- meanShiftClinic(ClinicRFU,"Injury")[[1]]
meanRFU2 <- meanShiftClinic(ClinicRFU,"Injury")[[2]]
hist((meanRFU2-meanRFU1)/meanRFU1)

#Investigate plate-based covariates to see what drives the difference between plates within tranche 2;
# - Date of run
# - Scanner ID
ClinicRFU <- addCLinicMeta(ClinicMeta2,RawM2)
RFUdone <- filterHM2(ClinicRFU)[[1]]
clinicMetaDone <- filterHM2(ClinicRFU)[[2]]
pc_norm <- prcomp(as.matrix(RFUdone),scale = TRUE)
PlotDat = data.frame(cbind(pc_norm$x[,1:3],clinicMetaDone))

ggpairs(PlotDat, columns=1:3, aes(color= as.factor(ScannerID)),
        title= paste("RawM2 from ",titleMessage[[Dd]]),
        diag=list(continuous=wrap("densityDiag",alpha=0.4)),
        lower=list(continuous = wrap("points",alpha=0.9,size=0.1)),
        upper = list(continuous = "blank"),
        legend = c(1,1)) + labs(fill = "ScannerID",theme(legend.title=element_text(size=18),legend.text=element_text(size=18),axis.title=element_text(size=18,face="bold")))

ggpairs(PlotDat, columns=1:3, aes(color= as.factor(PlateRunDate)),
        title= paste("RawM2 from ",titleMessage[[Dd]]),
        diag=list(continuous=wrap("densityDiag",alpha=0.4)),
        lower=list(continuous = wrap("points",alpha=0.9,size=0.1)),
        upper = list(continuous = "blank"),
        legend = c(1,1)) + labs(fill = "Date of Run",theme(legend.title=element_text(size=18),legend.text=element_text(size=18),axis.title=element_text(size=18,face="bold")))

ggpairs(PlotDat, columns=1:3, aes(color= as.factor(diseaseGroup)),
        title= paste("RawM2 from ",titleMessage[[Dd]]),
        diag=list(continuous=wrap("densityDiag",alpha=0.4)),
        lower=list(continuous = wrap("points",alpha=0.9,size=0.1)),
        upper = list(continuous = "blank"),
        legend = c(1,1)) + labs(fill = "Date of Run",theme(legend.title=element_text(size=18),legend.text=element_text(size=18),axis.title=element_text(size=18,face="bold")))


### it is convenient to make plot PCA referring to any confounders included in PlotDat here.
### code is easy to check tranche1 data as well, just change the following one line
ClinicRFU <- addCLinicMeta(ClinicMeta1,RawM1)
RFUdone <- filterHM2(ClinicRFU)[[1]]
clinicMetaDone <- filterHM2(ClinicRFU)[[2]]
pc_norm <- prcomp(as.matrix(RFUdone),scale = TRUE)
PlotDat = data.frame(cbind(pc_norm$x[,1:3],clinicMetaDone))

### plate batch effect of tranche2 reason:
library(mclust)
similairyKmHc <- mclust::adjustedRandIndex(ClinicRFU$PlateID,ClinicRFU$PlateRunDate)
similairyKmHc <- mclust::adjustedRandIndex(ClinicRFU$PlateID,ClinicRFU$ScannerID)
similairyKmHc <- mclust::adjustedRandIndex(ClinicRFU$PlateID,ClinicRFU$diseaseGroup)
similairyKmHc

library(ggridges)
MySomaShift = Full2mean/Full1mean
SomaShift = Soma1mean/Soma2mean
plotMean = data.frame(t(rbind(c(rep("MySomaShift",length(MySomaShift)),rep("SomaShift",length(SomaShift))),
                              c(MySomaShift,SomaShift),seq(1:(length(MySomaShift)+length(SomaShift))))))
colnames(plotMean)=c("Type","Expression","seqID")
ggplot(plotMean, aes(x = as.numeric(Expression), y = Type)) +
  geom_density_ridges() +
  theme_ridges() + 
  theme(legend.position = "none")

sPlate <- unique(MetaWhole$PlateId)[[7]]
sPlateFrame <- MetaWhole[which(MetaWhole$PlateId==sPlate),]
sPlatePlot <- cbind(rep(1:12,each=8),rep(1:8,12),sPlateFrame[,1:24])
colnames(sPlatePlot)[1:2] <-c("y","x")

ggplot(sPlatePlot, aes(x = x, y = y)) + 
  geom_point(aes(color = Corhort), alpha = 0.5,size=14) 
