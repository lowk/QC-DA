library(readxl) 
library(SomaDataIO)
source("QCnormTranche2.R")
library(sva)

### change 5 input parameters here: Funlist (select from "HYBNORM,PLATESCALE,MIDNORM/MIDNORMcali/MIDNORMsamp,CALIBRATION". The order is the aimed normalisation procedure)
# 1.Funlist = list(RawM)
# 2.Funlist = list(HYBNORM)
# 3.Funlist = list(HYBNORM,MIDNORMcali,PLATESCALE) 
# 4.Funlist = list(HYBNORM,MIDNORMcali,PLATESCALE,MIDNORMsamp)
# 5.Funlist = list(HYBNORM,MIDNORMcali,PLATESCALE,MIDNORMsamp,CALIBRATION)
# 6.Funlist = list(HYBNORM,MIDNORMcali,PLATESCALE,CALIBRATION)

### tranche1 and 2, use the following normalisation steps within each tranche
Funlist = list(HYBNORM,MIDNORMcali,PLATESCALE,CALIBRATION)

############ Normalise two tranches data seperately; extract col table information from adat file#############################
inputfile11 <- paste(myFilePath,"SS-200008.ADat",sep="")
inputfile12 <- paste(myFilePath,"SS-200008.hybNorm.medNormInt.plateScale.medNormRefSMP.ADat",sep="")
inputfile13 <- paste(myFilePath,"SS-200008.hybNorm.medNormInt.plateScale.ADat",sep="")
inputfile21 <- paste(myFilePath,"SS-205086.adat",sep="")
inputfile22 <- paste(myFilePath,"SS-205086.hybNormRef.medNormInt.medNormRefSMP.adat",sep="")
inputfile23 <- paste(myFilePath,"SS-205086.hybNormRef.medNormInt.adat",sep="")
inputfile3 <- paste(myFilePath,"/clinic/STEpUP_QCData_Tranche1.xlsx",sep="")
inputfile4 <- paste(myFilePath,"/clinic/STEpUP_QCData_Tranche2_09Mar2021.xlsx",sep="")

RawM1 <- read.adat(inputfile11)
ColTable1 <- attr(read.adat(inputfile12),"Col.Meta") ### ColTable for the convenience of further enrichment and network analysis
MySoma1 = UserNorm(Funlist,RawM1) ## user select which normlisation methods.

RawM2 <- read.adat(inputfile21)
ColTable2 <- attr(read.adat(inputfile22),"Col.Meta")
MySoma2 = UserNorm(Funlist,RawM2)

### check code accuracy by comparing our normalized data VS somalogic provided data
Flist = list(HYBNORM,MIDNORMcali)
checkMyCode(RawM2,inputfile23,Flist)
Flist = list(HYBNORM,MIDNORMcali,PLATESCALE,MIDNORMsamp)
checkMyCode(RawM1,inputfile12,Flist)
checkMyCode(RawM2,inputfile22,Flist) ### only 9 HybControlElution has correlation coefficient < 0.7

############ initialise information from excel files ############################
BioMeta1 <- data.frame(ExtractClinicG(RawM1,inputfile3,1))
BioMeta2 <- data.frame(ExtractClinicG(RawM2,inputfile4,2))

###adjust BioMeta row order to match MySoma
idADJ = vector("numeric",length=nrow(MySoma2))
for (kk in 1:nrow(MySoma2)){
  idADJ[kk] <- which(rownames(BioMeta2) == rownames(MySoma2)[kk])
}
BioMeta2 <- BioMeta2[idADJ,]
BioMeta2[,"patientGender"] <- sub("Male","M",BioMeta2[,"patientGender"])
BioMeta2[,"patientGender"] <- sub("Female","F",BioMeta2[,"patientGender"])
##restored BioMeta has completed the above order match and gender symbol replacement.
## check rownames and colnames:
all(rownames(BioMeta1)==rownames(MySoma1))
all(rownames(BioMeta2)==rownames(MySoma2))
all(colnames(MySoma1) == colnames(RawM1))
all(colnames(MySoma2) == colnames(RawM2))


############ combat correction for plate batch effect and tranche batch effect #############################
###Include only the human samples and human proteins
MySoma1Done <- filterHM(MySoma1,BioMeta1)[[1]]
MySoma2Done <- filterHM(MySoma2,BioMeta2)[[1]]
BioMeta1Done <- data.frame(filterHM(MySoma1,BioMeta1)[[2]])
BioMeta2Done <- data.frame(filterHM(MySoma2,BioMeta2)[[2]])

MySomaAllDat <- rbind(MySoma1Done[,which(colnames(MySoma1Done) == "CRYBB2.10000.28"):ncol(MySoma1Done)],MySoma2Done[,which(colnames(MySoma2Done) == "CRYBB2.10000.28"):ncol(MySoma2Done)])
tranche <- c(rep(1,nrow(MySoma1Done)),rep(2,nrow(MySoma2Done)))
MySomaPlate <- as.matrix(c(MySoma1Done[,"PlateId"],MySoma2Done[,"PlateId"]),ncol=1)
PlateBatch <- GetPlateBatch(MySomaPlate)
MetaRaw = cbind(tranche,PlateBatch)
MySomaAll <- rbind(MySoma1Done[,c(1:13,which(colnames(MySoma1Done) == "CRYBB2.10000.28"):ncol(MySoma1Done))],MySoma2Done[,c(1:13,which(colnames(MySoma2Done) == "CRYBB2.10000.28"):ncol(MySoma2Done))])
rownames(MetaRaw)=rownames(MySomaAll)
colnames(MetaRaw) = c("Tranche Batch","Plate Batch")
PlateBatch <- as.vector(PlateBatch) ### Combat argument requirement

combat_all_batch <- sva::ComBat(t(log(MySomaAllDat)), PlateBatch, mod=NULL, par.prior = TRUE, prior.plots = FALSE)
combat_all = data.frame(t(combat_all_batch))

### k nearest neiborhood test on combatted data
BioMetaDone <- data.frame(cbind((rbind(BioMeta1Done,BioMeta2Done)),tranche,PlateBatch))
TrancheEffect <- KNNtest(t(combat_all_batch),BioMetaDone[,9:10],2,2)
TrancheEffect

############ now save Raw RFU, normliased RFU, and combatted RFU in "sourceData.rds"
sourceData = list("RawM1"=RawM1,"RawM2"=RawM2,"MySoma1"=MySoma1,"MySoma2"=MySoma2,"combat_all"=combat_all,"BioMeta1"=BioMeta1,"BioMeta2"=BioMeta2)
saveRDS(sourceData,file="sourceData.rds")
