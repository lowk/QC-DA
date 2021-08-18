### script to filter unquialified proteins and samples before downstream analysis.
library(diptest)
source("QCnormTranche2.R")
source("QCassessTranche2.R")

### filter to exclude buffer|calibrator, select only clinical samples and human proteins
### input RFU matrix after normalisation steps, output filtered RFU matrix only with human data 

filterHM <- function(MySoma,BioMeta){
  HMpro <- which(!grepl("HybControlElution|Non",colnames(MySoma)))
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

### read in prestroed source data: Raw RFU, normalised RFU, and combatted RFU
sourceData <- readRDS("sourceData.rds", refhook = NULL)
RawM1 = sourceData$RawM1
RawM2 = sourceData$RawM2
MySoma1 = sourceData$MySoma1
MySoma2 = sourceData$MySoma2
combat_all = sourceData$combat_all
BioMeta1 = sourceData$BioMeta1
BioMeta2 = sourceData$BioMeta2

###Include only the human samples and huma proteins
MySoma1Done <- filterHM(MySoma1,BioMeta1)[[1]]
MySoma2Done <- filterHM(MySoma2,BioMeta2)[[1]]
BioMeta1Done <- data.frame(filterHM(MySoma1,BioMeta1)[[2]])
BioMeta2Done <- data.frame(filterHM(MySoma2,BioMeta2)[[2]])
tranche <- c(rep(1,nrow(MySoma1Done)),rep(2,nrow(MySoma2Done)))
MySomaPlate <- as.matrix(c(MySoma1Done[,"PlateId"],MySoma2Done[,"PlateId"]),ncol=1)
PlateBatch <- GetPlateBatch(MySomaPlate)
BioMetaM <- data.frame(cbind((rbind(BioMeta1Done,BioMeta2Done)),tranche,PlateBatch))

MySomaRaw <- rbind(MySoma1[,c(1:13,25:ncol(MySoma1))],MySoma2[,c(1:13,26:ncol(MySoma2))])
MySomaAll <- rbind(MySoma1Done[,c(1:13,25:ncol(MySoma1Done))],MySoma2Done[,c(1:13,26:ncol(MySoma2Done))])
exprDat_normStore <- combat_all
exprDat_normMr <- data.frame(cbind(BioMetaM,MySomaAll[,c(1:13)],exprDat_normStore)) ### exprDat_normMr: combined data set with meta clinical information

### OA;disease group filter: set the "clinicType and clinic"
clinicType="OA"
clinic = "OA"
# clinicType="Injury"
# clinic = "INJ"
BioMeta = BioMetaM[which(BioMetaM$diseaseGroup==clinicType),]
exprDat_normM <- exprDat_normMr[which(exprDat_normMr$diseaseGroup==clinicType),]
calib_normM <- exprDat_normMr[grep(paste(clinic,"POOL"),exprDat_normMr$SampleId),]
calib_norm <- as.matrix(calib_normM[,-c(1:(which(colnames(calib_normM)=="CRYBB2.10000.28")-1))])
calibPlates <-  calib_normM$PlateId
exprDat_norm <- exprDat_normM[,-c(1:(which(colnames(exprDat_normM)=="CRYBB2.10000.28")-1))]

##1: Total protein check. Limit of detection.
totalProtein_norm = TotalProCheck(exprDat_norm,BioMeta)
CompM = LoDdetection(MySomaRaw,exprDat_norm) ###CompM > 0, RFU above limit of Detection
# writeLines(paste("total RFUs below limit of Detection is ", length(which(CompM<0))," percentage of total RFUs ",length(which(CompM<0))/(nrow(CompM)*ncol(CompM)),"\n"),conTxt)
paste("total RFUs below limit of Detection is ", length(which(CompM<0))," percentage of total RFUs ",length(which(CompM<0))/(nrow(CompM)*ncol(CompM)),"\n")

##2: Checks against calibrators.
CalibratorCheck(calib_norm,"OA")
CalibratorCheck(calib_norm,"INJ")

# pool variance explained
R2_norm1 = VarExp(calib_normM, calib_norm,"OA",exprDat_norm)
R2_norm2 = VarExp(calib_normM, calib_norm,"INJ",exprDat_norm)

###Check 3: PCA. 
pc_norm <- prcomp(as.matrix(exprDat_norm),scale = TRUE)

###Check 4: Techical confounders
ConfounderTable2 = ConfouderCheck(totalProtein_norm,exprDat_norm,BioMeta) ### confounders against each proteins 
ps_norm = ConfounderTable2[[2]]
apply(ConfounderTable2[[2]],2,function(x){(length(which(x<0.05/5004)))/5004})

###Check 5: multimodel dip test
dipP=vector()
for (dipcounter in 1:ncol(exprDat_norm)){
  dipP[dipcounter] = dip.test(exprDat_norm[,dipcounter], simulate.p.value = FALSE, B = 2000)$p.value
}
length(which(dipP<1e-5))

### summary of remover:
### protein removed: 
Remove1 = which(R2_norm1<0.4)  ### OA/injury repeats
length(Remove1)
Remove2 = which(R2_norm2<0.4)
length(Remove2)

### VarExp for Freeze/Thraw Pool
tempFT <- apply(calib_norm[calib_normM$SampleId %in% paste0("OA POOL-HT-",c(6,25),"/29"),],2,var)/apply(exprDat_norm,2,var)
tempFT[tempFT > 1] <- 1
R2_normFT <- (1 - tempFT)^2
Remove3 = which(R2_normFT<0.4) ### Freeze/Thaw repeats
length(Remove3)
### about R2 explained: !(A | B) = !A & !B
remove1.3 = intersect(Remove1,Remove3)
length(remove1.3)

Remove4 = which(ConfounderTable2[[2]][,"SampleAge"]<1e-5)  ### Anova sample age p
length(Remove4)
Remove5 = which(ConfounderTable2[[2]][,"Plate"]<1e-5)  ### Anova sample age p
length(Remove5)
Remove6 = which(ConfounderTable2[[2]][,"tranche"]<1e-5)  ### Anova sample age p
length(Remove6)

remove7 = which(dipP<1e-5)
length(remove7)

ProteinRatio = apply(CompM,2,function(x){length(which(x<0))/nrow(CompM)})
removeTotalPro = which(ProteinRatio>0.25)
length(removeTotalPro)

removeNameP = colnames(exprDat_norm)[removeTotalPro]
removeName = removeNameP[!grepl("HybControlElution|Non|Spuriomer",removeNameP)]
removeIDprotein = unique(c(remove1.3,Remove4,Remove5,Remove6,remove7,removeTotalPro))
removeIDprotein = unique(c(Remove2,Remove4,Remove5,Remove6,remove7,removeTotalPro))
length(removeIDprotein)
removeProteinName = colnames(exprDat_norm)[removeIDprotein]

###sample removed:
###PCA outliers
myPoints <- pc_norm$x[,1:10]
centroid <- colMeans(myPoints)
distanceToCenter = matrix(NA,ncol=1,nrow=nrow(myPoints))
for (mK in 1: nrow(myPoints)){
  distanceToCenter[mK] <- sqrt(sum((myPoints[mK,]-centroid)^2))
}
PCAout <- which(distanceToCenter>(mean(distanceToCenter) + 5*sd(distanceToCenter)))
length(PCAout)

### total protein outliers
TotalProteinOut <- which(totalProtein_norm>(mean(totalProtein_norm) +5*sd(totalProtein_norm)))
length(TotalProteinOut)
removeIDsample = unique(c(PCAout,TotalProteinOut))
length(removeIDsample)
removedSampName = MySomaAll[rownames(exprDat_norm)[removeIDsample],]$SampleId

filterOADat = exprDat_norm[-removeIDsample,-removeIDprotein]
filterOAMatrix = cbind(exprDat_normM[-removeIDsample,1:23],filterOADat)

filterINJDat = exprDat_norm[-removeIDsample,-removeIDprotein]
filterINJMatrix = cbind(exprDat_normM[-removeIDsample,1:23],filterINJDat)


