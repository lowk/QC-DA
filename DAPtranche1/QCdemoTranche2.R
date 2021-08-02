### mainbody
clr  = function()
{
  rm(list = ls(envir = .GlobalEnv),
     envir = .GlobalEnv)
}
clr()

###main 
library(readxl) 
library(factoextra)
library(gridExtra)
library(SomaDataIO)
library(GGally)
library(e1071)
library(Rtsne)
library(Hmisc)
library(moments)
source("QCnormTranche2.R")
source("QCassessTranche2.R")


### change 5 input parameters here: Funlist (select from "HYBNORM,PLATESCALE,MIDNORM/MIDNORMcali/MIDNORMsamp,CALIBRATION". The order is the aimed normalisation procedure)
### choose the output files with remindable names, such as "tableOut.txt","MySomaPlot.pdf","ps_pca_norm.txt","ps_norm.txt"
# 1.Funlist = list(RawM)
# 2.Funlist = list(HYBNORM)
# 3.Funlist = list(HYBNORM,MIDNORMcali,PLATESCALE) 
# 4.Funlist = list(HYBNORM,MIDNORMcali,PLATESCALE,MIDNORMsamp)
# 5.Funlist = list(HYBNORM,MIDNORMcali,PLATESCALE,MIDNORMsamp,CALIBRATION)
# 6.Funlist = list(HYBNORM,MIDNORMcali,PLATESCALE,CALIBRATION)

Funlist = list(HYBNORM,MIDNORMcali,PLATESCALE,CALIBRATION)
# MySomaPlot = "MySomaPlot1.pdf"
# tableOut = "tableOutTestxxxx.csv"
# ps_pca_normOut = "ps_pca_normTest1.csv"
# ps_normOut = "ps_normTest1.csv"


#part1#######################################################################################################
myFilePath <- "/Users/ydeng/Documents/QCstepOA/"  ### set required file directory 
inputfile1 <- paste(myFilePath,"SS-205086.adat",sep="")
inputfile2 <- paste(myFilePath,"SS-205086.hybNormRef.medNormInt.medNormRefSMP.adat",sep="")
inputfile3 <- paste(myFilePath,"/clinic/STEpUP_QCData_Tranche2_09Mar2021.xlsx",sep="")
initiaList<- initQCnorm(inputfile1,inputfile2) ###raw RFUs from Adat, ajust SampleType according to tranch Excel
RawM2 <- initiaList[[1]]
ColTable <- initiaList[[2]] ### ColTable for the convenience of further enrichment and network analysis

metadata <- read_excel(inputfile3,range="A2:O618",col_names=TRUE)
### user select which normlisation methods. 
MySoma2 = UserNorm(Funlist,RawM2)
# SomaM <- read.adat(inputfile2)
# corShip = CompTWO(SomaM,MySoma2) ### comparison between our normalisation and Adat
# DatStartId2 <- which(colnames(MySoma2)=="RMA")+1
# colnames(MySoma2[DatStartId2:ncol(MySoma2)])[which(corShip<0.6)] ##HybControlElution highly different

#part2######################################################################################################
### initialise information from excel files
RawMList <- ExtractClinicG(RawM2,inputfile3)
RawM2 <- RawMList[[1]]
BioMeta2 <- RawMList[[2]]



###MySomaAll: combined tranches samples with plate id ect. MySomaRaw:combined tranches with calibrator and buffer; exprDat_norm: = combined tranches sample hunman protein expression data 
### exprDat_normM:combat human samples with plateid ect.
# koad("BatchCorrection.RData")
# load("BioMeta.RData") ### this file stores BioMeta1 and BioMeta2, rownames matching MySoma1 and MySoma2.
###adjust BioMeta row order to match MySoma
idADJ = vector("numeric",length=nrow(MySoma2))
for (kk in 1:nrow(MySoma2)){
  idADJ[kk] <- which(rownames(BioMeta2) == rownames(MySoma2)[kk])
}
BioMeta2 <- BioMeta2[idADJ,]
BioMeta2[,"patientGender"] <- sub("Male","M",BioMeta2[,"patientGender"])
BioMeta2[,"patientGender"] <- sub("Female","F",BioMeta2[,"patientGender"])
###Include only the human samples and huma proteins, and perform the same analysis from main body line 132
MySoma1Done <- filterHM(MySoma1,BioMeta1)[[1]]
MySoma2Done <- filterHM(MySoma2,BioMeta2)[[1]]
BioMeta1Done <- filterHM(MySoma1,BioMeta1)[[2]]
BioMeta2Done <- filterHM(MySoma2,BioMeta2)[[2]]

MySomaAllDat <- rbind(MySoma1Done[,25:ncol(MySoma1Done)],MySoma2Done[,26:ncol(MySoma2Done)])
tranche <- c(rep(1,nrow(MySoma1Done)),rep(2,nrow(MySoma2Done)))
MySomaPlate <- as.matrix(c(MySoma1Done[,"PlateId"],MySoma2Done[,"PlateId"]),ncol=1)
PlateBatch <- GetPlateBatch(MySomaPlate)
MetaRaw = cbind(tranche,PlateBatch)
rownames(MetaRaw)=rownames(MySomaAll)
colnames(MetaRaw) = c("Tranche Batch","Plate Batch")
PlateBatch <- as.vector(PlateBatch) ### Combat argument requirement

MySomaRaw <- rbind(MySoma1[,c(1:13,25:ncol(MySoma1))],MySoma2[,c(1:13,26:ncol(MySoma2))])
MySomaAll <- rbind(MySoma1Done[,c(1:13,25:ncol(MySoma1Done))],MySoma2Done[,c(1:13,26:ncol(MySoma2Done))])
exprDat_norm <- t(combat_all_batch6)
exprDat_normM <- cbind(MySomaAll[,c(1:13)],exprDat_norm) 
calib_normM <- data.frame(exprDat_normM[grep("POOL",exprDat_normM$SampleId),])
calib_norm <- as.matrix(calib_normM[,-c(1:(which(colnames(calib_normM)=="CRYBB2.10000.28")-1))])
calibIDs <-  calib_normM$SampleId[grep("POOL",calib_normM$SampleId)]  ###calibID corresponding to calib_norm
calibPlates <-  calib_normM$PlateId[grep("POOL",calib_normM$SampleId)]

BioMeta <- data.frame(rbind(BioMeta1Done,BioMeta2Done))

### adjust BioMeta row order matching exprDat_norm. RawM2 row order different from MySoma2, due to ###Tranche2: levels(as.factor): plate order is this: "P0028435" "P0028436" "P0028443" "P0028444" "P0028445" "P0028446" "P0028453"
###But in the adat file, plate order is: "P0028435" "P0028436" "P0028443" "P0028444" "P0028453" "P0028445" "P0028446"
### always check the row order of meta data and raw and MySoma

# ### calibrator tranche2 only
# calib_norm2 <- data.frame(MySoma2[,!grepl("HybControlElution|Non",colnames(MySoma2))])
# calib_norm2<- as.matrix(calib_norm2[grep("POOL",MySoma2$SampleId),-c(1:which(colnames(MySoma2)=="CRYBB2.10000.28"))])
# calibIDs2 <-  MySoma2$SampleId[grep("POOL",MySoma2$SampleId)]  ###calibID corresponding to calib_norm
# calibPlates2 <-  MySoma2$PlateId[grep("POOL",MySoma2$SampleId)]

#part3########################################################################################################
### QC assess begins. Based on different normalised MySoma.
# pdf(MySomaPlot)
# par(mfrow=c(1,3))
# conTxt = file(tableOut,"append")

##1: Total protein check. Limit of detection.
totalProtein_norm = TotalProCheck(exprDat_norm,BioMeta)
CompM = LoDdetection(MySomaRaw) ###CompM > 0, RFU above limit of Detection
# writeLines(paste("total RFUs below limit of Detection is ", length(which(CompM<0))," percentage of total RFUs ",length(which(CompM<0))/(nrow(CompM)*ncol(CompM)),"\n"),conTxt)
paste("total RFUs below limit of Detection is ", length(which(CompM<0))," percentage of total RFUs ",length(which(CompM<0))/(nrow(CompM)*ncol(CompM)),"\n")

##2: Checks against calibrators.
CalibratorCheck(calib_norm,"OA")
CalibratorCheck(calib_norm,"INJ")

# pool variance explained
R2_norm1 = VarExp(calib_norm,"OA")
R2_norm2 = VarExp(calib_norm,"INJ")

##CV breakdowns:
printListDone1 = CVbreak(calib_norm,"OA")
writeLines("CVbreakdown for OA group:",conTxt)
write.table(printListDone1,conTxt,row.names = FALSE,col.names=FALSE)
printListDone2 = CVbreak(calib_norm,"INJ")
writeLines("CVbreakdown for INJ group:",conTxt)
write.table(printListDone2,conTxt,row.names = FALSE,col.names=FALSE)

###Check 3: PCA. 
listPCA = PCAglob(MetaRawM,exprDat_norm)
pc_norm = listPCA[[1]] ###pc_norm: coordinates in first 10 principles 
listPCA[[2]] ###plot: PCs against corresponding eigen values 
listPCA[[3]] ###plot: variance explained by PCs


###Check 4: Techical confounders
ConfounderTable1 = ConfouderCheck(totalProtein_norm,pc_norm$x[,1:10],BioMeta) ### confounders against PCs 
ps_pca_norm = ConfounderTable1[[2]]
writeLines("\ntechnical counfoundr check in order: total protein,ps_pca_norm,varExp_norm",conTxt)
for (ct in 1:3) {write.table(ConfounderTable1[[ct]], conTxt)}
write.csv(ps_pca_norm, file = ps_pca_normOut) ### restore a ps_norm formal file, for direct comparison within one graph between two MySoma

ConfounderTable2 = ConfouderCheck(totalProtein_norm,exprDat_norm,BioMeta) ### confounders against each proteins 
ps_norm = ConfounderTable2[[2]]
writeLines("\ntechnical counfoundr check in order: ps_norm,varExp_norm",conTxt)
for (ct in 2:3) {write.table(ConfounderTable2[[ct]], conTxt)} 
write.csv(ps_norm, file = ps_normOut) 

confounderPlot1(ConfounderTable1)  ###confounder vs total protein
confounderPlot2(ConfounderTable1,0) ### confounder vs per PC
confounderPlot2(ConfounderTable2,1) ### confounder vs per protein

apply(ConfounderTable2[[2]],2,function(x){(length(which(x<0.05/5004)))/5004})

###Check 5: External validation
toTest1 <- c("mcp1bl","il6bl","il8bl","mmp3bl","activinabl","tsg6bl","timp1bl","tgfb1bl","fgf2bl")
toTest2 <- c("CCL2.2578.67","IL6.4673.13","CXCL8.3447.64","MMP3.2788.55","INHBA.13738.8","TNFAIP6.5036.50","TIMP1.2211.9","TGFB1.2333.72","FGF2.3025.50")

###check Ben (this is all OA) / Historic data (this is all knee injury). Return correlation between RFUs and datas from another equipment
CorData_norm1 = ExtVal1(sandwich_master_Ben,exprDat_norm,toTest1,toTest2,conTxt)
writeLines("\nsandwich_master_Ben:",conTxt)
write.table(CorData_norm1,conTxt)
CorData_norm2 = ExtVal1(sandwich_master_Historic,exprDat_norm,toTest1,toTest2,conTxt)
writeLines("\nsandwich_master_Historic:",conTxt)
write.table(CorData_norm2,conTxt)

CorData_plate1 = ExtVal2(sandwich_master_Ben,"sandwich_Ben",exprDat_norm,metadata_reord,plateID,toTest1,toTest2,conTxt) 
writeLines("\nPlate sandwich_master_Ben:",conTxt)
write.table(CorData_plate1,conTxt,row.names = FALSE,col.names=FALSE)
CorData_plate2 = ExtVal2(sandwich_master_Historic,"sandwich_Historic",exprDat_norm,metadata_reord,plateID,toTest1,toTest2,conTxt) 
writeLines("\nPlate sandwich_master_Historic:",conTxt)
write.table(CorData_plate2,conTxt,row.names = FALSE,col.names=FALSE)

###check the HT treated data (5 OA and 5 injury)
toTest3 = c("CD14","CD14","VEGFA","VEGFA","VEGFA")
toTest4 = c("CD14.16914.104","CD14.8969.49","VEGFA.19437.61","VEGFA.2597.8","VEGFA.4867.15")
HAse_UT = HTcheck(sandwich_master_HAse_UT,exprDat_norm,toTest3,toTest4) 
writeLines("\nsandwich_master_HAse_UT:",conTxt)
write.table(HAse_UT,conTxt)
HAse_HT = HTcheck(sandwich_master_HAse_HT,exprDat_norm,toTest3,toTest4)
writeLines("\nsandwich_master_HAse_HT:",conTxt)
write.table(HAse_HT,conTxt)
HAse_HTF = HTcheck(sandwich_master_HAse_HTF,exprDat_norm,toTest3,toTest4) 
writeLines("\nsandwich_master_HAse_HTF:",conTxt)
write.table(HAse_HTF,conTxt)


###check against predicted R2 from the repeats
R2repeats(R2_norm1,CorData_norm1,"OA")
R2repeats(R2_norm2,CorData_norm2,"Injury")


###Check 7: Blood staining markers. Return correlation between bloodmarker and stains
BMarkerCor = BloodMarker(exprDat_norm)  
writeLines("\nblood marker check",conTxt)
write.table(BMarkerCor,conTxt,row.names = FALSE)

###Check 8: Sex markers
pre_norm_list = SexCheck1(exprDat_norm,metadata)
SexCheck2(pre_norm_list)

dev.off()
close(conTxt)


# ### supplement plots comparison from two MySoma directly
# ps_pca_norm1 = read.csv("ps_pca_norm1.csv")
# ps_pca_norm2 = read.csv("ps_pca_norm2.csv")
# ps_norm1 = read.csv("ps_norm1.csv")
# ps_norm2 = read.csv("ps_norm2.csv")
# pdf("MySomoPlotDirectComp.pdf")
# plotDirectComp(ps_pca_norm1,ps_pca_norm2)
# plotDirectComp(ps_norm1,ps_norm2)
# dev.off()

### summary of remover:
### protein removed: 
Remove1 = which(R2_norm1<0.1)  ### OA repeats
length(Remove1)
Remove2 = which(R2_norm2<0.1)  ### INJ repeats
length(Remove2)
### VarExp for Freeze/Thraw Pool
tempFT <- apply(calib_norm[calibIDs %in% paste0("OA POOL-HT-",c(6,25),"/29"),],2,var)/apply(exprDat_norm,2,var)
tempFT[tempFT > 1] <- 1
R2_normFT <- (1 - tempFT)^2
Remove3 = which(R2_normFT<0.1) ### Freeze/Thaw repeats
length(Remove3)
Remove4 = which(ConfounderTable2[[2]][,"SampleAge"]<(0.05/5004))  ### Anova sample age p
length(Remove4)
removeTotalPro = which(ProteinRatio>0.25)
length(removeTotalPro)
removeTotalSamp = which(SampRatio>0.25)
length(SampRatio)
removeNameP = colnames(DatSamp)[removeTotalPro]
removeName = removeNameP[!grepl("HybControlElution|Non|Spuriomer",removeNameP)]
removeIDprotein = unique(c(Remove1,Remove2,Remove3,Remove4,removeTotalPro))
length(removeIDprotein)
removeProteinName = colnames(exprDat_norm)[removeIDprotein]

###sample removed:
mismatchID50 ### gender mismatch
length(mismatchID50)
PCAout <- which(distanceToCenter>(mean(distanceToCenter) + 5*sd(distanceToCenter)))
length(PCAout)
TotalProteinOut <- which(totalProtein_norm>(mean(totalProtein_norm) +5*sd(totalProtein_norm)))
length(TotalProteinOut)
removeIDsample = unique(c(PCAout,mismatchID50,TotalProteinOut))
length(removeIDsample)
removedSampName = metadata_reord$`STEpUP Sample Identification Number (SIN)`[removeIDsample]


conTxtT = file("removeList2.csv","append")
writeLines("Summary of Removed Sample Identification Number(SIN)",conTxtT)
writeLines(removedSampName,conTxtT,sep=",")
writeLines("\nSummary of Removed Protein Name",conTxtT)
writeLines(removeProteinName,conTxtT,sep=",")

writeLines("\nRemoved protein because of R2OA<0.8",conTxtT)
writeLines(colnames(exprDat_norm)[Remove1],conTxtT,sep=",")
writeLines("\nRemoved protein because of R2INJ<0.8",conTxtT)
writeLines(colnames(exprDat_norm)[Remove2],conTxtT,sep=",")
writeLines("\nRemoved protein because of R2Freeze/Thaw<0.8",conTxtT)
writeLines(colnames(exprDat_norm)[Remove3],conTxtT,sep=",")
writeLines("\nRemoved protein because of R2Freeze/Thaw<0.8",conTxtT)
writeLines(colnames(exprDat_norm)[Remove3],conTxtT,sep=",")
writeLines("\nRemoved protein because of Anova test for plate p<1e-5",conTxtT)
writeLines(colnames(exprDat_norm)[Remove4],conTxtT,sep=",")
writeLines("\nRemoved protein because of Anova test for sample age p<1e-5",conTxtT)
writeLines(colnames(exprDat_norm)[Remove5],conTxtT,sep=",")

writeLines("\nRemoved sample because of gender mismatch",conTxtT)
writeLines(metadata_reord$`STEpUP Sample Identification Number (SIN)`[mismatchID50],conTxtT,sep=",")
writeLines("\nRemoved sample because of PCA outliers",conTxtT)
writeLines(metadata_reord$`STEpUP Sample Identification Number (SIN)`[PCAout],conTxtT,sep=",")
writeLines("\nRemoved sample because of total protein distribution",conTxtT)
writeLines(metadata_reord$`STEpUP Sample Identification Number (SIN)`[TotalProteinOut],conTxtT,sep=",")
writeLines("\nRemoved sample because of limit of detection",conTxtT)
writeLines(removeName,conTxtT,sep=",")
close(conTxtT)

write.csv(MySoma,"NormalisedSomaScanData.csv")

length(which(SampRatio>0.25))
apply(ConfounderTable2[[2]],2,function(x){length(which(x<1e-5))})

library(diptest)
dipP=vector()
for (dipcounter in 1:ncol(exprDat_norm)){
 dipP[dipcounter] = dip(exprDat_norm[,dipcounter], full.result = FALSE, min.is.0 = FALSE, debug = FALSE)
}
length(which(dipP<(0.05/5004)))

       