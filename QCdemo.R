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
source("QCnormV1.R")
source("QCassess.R")


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
inputfile1 <- paste(myFilePath,"SS-200008.ADat",sep="")
inputfile2 <- paste(myFilePath,"SS-200008.hybNorm.medNormInt.plateScale.medNormRefSMP.ADat",sep="")
inputfilex <-paste(myFilePath,"SS-200008.hybNorm.medNormInt.plateScale.ADat",sep="")

initiaList<- initQCnorm(inputfile1,inputfile2) ###raw RFUs from Adat, ajust SampleType according to tranch Excel
RawM <- initiaList[[1]]
ColTable <- initiaList[[2]]
RawMList = ExtractClinicG(RawM,myFilePath) ### add columns information for our own case
RawM = RawMList[[1]]
MetaRawM = RawMList[[2]]
# SomaM <- read.adat(inputfile2) ### Normalised RFUs from Adat


### user select which normlisation methods. 
MySoma = UserNorm(Funlist,RawM)
# MySoma = RawM
# corShip = CompTWO(SomaM,MySoma) ### comparison between our normalisation and Adat
#colnames(SomaM[30:5308])[which(corShip<0.92)] ##HybControlElution highly different

#part2######################################################################################################
### initialise information from excel files
ExelDat(MySoma,MetaRawM,myFilePath)


#part3########################################################################################################
### QC assess begins. Based on different normalised MySoma.
# pdf(MySomaPlot)
# par(mfrow=c(1,3))
# conTxt = file(tableOut,"append")

##1: Total protein check. Limit of detection.
totalProtein_norm = TotalProCheck(exprDat_norm)
CompM = LoDdetection(MySoma) ###CompM > 0, RFU above limit of Detection
writeLines(paste("total RFUs below limit of Detection is ", length(which(CompM<0))," percentage of total RFUs ",length(which(CompM<0))/nrow(CompM)*ncol(CompM),"\n"),conTxt)

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
# 
# exprDat_norm_tempM <- MySoma[grep("STEP",MySoma$SampleId),]
# exprDat_norm_OA = exprDat_norm[grep("OA",exprDat_norm_tempM$SampleType),]
# listPCAOA = PCAglob(exprDat_norm_OA)
# pc_norm_OA = listPCAOA[[1]]
# listPCAOA[[2]]
# listPCAOA[[3]]
# exprDat_norm_INJ = exprDat_norm[grep("INJ",exprDat_norm_tempM$SampleType),]
# listPCAINJ = PCAglob(exprDat_norm_INJ)
# pc_norm_INJ = listPCAINJ[[1]]
# listPCAINJ[[2]]
# listPCAINJ[[3]]

###Check 4: Techical confounders
ConfounderTable1 = ConfouderCheck(totalProtein_norm,pc_norm$x[,1:5]) ### confounders against PCs 
ps_pca_norm = ConfounderTable1[[2]]
writeLines("\ntechnical counfoundr check in order: total protein,ps_pca_norm,varExp_norm",conTxt)
for (ct in 1:3) {write.table(ConfounderTable1[[ct]], conTxt)}
write.csv(ps_pca_norm, file = ps_pca_normOut) ### restore a ps_norm formal file, for direct comparison within one graph between two MySoma

ConfounderTable2 = ConfouderCheck(totalProtein_norm,exprDat_norm) ### confounders against each proteins 
ps_norm = ConfounderTable2[[2]]
writeLines("\ntechnical counfoundr check in order: ps_norm,varExp_norm",conTxt)
for (ct in 2:3) {write.table(ConfounderTable2[[ct]], conTxt)} 
write.csv(ps_norm, file = ps_normOut) 

confounderPlot1(ConfounderTable1)  ###confounder vs total protein
confounderPlot2(ConfounderTable1,0) ### confounder vs per PC
confounderPlot2(ConfounderTable2,1) ### confounder vs per protein

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
pre_norm_list = SexCheck1(exprDat_norm,metadata_reord)
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
Remove1 = which(R2_norm1<0.8)  ### OA repeats
Remove2 = which(R2_norm2<0.8)  ### INJ repeats
### VarExp for Freeze/Thraw Pool
tempFT <- apply(calib_norm[calibIDs %in% paste0("OA POOL-HT-",c(6,25),"/29"),],2,var)/apply(exprDat_norm,2,var)
tempFT[tempFT > 1] <- 1
R2_normFT <- (1 - tempFT)^2
Remove3 = which(R2_normFT<0.8) ### Freeze/Thaw repeats
Remove4 = which(temp1<1e-5)  ### Anova Plate p
Remove5 = which(temp2<1e-5)  ### Anova sample age p
length(unique(c(Remove1,Remove2,Remove3,Remove4,Remove5)))
removeIDprotein = unique(c(Remove1,Remove2,Remove3,Remove4,Remove5))
removeProteinName = colnames(exprDat_norm)[removeIDprotein]
###sample removed:
mismatchID50 ### gender mismatch
PCAout <- which(distanceToCenter>(mean(distanceToCenter) + 2*sd(distanceToCenter)))
TotalProteinOut <- which(totalProtein_norm>(mean(totalProtein_norm) +2*sd(totalProtein_norm)))
removeIDsample = unique(c(PCAout,mismatchID50,TotalProteinOut))
removedSampName = metadata_reord$`STEpUP Sample Identification Number (SIN)`[removeIDsample]

removeTotalPro = which(ProteinRatio>0.25)
removeNameP = colnames(DatSamp)[removeTotalPro]
removeName = removeNameP[!grepl("HybControlElution|Non|Spuriomer",removeNameP)]

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

