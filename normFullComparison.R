### Comparisons among different normalisations
library(readxl) 
library(SomaDataIO)
source("QCnormTranche2.R")
library(sva)
library(umap)
library(GGally)

myFilePath <- "/Users/ydeng/Documents/QCstepOA/normComp/"  ### set required file directory
inputfile1 <- paste(myFilePath,"SS-200008.ADat",sep="") ### RawM1
inputfile2 <- paste(myFilePath,"SS-205086.adat",sep="") ###RawM2
inputfile12 <- paste(myFilePath,"SS-200008.hybNorm.medNormInt.plateScale.medNormRefSMP.ADat",sep="")
inputfile3 <- paste(myFilePath,"SS-200008_v4_SynovialFluid.hybNormRef.medNormInt.plateScale.calibrate.20211007.adat",sep="") ### Soma1
inputfile4 <- paste(myFilePath,"SS-205086_v4_SynovialFluid.hybNormRef.medNormInt.plateScale.calibrate.20211007.adat",sep="")  ###Soma2
pdf(file = "/Users/ydeng/Documents/QCstepOA/normComp/normPlot.pdf")


### retrieve soma data
RawM1 = read.adat(inputfile1)
RawM2 = read.adat(inputfile2)
Soma1 = read.adat(inputfile3)
Soma2 = read.adat(inputfile4)
ColTable = initQCnorm(inputfile1,inputfile12) ### change to use attr(read.adat(inputfile12),"Col.Meta") as Luke suggested

### full normalisation
Funlist1 = list(HYBNORM,MIDNORMcali,PLATESCALE,CALIBRATION)
MySoma1Full = UserNorm(Funlist1,RawM1)
MySoma2Full = UserNorm(Funlist1,RawM2)
MySoma1FullLog = UserNorm(Funlist1,log(RawM1)) ### log transformed RFU before normalisation steps
MySoma2FullLog = UserNorm(Funlist1,log(RawM2))

### hybnorm only normalisation
Funlist2 = list(HYBNORM)
MySoma1Hyb = UserNorm(Funlist2,RawM1)
MySoma2Hyb = UserNorm(Funlist2,RawM2)

### compare our normalisation and soma new normalisation (2021.10.08)
CompTWO(Soma1,MySoma1Full,0.9)
CompTWO(Soma2,MySoma2Full,0.9)
CompTWO(Soma1,MySoma1Hyb,0.7)
CompTWO(Soma2,MySoma2Hyb,0.7)

# SomaId of new normliased version (2021.10.08) are not consistent with previous version(2020.07.08) either.
# SomaId1 = attributes(Soma1)$Col.Meta["SomaId"]
# SomaId2 = attributes(RawM1)$Col.Meta["SomaId"]
# identical(SomaId1,SomaId2)
# changeOrder = vector()
# SomaIDKK=1
# notMatchSomaID=vector()
# for (k in 1:nrow(SomaId2)){
#   if(!any(which(as.matrix(SomaId1[,1]) %in% SomaId2[k,1]))){notMatchSomaID[SomaIDKK] = k
#   SomaIDKK=SomaIDKK+1}
#   else{changeOrder[k] = which(as.matrix(SomaId1[,1]) %in% SomaId2[k,1])}
# }
# SomaId2[notMatchSomaID,1]

### combat on combined two tranches -> PCA|UMAP regarding to batches -> KNN test
# Codes in order for Types of normalization:
# Raw data
# Hybnorm only
# Hybnorm only + ComBat
# “Full norm”, i.e. using tranche 1 calibrators for both tranches
# Full norm + ComBat
# Darryl norm (i.e. SomaLogic’s new normalization)
# Darryl norm + ComBat

# Raw data
TestRaw <- MyCombat(RawM1,RawM2,1)
CombinedRaw = TestRaw[[1]] ###combined expression data
TestRaw[[2]] ###KNN test results
batchMeta_Raw = TestRaw[[3]]

PlotPCA(CombinedRaw,batchMeta_Raw,3,1,"PCA on combined raw data")
PlotPCA(CombinedRaw,batchMeta_Raw,3,2,"PCA on combined raw data")
PlotUmap(CombinedRaw,batchMeta_Raw,1,"UMAP on combined raw data")
PlotUmap(CombinedRaw,batchMeta_Raw,2,"UMAP on combined raw data")

CVbreak(RawM1,RawM2,"OA",CombinedRaw,"combined raw data")
CVbreak(RawM1,RawM2,"INJ",CombinedRaw,"combined raw data")
print("RawM1 accuracy against external controls")
ExtVal(RawM1)

# Hybnorm only
TestHyb <- MyCombat(MySoma1Hyb,MySoma2Hyb,1) ### without combat
CombinedHyb = TestHyb[[1]] ###combined expression data
TestHyb[[2]]###KNN test results
batchMeta_HybOnly = TestHyb[[3]]

PlotPCA(CombinedHyb,batchMeta_HybOnly,3,1,"PCA on combined hyb only data")
PlotPCA(CombinedHyb,batchMeta_HybOnly,3,2,"PCA on combined hyb only data")
PlotUmap(CombinedHyb,batchMeta_HybOnly,1,"UMAP on combined hyb only data")
PlotUmap(CombinedHyb,batchMeta_HybOnly,2,"UMAP on combined hyb only data")

CVbreak(MySoma1Hyb,MySoma2Hyb,"OA",CombinedHyb,"combined hyb only data")
CVbreak(MySoma1Hyb,MySoma2Hyb,"INJ",CombinedHyb,"combined hyb only data")

print("MySoma1Hyb accuracy against external controls")
ExtVal(MySoma1Hyb)


# Hybnorm only + ComBat
TestMyHyb <- MyCombat(MySoma1Hyb,MySoma2Hyb,0)
combat_MySomaHyb <- TestMyHyb[[1]] ###combined expression data
TestMyHyb[[2]]###KNN test results
batchMeta_MySomaHyb = TestMyHyb[[3]]

PlotPCA(combat_MySomaHyb,batchMeta_MySomaHyb,3,1,"PCA on combat_MySomaHyb")
PlotPCA(combat_MySomaHyb,batchMeta_MySomaHyb,3,2,"PCA on combat_MySomaHyb")
PlotUmap(combat_MySomaHyb,batchMeta_MySomaHyb,1,"UMAP on combat_MySomaHyb")
PlotUmap(combat_MySomaHyb,batchMeta_MySomaHyb,2,"UMAP on combat_MySomaHyb")

CVbreak(MySoma1Hyb,MySoma2Hyb,"OA",combat_MySomaHyb,"combat_MySomaHyb")
CVbreak(MySoma1Hyb,MySoma2Hyb,"INJ",combat_MySomaHyb,"combat_MySomaHyb")

print("combat_MySomaHyb accuracy against external controls")
cutMySoma1Hyb = MySoma1Hyb[which(grepl("Sample",MySoma1Hyb[,"SampleType"])),1:(which(colnames(MySoma1Hyb)=="CRYBB2.10000.28")-1)]
ExtVal(cbind(cutMySoma1Hyb,combat_MySomaHyb[1:nrow(cutMySoma1Hyb),]))

# “Full norm”, i.e. using tranche 1 calibrators for both tranches
TestMyFullOnly <- MyCombat(MySoma1Full,MySoma2Full,1)
MyFullOnly <- TestMyFullOnly[[1]] ###combined expression data
TestMyFullOnly[[2]]###KNN test results
batchMeta_MyFullOnly = TestMyFullOnly[[3]]

PlotPCA(MyFullOnly,batchMeta_MyFullOnly,3,1,"PCA on combined my full norm")  ###PCA and UMAP plot
PlotPCA(MyFullOnly,batchMeta_MyFullOnly,3,2,"PCA on combined my full norm")
PlotUmap(MyFullOnly,batchMeta_MyFullOnly,1,"UMAP on combined my full norm")
PlotUmap(MyFullOnly,batchMeta_MyFullOnly,2,"UMAP on combined my full norm")

CVbreak(MySoma1Full,MySoma2Full,"OA",log(MyFullOnly),"combined my full norm")
CVbreak(MySoma1Full,MySoma2Full,"INJ",log(MyFullOnly),"combined my full norm")

print("MySoma1Full accuracy against external controls")
ExtVal(MySoma1Full)

# Full norm + ComBat
TestMyFull <- MyCombat(MySoma1Full,MySoma2Full,0)
combat_MySomaFull <- TestMyFull[[1]] ###combined expression data
TestMyFull[[2]]###KNN test results
batchMeta_MySomaFul = TestMyFull[[3]]

PlotPCA(combat_MySomaFull,batchMeta_MySomaFul,3,1,"PCA on combat_MySomaFul")
PlotPCA(combat_MySomaFull,batchMeta_MySomaFul,3,2,"PCA on combat_MySomaFul")
PlotUmap(combat_MySomaFull,batchMeta_MySomaFul,1,"UMAP on combat_MySomaFul")
PlotUmap(combat_MySomaFull,batchMeta_MySomaFul,2,"UMAP on combat_MySomaFul")

CVbreak(MySoma1Full,MySoma2Full,"OA",combat_MySomaFull,"combat_MySomaFul")
CVbreak(MySoma1Full,MySoma2Full,"INJ",combat_MySomaFull,"combat_MySomaFul")

print("combat_MySomaFul accuracy against external controls")
cutMySoma1Full = MySoma1Full[which(grepl("Sample",MySoma1Full[,"SampleType"])),1:(which(colnames(MySoma1Full)=="CRYBB2.10000.28")-1)]
ExtVal(cbind(cutMySoma1Full,combat_MySomaFull[1:nrow(cutMySoma1Hyb),]))


# Darryl norm (i.e. SomaLogic’s new normalization)
TestSomaOnly <- MyCombat(Soma1,Soma2,1)
SomaOnly = TestSomaOnly[[1]] ###combined expression data
TestSomaOnly[[2]]###KNN test results
batchMeta_SomaOnly = TestSomaOnly[[3]]

PlotPCA(SomaOnly,batchMeta_SomaOnly,3,1,"PCA on combined Soma Only")
PlotPCA(SomaOnly,batchMeta_SomaOnly,3,2,"PCA on combined Soma Only")
PlotUmap(SomaOnly,batchMeta_SomaOnly,1,"UMAP on combined Soma Only")
PlotUmap(SomaOnly,batchMeta_SomaOnly,2,"UMAP on combined Soma Only")

CVbreak(Soma1,Soma2,"OA",SomaOnly,"combined Soma Only")
CVbreak(Soma1,Soma2,"INJ",SomaOnly,"combined Soma Only")

print("SomaOnly accuracy against external controls")
ExtVal(Soma1)


# Darryl norm + ComBat
TestSoma <- MyCombat(Soma1,Soma2,0)
combat_Soma = TestSoma[[1]] ###combined expression data
TestSoma[[2]]###KNN test results
batchMeta_Soma = TestSoma[[3]]

PlotPCA(combat_Soma,batchMeta_Soma,3,1,"PCA on combat_Soma")
PlotPCA(combat_Soma,batchMeta_Soma,3,2,"PCA on combat_Soma")
PlotUmap(combat_Soma,batchMeta_Soma,1,"UMAP on combat_Soma")
PlotUmap(combat_Soma,batchMeta_Soma,2,"UMAP on combat_Soma")

CVbreak(Soma1,Soma2,"OA",combat_Soma,"combat_Soma")
CVbreak(Soma1,Soma2,"INJ",combat_Soma,"combat_Soma")

print("combat_Soma accuracy against external controls")
cutCombat_Soma = Soma1[which(grepl("Sample",Soma1[,"SampleType"])),1:(which(colnames(Soma1)=="CRYBB2.10000.28")-1)]
ExtVal(cbind(cutCombat_Soma,combat_Soma[1:nrow(cutCombat_Soma),]))

# log + “Full norm”, i.e. using tranche 1 calibrators for both tranches
TestMyFullLogOnly <- MyCombat(MySoma1FullLog,MySoma2FullLog,1)
MyFullLogOnly <- TestMyFullLogOnly[[1]] ###combined expression data
TestMyFullLogOnly[[2]]###KNN test results
batchMeta_MyFullLogOnly = TestMyFullLogOnly[[3]]

PlotPCA(MyFullLogOnly,batchMeta_MyFullLogOnly,3,1,"PCA on combined my full norm on log RFU")  ###PCA and UMAP plot
PlotPCA(MyFullLogOnly,batchMeta_MyFullLogOnly,3,2,"PCA on combined my full norm on log RFU")
PlotUmap(MyFullLogOnly,batchMeta_MyFullLogOnly,1,"UMAP on combined my full norm on log RFU")
PlotUmap(MyFullLogOnly,batchMeta_MyFullLogOnly,2,"UMAP on combined my full norm on log RFU")

CVbreak(MySoma1FullLog,MySoma2FullLog,"OA",MyFullLogOnly,"combined my full norm on log RFU")
CVbreak(MySoma1FullLog,MySoma2FullLog,"INJ",MyFullLogOnly,"combined my full norm on log RFU")

print("MySoma1FullLog accuracy against external controls")
ExtVal(MySoma1FullLog)

# log + Full norm + ComBat
TestFullLogCombat <- MyCombat(MySoma1FullLog,MySoma2FullLog,0)
combat_MySomaFullLog <- TestFullLogCombat[[1]] ###combined expression data
TestFullLogCombat[[2]]###KNN test results
batchMeta_MySomaFullLog = TestFullLogCombat[[3]]

PlotPCA(combat_MySomaFullLog,batchMeta_MySomaFullLog,3,1,"PCA on combatted my full norm on log RFU")  ###PCA and UMAP plot
PlotPCA(combat_MySomaFullLog,batchMeta_MySomaFullLog,3,2,"PCA on combatted my full norm on log RFU")
PlotUmap(combat_MySomaFullLog,batchMeta_MySomaFullLog,1,"UMAP on combatted my full norm on log RFU")
PlotUmap(combat_MySomaFullLog,batchMeta_MySomaFullLog,2,"UMAP on combatted my full norm on log RFU")

CVbreak(MySoma1FullLog,MySoma2FullLog,"OA",combat_MySomaFullLog,"combatted my full norm on log RFU")
CVbreak(MySoma1FullLog,MySoma2FullLog,"INJ",combat_MySomaFullLog,"combatted my full norm on log RFU")

print("Combatted MySoma1FullLog accuracy against external controls")
cutCombat_MyFullCombat = MySoma1FullLog[which(grepl("Sample",MySoma1FullLog[,"SampleType"])),1:(which(colnames(MySoma1FullLog)=="CRYBB2.10000.28")-1)]
ExtVal(cbind(cutCombat_MyFullCombat,combat_MySomaFullLog[1:nrow(cutCombat_MyFullCombat),]))


dev.off()

# Raw data: RawM1,RawM2, CombinedRaw is the combined raw data.
# Hybnorm only: HubNorm only on seperate tranche individually: MySoma1Hyb, MySoma2Hyb; CombinedHyb is the expression profile of combined hyb normalised data.
# Hybnorm only + ComBat: combat_MySomaHyb is the expression profile by further combatted  on CombinedHyb.
# “Full norm”: full norm on seprate tranche individually: MySoma1Full,MySoma2Full; MyFullOnly is the expression profile of combined full normalised data.
# Full norm + ComBat: combat_MySomaFull is the expression profile by further combatted on MyFullOnly.
# Darryl norm: Daryyl normliased: Soma1,Soma2; SomaOnly is the expression profile of combined Soma1 and Soma2.
# Darryl norm + ComBat: combat_Soma is the expression profile by further combatted on SomaOnly.
# note that: combined data sets and combatted data sets only include human samples and human proteins.

normList = list("RawM1"=RawM1,"RawM2"=RawM2,"CombinedRaw"=CombinedRaw,
                "MySoma1Hyb"=MySoma1Hyb,"MySoma2Hyb"=MySoma2Hyb,"CombinedHyb"=CombinedHyb,"combat_MySomaHyb"=combat_MySomaHyb,
                "MySoma1Full"=MySoma1Full,"MySoma2Full"=MySoma2Full,"MyFullOnly"=MyFullOnly,"combat_MySomaFull"=combat_MySomaFull,
                "Soma1"=Soma1,"Soma2"=Soma2,"SomaOnly"=SomaOnly,"combat_Soma"=combat_Soma)
saveRDS(normList,file="normList.rds")

batchMeta = list()batchMeta_Raw

batchMeta_HybOnly

batchMeta_MySomaHyb

batchMeta_SomaOnly

batchMeta_Soma
