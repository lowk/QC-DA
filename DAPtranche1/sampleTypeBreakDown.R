## Tranche1 sample type breakdown
inputfile1 <- paste(myFilePath,"SS-200008.ADat",sep="")
RawM1 <- read.adat(inputfile1)

nonSampleID1 = which(!grepl("STEP",RawM1[,"SampleId"])==TRUE)
length(nonSampleID1) ### 45
nonSampleDat1 = RawM1[,"SampleId"][nonSampleID1]
levels(as.factor(nonSampleDat1))

length(which(grepl("INJ POOL-HT",RawM1[,"SampleId"])==TRUE))
# OA POOL-HT = 7
# INJ POOL-H = 7
# UNSP POOL = 0

### longitudinal sample check
SampleID1 = grep("STEP",RawM1[,"SampleId"])
length(SampleID1)  ###435
SampleDat1 = RawM1[,"SampleId"][SampleID1]
SampleDat11 = sub("-F-V[[:alnum:]]*-HT*","",SampleDat1)
length(which(duplicated(SampleDat11)==TRUE))  ###no duplicated samples

## Tranche2 sample type breakdown
clinic <- read_excel("/Users/ydeng/Documents/QCstepOA/clinic/STepUP_QCData_Tranche2_09MAR2021.xlsx", col_names = TRUE)
clinicSTEP = clinic$...1

inputfile2 <- paste(myFilePath,"SS-205086.adat",sep="")
RawM2 <- read.adat(inputfile2)

nonSampleID2 = which(!grepl("STEP",RawM2[,"SampleId"])==TRUE)
length(nonSampleID2) ### 61
nonSampleDat2 = RawM2[,"SampleId"][nonSampleID2]
levels(as.factor(nonSampleDat2))

length(which(grepl("UNSP POOL",RawM2[,"SampleId"])==TRUE))
# OA POOL-HT = 8
# INJ POOL-H = 7
# UNSP POOL = 4
unPoolID = grep("UNSP POOL",nonSampleDat2)

###longitudinal sample check
SampleID2 = grep("STEP",RawM2[,"SampleId"])
length(SampleID2)  ###611
SampleDat2 = RawM2[,"SampleId"][SampleID2]

### remove spuns
SampleDat22 =SampleDat2[!grepl("SP",SampleDat2)==TRUE]
length(SampleDat22)  ###593
length(which(clinic$...3=="OA"))###532
length(which(clinic$...3=="Injury"))###72
length(which(clinic$...3=="RA"))###7
length(which(clinic$...3=="Control"))###3

length(grep("UN",SampleDat2))
SampleDat2[grep("UN",SampleDat2)] ###18
length(grep("SP",SampleDat2))
SampleDat2[grep("SP",SampleDat2)] ###18

SampleDat222 = sub("-F-V.","",SampleDat22)
SampleDat222 = sub("-HT.","",SampleDat222)
SampleDat222 = sub("-UN","",SampleDat222)

###within tranche2, duplicated samples
length(which(duplicated(SampleDat222)==TRUE)) ###32
tranche2DupSTEP = which(duplicated(SampleDat222)==TRUE)
tranche2DupDat = SampleDat222[tranche2DupSTEP]

#IDclinicT2 = vector(mode="list",length=length(tranche2DupDat))
IDclinicT2 = vector(mode="integer",length=length(tranche2DupDat))
for (longCounter in 1:length(tranche2DupDat)){
  IDclinicT2[longCounter] = grep(tranche2DupDat[longCounter],clinicSTEP)[1]
} ### 3 injury samples measured 3 times.

longType = clinic$...3[IDclinicT2]
###longitudinal OA
length(grep("OA",longType))  ###18
###longitudinal Injury
length(grep("Injury",longType) )###14
### if tranche2 only, baseline OA: 532-18 =514
### if tranche2 only, baseline injury: 72-14 -3 =55


###combine two tranches data
RawM12SampleType = c(RawM1[,"SampleId"],RawM2[,"SampleId"])
length(RawM12SampleType) ##1152
STEP = grep("STEP",RawM12SampleType)
STEPdat = RawM12SampleType[STEP]
#STEPdat = sub("-F-V[:alnum:]-HT*","",RawM12SampleType[STEP])
STEPdat2 = sub("-..*","",STEPdat)
length(STEPdat2) ##1046
length(unique(STEPdat2))

### duplicated STEPIDs 
longSamp = unique(STEPdat2[duplicated(STEPdat2)]) 

duplicatedTimes=vector()
for (ii in 1:length(longSamp)){
  duplicatedTimes[ii] = length(grep(longSamp[ii],STEPdat))
} 
length(which(duplicatedTimes==3))  ### 6 samples measured 3 times
duplicated3ID = which(duplicatedTimes==3)
duplicated3STEP = longSamp[duplicated3ID] ### sample name measured for 3 times 

IDclinic = vector(mode="integer",length=length(longSamp))
for (longCounter in 1:length(longSamp)){
  IDclinic[longCounter] = grep(longSamp[longCounter],clinicSTEP)[1]
}


longType = clinic$...3[IDclinic]
###longitudinal OA
length(grep("OA",longType))  ###36
###longitudinal Injury
length(grep("Injury",longType) )###25

### extract sample type for duplicated sample of 3 times
for (longCounter in 1:length(duplicated3STEP)){
  IDclinic[longCounter] = grep(duplicated3STEP[longCounter],clinicSTEP)[1]
}
longType = clinic$...3[IDclinic]
length(grep("OA",longType))  ###0
###longitudinal Injury
length(grep("Injury",longType) )###3

### if tranche1 and 2 together, baseline OA 258 + 532 -36 = 754
### if tranche1 and 2 together, baseline Injury 175 + 72 -25 -3*2  = 216

