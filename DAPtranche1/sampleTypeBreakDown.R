## Tranche1 sample type breakdown
inputfile1 <- paste(myFilePath,"SS-200008.ADat",sep="")
RawM1 <- read.adat(inputfile1)

SampleID1 = grep("STEP",RawM1[,"SampleId"])
length(SampleID1)  ###435
SampleDat1 = RawM1[,"SampleId"][SampleID1]
length(unique(SampleDat1))  ###435

nonSampleID1 = which(!grepl("STEP",RawM1[,"SampleId"])==TRUE)
length(nonSampleID1) ### 45
nonSampleDat1 = RawM1[,"SampleId"][nonSampleID1]
levels(as.factor(nonSampleDat1))

length(which(grepl("INJ POOL-HT",RawM1[,"SampleId"])==TRUE))
# OA POOL-HT = 7
# INJ POOL-H = 7
# UNSP POOL = 0


## Tranche2 sample type breakdown
inputfile2 <- paste(myFilePath,"SS-205086.adat",sep="")
RawM2 <- read.adat(inputfile2)

SampleID2 = grep("STEP",RawM2[,"SampleId"])
length(SampleID2)  ###611
SampleDat2 = RawM2[,"SampleId"][SampleID]
length(unique(SampleDat2))  ###611
length(grep("UN",SampleDat2))
SampleDat2[grep("UN",SampleDat2)] ###18
length(grep("SP",SampleDat2))
SampleDat2[grep("SP",SampleDat2)] ###18
length(grep("V10",SampleDat2))
V10samp = SampleDat2[grep("V10",SampleDat2)]

sampleDat3 = sub("-F-V[[:alnum:]]*-HT*","",SampleDat2)### don't forget [] wildcard
sampleDat33 = sampleDat33[!grepl("UN|SP",sampleDat3)]
length(sampleDat33)
length(unique(sampleDat33)) ### 575 non duplicated sample for longitudinal

nonSampleID2 = which(!grepl("STEP",RawM2[,"SampleId"])==TRUE)
length(nonSampleID2) ### 61
nonSampleDat2 = RawM2[,"SampleId"][nonSampleID2]
levels(as.factor(nonSampleDat2))

length(which(grepl("UNSP POOL",RawM2[,"SampleId"])==TRUE))
# OA POOL-HT = 8
# INJ POOL-H = 7
# UNSP POOL = 4
unPoolID = grep("UNSP POOL",nonSampleDat2)

### check whether UNSP POOL the same as ...-F-V1-HT1-UN


RawM12SampleType = c(RawM1[,"SampleId"],RawM2[,"SampleId"])
length(RawM12SampleType) ##1152
STEP = grep("STEP",RawM12SampleType)
STEPdat = sub("-F-V[[:alnum:]]*-HT*","",RawM12SampleType[STEP])
STEPdat2 = sub("-..*","",STEPdat)
length(STEPdat2) ##1046
length(unique(STEPdat2))
### duplicated STEPIDs 
STEPdat2[duplicated(STEPdat2)] ###64
