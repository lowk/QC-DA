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

### check my code accuracy by comparing our normalisation and somalogic provided normalised data
checkMyCode <- function(RawM,SomaFile,Flist){
  SomaM <- read.adat(SomaFile)
  MySoma = UserNorm(Flist,RawM)
  corShip = CompTWO(SomaM,MySoma) 
  DatStartId <- which(colnames(MySoma)=="CRYBB2.10000.28")
  colnames(MySoma[DatStartId:ncol(MySoma)])[which(corShip<0.7)] ### different analyses
}


### RFUs Normalisation, across all the plates
getMySoma <- function(Platelist){
  MySoma = data.frame(Platelist[[1]]) ###MySoma: the whole dataframe for all the plates
  for (plateC in 2:length(Platelist)){
    MySomaT = data.frame(Platelist[[plateC]])
    MySoma = rbind(MySoma,MySomaT)
  }
  return(MySoma)
}
### input raw adat file. undergo selected normalisation. output validated RFUs for further analysis

### 1. Hybridization normalization performed on one plate. Input RawM, output RawM

HYBNORM <- function(RawM){
  
  PlateIdUni = levels(factor(RawM$PlateId))
  
  Platelist = list()
  
  DatStartId <- which(colnames(RawM)=="CRYBB2.10000.28") ###for calculation convenience, extract data zone only. Only tranche 2 has "RMA" column
  
  HybId = which(grepl("HybControlElution",colnames(RawM))) -(DatStartId-1)
  
  for (plateCounter in 1:length(PlateIdUni)){
    
    PlateIdSg = which(RawM$PlateId == PlateIdUni[plateCounter])
    
    RawMS = RawM[PlateIdSg,] ### single plate
    
    datZone = RawMS[,DatStartId:ncol(RawMS)] ###RFUs zone only
    
    rCmedian = matrix(apply(datZone[,HybId],2,median),nrow=1)
    
    HybNorm1 = t(apply(datZone[,HybId],1,function(x){rCmedian/x}))
    
    rRmedian = apply(HybNorm1,1,median)
    
    HybNorm = apply(datZone,2,function(x){x*rRmedian})
    
    RawMSDone = cbind(RawMS[,1:(DatStartId-1)],HybNorm)
    
    Platelist[[plateCounter]] = RawMSDone 
  }
  
  MySoma = getMySoma(Platelist)  
  
  return(MySoma)
}



### 2.Plate scaling

PLATESCALE <- function(RawM){
  
  PlateIdUni = levels(factor(RawM$PlateId))
  
  Platelist = list()
  
  DatStartId <- which(colnames(RawM)=="CRYBB2.10000.28")  ###for calculation convenience, extract data zone only
  
  for (plateCounter in 1:length(PlateIdUni)){
    
    PlateIdSg = which(RawM$PlateId == PlateIdUni[plateCounter])
    
    RawMS = RawM[PlateIdSg,] ### single plate
    
    idCaliborator = which(RawMS$SampleType=="Calibrator")
    
    datZone = RawMS[,DatStartId:ncol(RawMS)]
    
    CaliboratorM = datZone[idCaliborator,]
    
    calibratorMedian = apply(CaliboratorM,2,median)
    
    PlateScaleRatio = PlateScale_Reference/calibratorMedian
    
    PlateScaleScalar = median(PlateScaleRatio)
    
    datZone2 = datZone*PlateScaleScalar
    
    RawMSDone = cbind(RawMS[,1:(DatStartId-1)],datZone2)
    
    Platelist[[plateCounter]] = RawMSDone 
  }
  
  MySoma = getMySoma(Platelist)  
  
  return(MySoma)
}


### 3. Calibration

CALIBRATION <- function(RawM){
  
  PlateIdUni = levels(factor(RawM$PlateId))
  
  Platelist = list()
  
  DatStartId <- which(colnames(RawM)=="CRYBB2.10000.28")  ###for calculation convenience, extract data zone only
  
  for (plateCounter in 1:length(PlateIdUni)){
    
    PlateIdSg = which(RawM$PlateId == PlateIdUni[plateCounter])
    
    RawMS = RawM[PlateIdSg,] ### single plate
    
    idCaliborator = which(RawMS$SampleType=="Calibrator")
    
    datZone = RawMS[,DatStartId:ncol(RawMS)]
    
    CaliboratorM = datZone[idCaliborator,]
    
    CalibratorMedian = apply(CaliboratorM,2,median)
    
    CalSet=PlateScale_Reference/CalibratorMedian
    
    datZone2 = t(apply(datZone,1,function(x) x*CalSet))
    
    RawMSDone = cbind(RawMS[,1:DatStartId-1],datZone2)
    
    Platelist[[plateCounter]] = RawMSDone 
  }
  
  MySoma = getMySoma(Platelist)  
  
  return(MySoma)
}


MIDNORM = function(RawM){ ###caliCase control which sample type to be applied MidNorm
  
  PlateIdUni = levels(factor(RawM$PlateId))
  
  Platelist = list()
  
  for (plateCounter in 1:length(PlateIdUni)){
    
    PlateIdSg = which(RawM$PlateId == PlateIdUni[plateCounter])
    
    RawMS = RawM[PlateIdSg,] ### single plate
    
    idHyb = which(grepl("HybControlElution",colnames(RawMS))==TRUE)
    idNonHyb = which(!grepl("HybControlElution",colnames(RawMS))==TRUE)
    RawMS1 = RawMS[,idNonHyb] ###RawM1 to track "Ratio of Normalization Median to Sample Value"
    DatStartId <- which(colnames(RawMS1)=="CRYBB2.10000.28")
    DatStartIdP <- which(colnames(RawMS)=="CRYBB2.10000.28")
    
    sampType = levels(factor(RawMS$SampleType))
    
    for (sampTypeCounter in 1:length(sampType)){
      
      idSamp = which(RawMS1$SampleType == sampType[[sampTypeCounter]])
      
      datZoneP = RawMS1[,DatStartId:ncol(RawMS1)] ### datazone excludes "HybControlElution"
      
      datZone = RawMS1[idSamp,DatStartId:ncol(RawMS1)]  ###dataZone only include interested sample type
      
      if(length(idSamp)==1) {SampTypeRFU = matrix(datZone,nrow=1)
      SampTypeMedian = SampTypeRFU
      SampMedianNorm = as.matrix(SampTypeMedian/SampTypeRFU)}
      
      else {SampTypeRFU = datZone
      SampTypeMedian = apply(SampTypeRFU,2,median)
      SampMedianNorm = t(apply(SampTypeRFU,1,function(x){SampTypeMedian/x}))}
      
      Dilute = Dilution[which(Dilution!="0")] ###SampMedianNorm & Dilute: ncol(SampMedianNorm)=length(Dilute)
      uniqDilute = levels(factor(Dilute))
      
      for (idDilute in (1:length(uniqDilute))){
        
        DataDiluteID = which(Dilute==uniqDilute[idDilute])
        
        MedianDiluteSingle = matrix(apply(SampMedianNorm[,DataDiluteID],1,median),ncol=1)
        
        DataDiluteNormT = apply(datZone[,DataDiluteID],2,function(x){MedianDiluteSingle*x})
        
        if (idDilute==1){DataDiluteNorm = DataDiluteNormT}
        else{DataDiluteNorm = cbind(DataDiluteNorm,DataDiluteNormT)}
      }
      
      if (sampTypeCounter==1) {RawMStemp = cbind(RawMS[idSamp,1:(DatStartIdP-1)],RawMS[idSamp,idHyb],DataDiluteNorm)}
      else{
        RawMStempT = cbind(RawMS[idSamp,1:(DatStartIdP-1)],RawMS[idSamp,idHyb],DataDiluteNorm)
        RawMStemp = rbind(RawMStemp,RawMStempT)
      }
    }
    
    Platelist[[plateCounter]] = RawMStemp
  }
  
  MySomaTemp = getMySoma(Platelist)
  
  ###up to date MySomaTemp has different indexing from Raw. We make them consistent
  rowOrderName = rownames(RawM)
  rowOrder = vector(mode="numeric",length=nrow(MySomaTemp))
  for (j in 1:nrow(MySomaTemp)){
    rowOrder[j] = which(rownames(MySomaTemp) %in% rowOrderName[j])
  }
  
  colOrderName = colnames(RawM)
  colOrder = vector(mode="numeric",length=ncol(MySomaTemp))
  for (k in 1:ncol(MySomaTemp)){
    colOrder[k] = which(colnames(MySomaTemp) %in% colOrderName[k])
  }
  
  MySoma = MySomaTemp[rowOrder,colOrder]
  return(MySoma)
}

MIDNORMcali = function(RawM){ ###caliCase control which sample type to be applied MidNorm
  
  PlateIdUni = levels(factor(RawM$PlateId))
  
  Platelist = list()
  
  for (plateCounter in 1:length(PlateIdUni)){
    
    PlateIdSg = which(RawM$PlateId == PlateIdUni[plateCounter])
    
    RawMS = RawM[PlateIdSg,] ### single plate
    
    idHyb = which(grepl("HybControlElution",colnames(RawMS))==TRUE)
    idNonHyb = which(!grepl("HybControlElution",colnames(RawMS))==TRUE)
    RawMS1 = RawMS[,idNonHyb] ###RawM1 to track "Ratio of Normalization Median to Sample Value"
    DatStartId = which(colnames(RawMS1) == "CRYBB2.10000.28")
    DatStartIdP = which(colnames(RawMS) == "CRYBB2.10000.28")
    
    sampType = "Calibrator"
    
    for (sampTypeCounter in 1:length(sampType)){
      
      idSamp = which(RawMS1$SampleType == sampType[[sampTypeCounter]])
      
      datZoneP = RawMS1[,DatStartId:ncol(RawMS1)] ### datazone excludes "HybControlElution"
      
      datZone = RawMS1[idSamp,DatStartId:ncol(RawMS1)]  ###dataZone only include interested sample type
      
      if(length(idSamp)==1) {SampTypeRFU = matrix(datZone,nrow=1)
      SampTypeMedian = SampTypeRFU
      SampMedianNorm = as.matrix(SampTypeMedian/SampTypeRFU)
      
      }else{SampTypeRFU = datZone
      SampTypeMedian = apply(SampTypeRFU,2,median)
      SampMedianNorm = t(apply(SampTypeRFU,1,function(x){SampTypeMedian/x}))}
      
      Dilute = Dilution[which(Dilution!="0")] ###SampMedianNorm & Dilute: ncol(SampMedianNorm)=length(Dilute)
      uniqDilute = levels(factor(Dilute))
      
      for (idDilute in (1:length(uniqDilute))){
        
        DataDiluteID = which(Dilute==uniqDilute[idDilute])
        
        MedianDiluteSingle = matrix(apply(SampMedianNorm[,DataDiluteID],1,median),ncol=1)
        
        DataDiluteNormT = apply(datZone[,DataDiluteID],2,function(x){MedianDiluteSingle*x})
        
        if (idDilute==1){DataDiluteNorm = DataDiluteNormT
        }else{DataDiluteNorm = cbind(DataDiluteNorm,DataDiluteNormT)}
      }
      
      if (sampTypeCounter==1) {RawMStemp = cbind(RawMS[idSamp,1:(DatStartIdP-1)],RawMS[idSamp,idHyb],DataDiluteNorm)}
      else{
        RawMStempT = cbind(RawMS[idSamp,1:(DatStartIdP-1)],RawMS[idSamp,idHyb],DataDiluteNorm)
        RawMStemp = rbind(RawMStemp,RawMStempT)
      }
      RawMStemp2 = rbind(RawMStemp,RawMS[which(RawMS$SampleType!="Calibrator"),])
    }
    
    Platelist[[plateCounter]] = RawMStemp2
  }
  
  MySomaTemp = getMySoma(Platelist)
  
  ###up to date MySomaTemp has different indexing from Raw. We make them consistent
  rowOrderName = rownames(RawM)
  rowOrder = vector(mode="numeric",length=nrow(MySomaTemp))
  for (j in 1:nrow(MySomaTemp)){
    rowOrder[j] = which(rownames(MySomaTemp) %in% rowOrderName[j])
  }
  
  colOrderName = colnames(RawM)
  colOrder = vector(mode="numeric",length=ncol(MySomaTemp))
  for (k in 1:ncol(MySomaTemp)){
    colOrder[k] = which(colnames(MySomaTemp) %in% colOrderName[k])
  }
  
  MySoma = MySomaTemp[rowOrder,colOrder]
  return(MySoma)
}

MIDNORMsamp = function(RawM){ 
  
  PlateIdUni = levels(factor(RawM$PlateId))
  
  Platelist = list()
  
  for (plateCounter in 1:length(PlateIdUni)){
    
    PlateIdSg = which(RawM$PlateId == PlateIdUni[plateCounter])
    
    RawMS = RawM[PlateIdSg,] ### single plate
    
    idHyb = which(grepl("HybControlElution",colnames(RawMS))==TRUE)
    idNonHyb = which(!grepl("HybControlElution",colnames(RawMS))==TRUE)
    RawMS1 = RawMS[,idNonHyb] ###RawM1 to track "Ratio of Normalization Median to Sample Value"
    DatStartId = which(colnames(RawMS1) == "CRYBB2.10000.28")
    DatStartIdP = which(colnames(RawMS) == "CRYBB2.10000.28")
    
    sampType = c("Buffer","Sample")
    
    for (sampTypeCounter in 1:length(sampType)){
      
      idSamp =grep(sampType[sampTypeCounter],RawMS1$SampleType) 
      
      datZoneP = RawMS1[,DatStartId:ncol(RawMS1)] ### datazone excludes "HybControlElution"
      
      datZone = RawMS1[idSamp,DatStartId:ncol(RawMS1)]  ###dataZone only include interested sample type
      
      if(length(idSamp)==1) {SampTypeRFU = matrix(as.matrix(datZone),nrow=1)
      SampTypeMedian = SampTypeRFU
      SampMedianNorm = as.matrix(SampTypeMedian/SampTypeRFU)}
      
      else {SampTypeRFU = datZone
      SampTypeMedian = apply(SampTypeRFU,2,median)
      SampMedianNorm = t(apply(SampTypeRFU,1,function(x){SampTypeMedian/x}))}
      
      Dilute = Dilution[which(Dilution!="0")] ###SampMedianNorm & Dilute: ncol(SampMedianNorm)=length(Dilute)
      uniqDilute = levels(factor(Dilute))
      
      for (idDilute in (1:length(uniqDilute))){
        
        DataDiluteID = which(Dilute==uniqDilute[idDilute])
        
        MedianDiluteSingle = matrix(apply(SampMedianNorm[,DataDiluteID],1,median),ncol=1)
        
        DataDiluteNormT = apply(datZone[,DataDiluteID],2,function(x){MedianDiluteSingle*x})
        
        if (idDilute==1){DataDiluteNorm = DataDiluteNormT}
        else{DataDiluteNorm = cbind(DataDiluteNorm,DataDiluteNormT)}
      }
      
      if (sampTypeCounter==1) {RawMStemp = cbind(RawMS[idSamp,1:(DatStartIdP-1)],RawMS[idSamp,idHyb],DataDiluteNorm)}
      else{
        RawMStempT = cbind(RawMS[idSamp,1:(DatStartIdP-1)],RawMS[idSamp,idHyb],DataDiluteNorm)
        RawMStemp = rbind(RawMStemp,RawMStempT)
      }
      RawMStemp2 = rbind(RawMStemp,RawMS[which(RawMS$SampleType =="Calibrator"),])
    }
    
    Platelist[[plateCounter]] = RawMStemp2
  }
  
  MySomaTemp = getMySoma(Platelist)
  
  ###up to date MySomaTemp has different indexing from Raw. We make them consistent
  rowOrderName = rownames(RawM)
  rowOrder = vector(mode="numeric",length=nrow(MySomaTemp))
  for (j in 1:nrow(MySomaTemp)){
    rowOrder[j] = which(rownames(MySomaTemp) %in% rowOrderName[j])
  }
  
  colOrderName = colnames(RawM)
  colOrder = vector(mode="numeric",length=ncol(MySomaTemp))
  for (k in 1:ncol(MySomaTemp)){
    colOrder[k] = which(colnames(MySomaTemp) %in% colOrderName[k])
  }
  
  MySoma = MySomaTemp[rowOrder,colOrder]
  return(MySoma)
}


### user define Funlist: any combination of normalisation in a interested order
UserNorm <- function(Funlist,RawM){
  for (FunCounter in 1:length(Funlist)){
    f <- Funlist[[FunCounter]]
    MySoma = f(RawM)
    RawM = MySoma
  }
  return(MySoma)
}


### compare correlation coefficiency among two RFUs
CompTWO <- function(SomaM,MySoma){
  
  DatStartId1 <- which(colnames(SomaM)=="CRYBB2.10000.28")
  DatStartId2 <- which(colnames(MySoma)=="CRYBB2.10000.28")
  
  rowOrderName = rownames(SomaM)
  rowOrder = vector(mode="numeric",length=nrow(SomaM))
  for (j in 1:nrow(SomaM)){
    rowOrder[j] = which(rownames(MySoma) %in% rowOrderName[j])
  }
  
  colOrderName = colnames(SomaM)[DatStartId1:ncol(SomaM)]
  colOrder = vector(mode="numeric",length=length(colOrderName))
  for (k in 1:length(colOrder)){
    colOrder[k] = which(colnames(MySoma)[DatStartId2:ncol(MySoma)] %in% colOrderName[k])
  }
  
  MySomaDone = MySoma[rowOrder,c(1:DatStartId2-1,colOrder+DatStartId2-1)]
  
  ### check row names colnames matching between two matrices
  # all(rownames(MySomaDone)==rownames(SomaM))
  # all(colnames(MySomaDone)[DatStartId2:ncol(MySomaDone)]==colnames(SomaM)[DatStartId1:ncol(SomaM)])
  
  corShip = vector(mode = "numeric", length=length(colOrderName)) ### correlation between my calculation and Adat
  for(dd in 1:length(corShip)){
    corShip[dd] = cor(MySomaDone[,dd+DatStartId2-1],SomaM[,dd+DatStartId1-1])
  }
  
  plot(corShip)
  return(corShip)
}

### divide RawM into two disease groups, analyse individually.clinicType:"OA"/"Injury"
ExtractClinicG = function(RawM,inputfile,trancheT){
  if(trancheT==1){metadata <- read_excel(inputfile,range="A2:Q438",col_names=TRUE)
  }else{metadata <- read_excel(inputfile,range="A2:O618",col_names=TRUE)}
  
  CohortInfor = metadata$`Cohort name`
  GroupInfor = metadata$Group
  bloodStainInfor = metadata$`Grading of SF bloodstaining prior to centrifugation  (if known)`
  bloodStainInfor[bloodStainInfor=="-"]<-NA
  
  ###sampleAge of tranche2 excel, we need to format: Data->text to column + general   
  ###01 JAN 2021 = 44197 as our up to date
  sampleAgeInfor <- 44197 - metadata$`Date of biological sampling`
  
  patientGenderInfor = metadata$`Patient gender`
  patientAgeInfor = metadata$`Patient age at Baseline SF sample`
  
  stepupIDTRaw = matrix(NA,ncol=1,nrow=nrow(metadata))
  CohortTRawM = matrix(NA,ncol=1,nrow=nrow(RawM))
  colnames(CohortTRawM) = "Corhort"
  GroupTRawM = matrix(NA,ncol=1,nrow=nrow(RawM))
  colnames(GroupTRawM) = "diseaseGroup"
  bloodStainTRawM = matrix(NA,ncol=1,nrow=nrow(RawM))
  colnames(bloodStainTRawM) = "bloodStain"
  sampleAgeTRaw = matrix(NA,ncol=1,nrow=nrow(RawM))
  colnames(sampleAgeTRaw) = "sampleAge"
  patientGenderTRaw = matrix(NA,ncol=1,nrow=nrow(RawM))
  colnames(patientGenderTRaw) = "patientGender"
  patientAgeTRaw = matrix(NA,ncol=1,nrow=nrow(RawM))
  colnames(patientAgeTRaw) = "patientAge"
  
  ### correct the SIN mistakes between clinic excel file and soma file
  clinicSIN <- metadata$`STEpUP Sample Identification Number (SIN)`
  ss = matrix(NA,ncol=1,nrow=length(clinicSIN))
  for (i in 1:length(ss)){
    hyphC <- length(which(grepl("-",strsplit(clinicSIN[i],"")[[1]])))
    sss <- strsplit(clinicSIN[i],"-")
    correct = vector(length =hyphC+1)
    for(j in 1:(hyphC+1)){
      if(j==2){correct[2] = sss[[1]][3]
      }else if(j==3){correct[3]=sss[[1]][2]
      }else{correct[j]=sss[[1]][j]}
    }
    ss[i] <- paste(correct,collapse="-")
  }
  STEPupName <- ss
  
  #special for tranche1 data:  
  STEPupName[STEPupName == "STEP1409-F-V1-HT1"] <- "STEP1409F-V1-HT1"
  
  
  for (spCounter in 1:nrow(metadata)){
    selPatient = which(RawM$SampleId == STEPupName[spCounter])
    CohortTRawM[selPatient] = CohortInfor[spCounter]
    GroupTRawM[selPatient]=GroupInfor[spCounter]
    bloodStainTRawM[selPatient] = bloodStainInfor[spCounter]
    sampleAgeTRaw[selPatient] = sampleAgeInfor[spCounter]
    patientGenderTRaw[selPatient] = patientGenderInfor[spCounter]
    patientAgeTRaw[selPatient] = patientAgeInfor[spCounter]
    stepupIDTRaw[selPatient] = STEPupName[spCounter]
  } ###RawM row order didn't change, just added more information from STEpUP_QCData_Tranche1.xlsx
  
  
  ### up to here integrate disease group type into RawM, RawM$SampleType including Pool type, disease sample type
  stepupIDTRaw[stepupIDTRaw == "STEP1409F-V1-HT1"] <- "STEP1409-F-V1-HT1"
  
  BioMeta = cbind(RawM$PlateId,stepupIDTRaw,GroupTRawM,CohortTRawM,bloodStainTRawM,sampleAgeTRaw,patientGenderTRaw,patientAgeTRaw)
  colnames(BioMeta)[1:2] = c("PlateID","STEpUPID")
  rownames(BioMeta) = rownames(RawM)
  
  return(BioMeta) 
}

### k batch for KNN testing: input expression data and batch matrix, return rejection rate matrix
KNNtest <- function(exprDat_norm,BioMeta,myRound,kBatch){
  
  ### k nearest neighborhood batch effect test, average among myRound, for kBatch. Per round includes different sampSizeS, and kNearS
  BatchEffectM = vector(mode="list",length=myRound)
  
  pc_norm <- prcomp(exprDat_norm,scale = TRUE)
  topPCn <- which(get_eigenvalue(pc_norm)$cumulative.variance.percent>80)[1]
  distPairs = pc_norm$x[,1:topPCn]
  DisM = as.matrix(dist(distPairs, method = "euclidean", diag = TRUE, upper = TRUE))
  
  sampSizeS = seq(10,500,by=100)  ### different percentage of 1%~25% samples for batch effect test (based on the paper)
  kNearS = seq(10,500,by=100)          ### different neiboughood sizes, based on testing, >250 no rejection at all
  
  for (roundCount in 1:myRound){
    positiveRate = matrix(NA,nrow=length(kNearS),ncol=length(sampSizeS))
    batchTest = vector(mode="list",length=kBatch)
    
    for (batchCf in 1:kBatch){ 
      tblWh = table(BioMeta[,batchCf])
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
            
            neiNameList = rownames(BioMeta)
            
            for (neiCounter in 1:(kNear-1)){ 
              idInMetaR[neiCounter] = which(neiNameList %in% neiborOriID[neiCounter]) ### meta for tested sample
            }
            
            neiDat = BioMeta[idInMetaR,batchCf]
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
            
            ### chi square based multinomial test, power is less than above one
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
    BatchEffectM[[roundCount]] = sapply(batchTest,mean)
  }
  
  BatchEffectMM = matrix(unlist(BatchEffectM),nrow=2,byrow=TRUE)
  meanM = apply(BatchEffectMM,2,mean)
  sdM = apply(BatchEffectMM,2,sd)
  BatchEffect = rbind(meanM,sdM)
  colnames(BatchEffect) = colnames(BioMeta)
  
  return(BatchEffect)
}
