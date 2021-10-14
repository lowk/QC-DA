### filter to exclude buffer|calibrator, select only clinical samples and human proteins
### input RFU matrix after normalisation steps, output filtered RFU matrix only with human data 
# filterHM <- function(MySoma,BioMeta){
#   HMpro <- which(!grepl("HybControlElution|Non",colnames(MySoma)))
#   HMsam <- which(grepl("Sample",MySoma[,"SampleType"]))
#   MySomaDone <- MySoma[HMsam,HMpro]
#   BioMetaDone <- BioMeta[HMsam,]
#   return(list(MySomaDone,BioMetaDone))
# }
filterHM <- function(MySoma){
  HMpro <- which(!grepl("HybControlElution|Non",colnames(MySoma)))
  HMsam <- which(grepl("Sample",MySoma[,"SampleType"]))
  MySomaDone <- MySoma[HMsam,HMpro]
  return(MySomaDone)
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



### user define Funlist: any combination of normalisation in a interested order
UserNorm <- function(Funlist,RawM){
  for (FunCounter in 1:length(Funlist)){
    f <- Funlist[[FunCounter]]
    MySoma = f(RawM)
    RawM = MySoma
  }
  return(MySoma)
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


MIDNORM = function(RawM,sampType){ ###sampType should be a vector subset of c("Calibrator","Buffer","Sample")
  
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
    
    for (sampTypeCounter in 1:length(sampType)){
      
      idSamp = which(RawMS1$SampleType == sampType[sampTypeCounter])
      
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


### function "CompTWO": conveniently compare any two Somalogic RFU matrix with the same dimensions. SomaM and MySoma are any two RUFs, with same dimensions. 
### output non matching rownames and colnames, and also use corThresh to set below which correlation coefficient, the protein names will be displayed. 

CompTWO <- function(SomaM,MySoma,corThresh){
  
  DatStartId1 <- which(colnames(SomaM)=="CRYBB2.10000.28")
  DatStartId2 <- which(colnames(MySoma)=="CRYBB2.10000.28")
  
  rowOrderName = rownames(SomaM)
  rowOrder = vector(mode="numeric",length=nrow(SomaM))
  rowKK=1
  notMatchRow=vector() ### rowKK, notMatchRow: extend code compatibility, when there are none matching records
  for (j in 1:nrow(SomaM)){
    if(!any(rownames(MySoma)==rowOrderName[j])){notMatchRow[rowKK]=j
    rowKK=rowKK+1}
    else{rowOrder[j] = which(rownames(MySoma) %in% rowOrderName[j])}
  }
  print(paste("row name not match: ",rowOrderName[notMatchRow]))
  
  colOrderName = colnames(SomaM)[DatStartId1:ncol(SomaM)]
  colOrder = vector(mode="numeric",length=length(colOrderName))
  colKK=1
  notMatchCol=vector()
  for (k in 1:length(colOrder)){
    if(!any(colnames(MySoma)[DatStartId2:ncol(MySoma)]==colOrderName[k])){notMatchCol[colKK] = k
    colKK=colKK+1}
    else{colOrder[k] = which(colnames(MySoma)[DatStartId2:ncol(MySoma)] %in% colOrderName[k])}
  }
  print(paste("column name not match in Soma: ",colOrderName[notMatchCol]))
  
  # MySomaDone = MySoma[rowOrder,c(1:DatStartId2-1,colOrder+DatStartId2-1)]
  ### check row names colnames matching between two matrices
  # all(rownames(MySomaDone)==rownames(SomaM))
  # all(colnames(MySomaDone)[DatStartId2:ncol(MySomaDone)]==colnames(SomaM)[DatStartId1:ncol(SomaM)])
  
  corShip = vector(mode = "numeric", length=length(colOrder)) ### correlation between my calculation and Adat
  for(dd in 1:length(corShip)){
    if(!(dd %in% notMatchCol)){
      if(length(notMatchRow)!=0){corShip[dd] = cor(MySoma[rowOrder,colOrder[dd]+DatStartId2-1],SomaM[!notMatchRow,dd+DatStartId1-1])}
      else{corShip[dd] = cor(MySoma[rowOrder,colOrder[dd]+DatStartId2-1],SomaM[,dd+DatStartId1-1])}
    }
    else{corShip[dd]=NA}
  }
  
  plot(corShip[!is.na(corShip)],ylab="correlation coefficient")
  print(paste("correlation coefficient between minimum value of ", min(corShip[!is.na(corShip)])," to user selected threshold of",corThresh, " :",colnames(SomaM)[DatStartId1:ncol(SomaM)][which(corShip<corThresh)]))
  return()
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
  if(trancheT==1){sampleAgeInfor <- 2021 - metadata$`Date of biological sampling`
  }else{sampleAgeInfor <- floor((44197 - metadata$`Date of biological sampling`)/365)}
  
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

### initialise required metadata for normalisation methods
initQCnorm <- function(inputfile1,inputfile2){
  
  RawM <- read.adat(inputfile1) 
  ### read in Col^MetaTable
  
  con1 = file(inputfile2, "r")
  
  while(TRUE) {
    sline = readLines(con1, n=1)
    
    if(length(sline) == 0){
      print("Meta data not intact")
      break}
    
    slineL = strsplit(sline,'\t')[[1]]
    
    if(length(slineL) > 28){
      getStartTerm = which(slineL!="")[1]
      getStartId = getStartTerm+1
      
      if(slineL[getStartTerm]=="TargetFullName"){ TargetFullName <<- slineL[getStartId:length(slineL)]}
      if(slineL[getStartTerm]=="Target"){ Target <<- slineL[getStartId:length(slineL)]}
      if(slineL[getStartTerm]=="UniProt"){ UniProt <<- slineL[getStartId:length(slineL)]}
      if(slineL[getStartTerm]=="EntrezGeneID"){ EntrezGeneID <<- slineL[getStartId:length(slineL)]}
      if(slineL[getStartTerm]=="EntrezGeneSymbol"){ EntrezGeneSymbol <<- slineL[getStartId:length(slineL)]}
      if(slineL[getStartTerm]=="Type"){Type <<- slineL[getStartId:length(slineL)]}
      if(slineL[getStartTerm]=="Dilution"){Dilution <<- slineL[getStartId:length(slineL)]}
      # if(slineL[getStartTerm]=="medNormRefSMP_ReferenceRFU"){
      #   medNormRefSMP <<- as.numeric(slineL[getStartId:length(slineL)])
      if(slineL[getStartTerm]=="PlateScale_Reference"){
        PlateScale_Reference <<- as.numeric(slineL[getStartId:length(slineL)])
        
        break} 
    }
  }
  close(con1)
  
  DatStartId = which(colnames(RawM) == "CRYBB2.10000.28")
  ColTable <- cbind(as.matrix(colnames(RawM)[DatStartId:ncol(RawM)]),as.matrix(UniProt),as.matrix(EntrezGeneID),as.matrix(EntrezGeneSymbol),as.matrix(TargetFullName),as.matrix(Target))
  colnames(ColTable) <- list("Protein Name", "UniPro ID", 'EntrezGeneID',"EntrezGeneSymbol","TargetFullName","Target")
  
  EntrezID = vector(mode="integer",length=nrow(ColTable)) ### ColTable: change multiple geneIDs to the first 1
  for (EntrezID in 1:nrow(ColTable)){
    IDsymlist = strsplit(ColTable[,"EntrezGeneSymbol"][EntrezID]," ")
    IDlist = strsplit(ColTable[,"EntrezGeneID"][EntrezID]," ")
    if(length(IDsymlist[[1]])>1){ColTable[,"EntrezGeneSymbol"][EntrezID]=IDsymlist[[1]][1]
    ColTable[,"EntrezGeneID"][EntrezID]=IDlist[[1]][1]}
  } 
  
  return(ColTable)
}

### combine tranches data and sva::combat function to remove batches
MyCombat <- function(RFU1,RFU2,noCombat){
  
  Done1 <- filterHM(RFU1)
  Done2 <- filterHM(RFU2)
  tranche <- c(rep(1,nrow(Done1)),rep(2,nrow(Done2)))
  MySomaPlate <- as.matrix(c(Done1[,"PlateId"],Done2[,"PlateId"]),ncol=1)
  PlateBatch <- GetPlateBatch(MySomaPlate)
  BatchMeta = cbind(tranche,PlateBatch)
  
  DoneAll <- rbind(Done1[,c(1:13,which(colnames(Done1) == "CRYBB2.10000.28"):ncol(Done1))],Done2[,c(1:13,which(colnames(Done2) == "CRYBB2.10000.28"):ncol(Done2))])
  DoneAllDat <- rbind(Done1[,which(colnames(Done1) == "CRYBB2.10000.28"):ncol(Done1)],Done2[,which(colnames(Done2) == "CRYBB2.10000.28"):ncol(Done2)])
  rownames(BatchMeta)=rownames(DoneAll)
  colnames(BatchMeta) = c("Tranche Batch","Plate Batch")
  PlateBatch <- as.vector(PlateBatch) ### Combat argument requirement
  
  if(noCombat==1){combat_all_batch = t(log(DoneAllDat))
  combat_all = data.frame(log(DoneAllDat))}
  else{combat_all_batch <- sva::ComBat(t(log(DoneAllDat)), PlateBatch, mod=NULL, par.prior = TRUE, prior.plots = FALSE)
  combat_all = data.frame(t(combat_all_batch))}
  
  batchTest <- KNNtest(t(combat_all_batch),BatchMeta,2,2)
  
  return(list(combat_all,batchTest,BatchMeta))
}


PlotPCA <- function(exprDat,BatchMeta,topPC,confounderC,titleMessage){
  pc_norm <- prcomp(as.matrix(exprDat),scale = TRUE)
  PlotDat = data.frame(cbind(pc_norm$x[,1:topPC],BatchMeta))
  
  gp <- ggpairs(PlotDat, columns=1:topPC, aes(color= as.factor(PlotDat[,topPC+confounderC])),
                title=titleMessage,
                diag=list(continuous=wrap("densityDiag",alpha=0.4)),
                lower=list(continuous = wrap("points",alpha=0.9,size=0.1)),
                upper = list(continuous = "blank"),
                legend = c(1,1)) + labs(fill = colnames(BatchMeta)[confounderC])
  return(gp)
}

PlotUmap <- function(exprDat,BatchMeta,confounderC,titleMessage){
  myUmap = umap(t(exprDat))
  df <- data.frame(x = myUmap$data[,"X1"],
                   y = myUmap$data[,"X2"],
                   WhichBatch = factor(BatchMeta[,confounderC]))
  colnames(df) = c("D1","D2",colnames(BatchMeta)[confounderC])
  ggplot(df, aes(x=D1, y=D2, color = df[,3])) +
    geom_point() + labs(title=titleMessage,color=colnames(BatchMeta)[confounderC])
}

CVbreak <- function(RFU1,RFU2,clinicType,exprDat,titleMessage){
  par(mfrow=c(1,2))
  
  Done1 <- filterHM(RFU1)
  Done2 <- filterHM(RFU2)
  
  DoneAll <- rbind(Done1[,c(1:13,which(colnames(Done1) == "CRYBB2.10000.28"):ncol(Done1))],Done2[,c(1:13,which(colnames(Done2) == "CRYBB2.10000.28"):ncol(Done2))])
  exprDat_Mr <- data.frame(cbind(DoneAll[,c(1:13)],exprDat)) 
  calib_normM <- exprDat_Mr[grep(paste(clinicType,"POOL"),exprDat_Mr$SampleId),]
  calib_norm <- as.matrix(calib_normM[,-c(1:(which(colnames(calib_normM)=="CRYBB2.10000.28")-1))])
  calibIDs <- calib_normM$SampleId
  
  if(clinicType == "OA"){suffix = "/29"
  freezeThaw <-  calibIDs %in% paste0("OA POOL-HT-",c(6,25,26),suffix)
  }else{suffix = "/25"
  freezeThaw <-  calibIDs %in% paste0("INJ POOL-HT-",c(6,25),suffix)}
  
  acrossPlates <- calibIDs %in% paste0(clinicType," POOL-HT-",c(1:5,7:13),suffix)
  withinPlates <- calibIDs %in% paste0(clinicType," POOL-HT-",c(1,6),suffix)
  
  
  temp1 <- (apply(calib_norm[acrossPlates,],2,sd)/apply(calib_norm[acrossPlates,],2,mean))
  temp2 <- (apply(calib_norm[withinPlates,],2,sd)/apply(calib_norm[withinPlates,],2,mean))
  
  plot(100*quantile(temp1,seq(0,1,length.out=100)),seq(0,1,length.out=100),type="l",xlim=c(0,10),xlab=paste("%CV",clinicType,"Group"),ylab="Cumulative total",main="",lwd=2,cex.lab=1.5)
  lines(100*quantile(temp2,seq(0,1,length.out=100)),seq(0,1,length.out=100),type="l",lwd=2,col="red")
  lines(c(0,quantile(temp1,0.8))*100,c(0.8,0.8),lty=2)
  lines(c(quantile(temp1,0.8),quantile(temp1,0.8))*100,c(0.8,0),lty=2)
  lines(c(0,quantile(temp2,0.8))*100,c(0.8,0.8),lty=2,col="red")
  lines(c(quantile(temp2,0.8),quantile(temp2,0.8))*100,c(0.8,0),lty=2,col="red")
  legend(8,0.4,c("Across Plates","Within Plates"),lwd=2,col=c("blue","red"),xpd=T,cex=0.75)
  
  ### across  tranche
  tempAct = matrix(0,nrow=5,ncol=ncol(calib_norm))
  
  for(i in 1:5){
    acrossTranche <- which(calibIDs %in% paste0(clinicType," POOL-HT-",c(i,i+6),suffix))
    tempAct[i,] <- apply(calib_norm[acrossTranche,],2,sd)/apply(calib_norm[acrossTranche,],2,mean)
  }
  tempAct <- apply(tempAct,2,mean)
  
  ### within tranche
  withinTranche1 <- calibIDs %in% paste0(clinicType," POOL-HT-",c(1:5),suffix)
  withinTranche2 <- calibIDs %in% paste0(clinicType," POOL-HT-",c(7:13),suffix)
  tempAct11= apply(calib_norm[withinTranche1,],2,sd)/apply(calib_norm[withinTranche1,],2,mean)
  tempAct22= apply(calib_norm[withinTranche2,],2,sd)/apply(calib_norm[withinTranche2,],2,mean)
  tempAct2 = (tempAct11+tempAct22)/2
  
  
  plot(100*quantile(tempAct,seq(0,1,length.out=100)),seq(0,1,length.out=100),type="l",xlim=c(0,10),xlab=paste("%CV",clinicType,"Group"),ylab="Cumulative total",main="",lwd=2,cex.lab=1.5)
  lines(c(0,quantile(tempAct,0.8))*100,c(0.8,0.8),lty=2)
  lines(c(quantile(tempAct,0.8),quantile(tempAct,0.8))*100,c(0.8,0),lty=2,)
  lines(100*quantile(tempAct2,seq(0,1,length.out=100)),seq(0,1,length.out=100),type="l",col="red")
  lines(c(0,quantile(tempAct2,0.8))*100,c(0.8,0.8),lty=2,col="red")
  lines(c(quantile(tempAct2,0.8),quantile(tempAct2,0.8))*100,c(0.8,0),lty=2,col="red")
  
  legend(6,0.4,c("Across Tranche","Within Tranche"),lwd=2,col=c("black","red"),xpd=T,cex=0.75)
  
  title(titleMessage)
  
  return()
}


### toTest1,toTest2 are biomarker lists from sandwich file and adat file individually
ExtVal <- function(exprDatM){
  metadata_xls <- read_excel("STEpUP_QCData_Tranche1.xlsx")
  temp1 <- data.frame(as.matrix(metadata_xls)[-1,1:17])
  names(temp1) <- as.matrix(metadata_xls)[1,1:17]
  metadata <- temp1
  rownames(metadata) <- metadata$`STEpUP Sample Identification Number (SIN)`
  stepupID <- exprDatM$SampleId[grep("STEP",exprDatM$SampleId)]
  stepupID[stepupID == "STEP1409F-V1-HT1"] <- "STEP1409-F-V1-HT1" #fix an apparant typo
  stepupID <- gsub("F-V1","V1-F",stepupID) #ID ordering seems to have changed?
  plateID <- exprDatM$PlateId[grep("STEP",exprDatM$SampleId)]
  exprDat_norm <- exprDatM[grep("STEP",exprDatM$SampleId),which(colnames(exprDatM)=="CRYBB2.10000.28"):ncol(exprDatM)]
  #reorder meta-data
  metadata_reord <- metadata[stepupID,]

  sandwich_master_xls_2 <- read_excel("Masterlist.xlsx",sheet=2)
  #process Ben data, averaging across replicates
  temp1 <- data.frame(sandwich_master_xls_2)
  temp2 <- (temp1[temp1$replicate == 1,-c(1:2)] + temp1[temp1$replicate == 2,-c(1:2)])/2
  sandwich_master <- data.frame(PIN=temp1[temp1$replicate == 1,1],temp2)
  rownames(sandwich_master) <- sandwich_master$PIN

  toTest1 <- c("mcp1bl","il6bl","il8bl","mmp3bl","activinabl","tsg6bl","timp1bl","tgfb1bl","fgf2bl")
  toTest2 <- c("CCL2.2578.67","IL6.4673.13","CXCL8.3447.64","MMP3.2788.55","INHBA.13738.8","TNFAIP6.5036.50","TIMP1.2211.9","TGFB1.2333.72","FGF2.3025.50")

  CorData_norm = data.frame(matrix(NA,nrow=length(toTest1),ncol=5))
  names(CorData_norm) <- c("SandwichName","SomaName","cor","Pvalue","N")

  for (compCounter in 1:length(toTest1)){
    CorData_norm[compCounter,] <- getpars(sandwich_master,exprDat_norm,toTest1[compCounter],toTest2[compCounter])
  }
  return(CorData_norm)
}

getpars <- function(sandwich_master,exprDat_norm,par1,par2) {
  temp1 <- sandwich_master[as.character(metadata_reord$`STEpUP Participant Identification Number (PIN)`),]

  if (sum(!(is.na(temp1[,par1] + exprDat_norm[,par2]))) == 0) return(c(par1,par2,NA,NA,0))
  md <- cor.test(temp1[,par1],exprDat_norm[,par2])
  c(par1,par2,unlist(md[c("estimate","p.value")]),2+md$parameter)

}

getpars2 <- function(temp1,par1,par2,keep) {
  if (sum(!(is.na(temp1[keep,par1] + exprDat_norm[keep,par2]))) == 0) return(c(par1,par2,NA,NA,0))
  md <- cor.test(temp1[keep,par1],exprDat_norm[keep,par2])
  c(par1,par2,unlist(md[c("estimate","p.value")]),2+md$parameter)
}

