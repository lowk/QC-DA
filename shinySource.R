setReactiveMySoma <- function(UserInput,RawM){
  
  userSt = strsplit(UserInput,",")[[1]]
  if(userSt =="RawM"){MySoma = RawM}
  else{
    numStep = length(userSt) 
    Funlist = vector(mode="list",length=numStep)
    for (counterStep in 1:numStep){
      Funlist[[counterStep]] = get(userSt[counterStep])
    }
    MySoma = UserNorm(Funlist,RawM)}
  
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
  
  DatStartId <- which(colnames(RawM)=="CLI")+1  ###for calculation convenience, extract data zone only
  
  HybId = which(Type == "Hybridization Control Elution")
  
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
  
  DatStartId <- which(colnames(RawM)=="CLI")+1  ###for calculation convenience, extract data zone only
  
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
  
  DatStartId <- which(colnames(RawM)=="CLI")+1  ###for calculation convenience, extract data zone only
  
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

### 4. Median Normalization. Input: RFU whole data frame, user defined sample type;
###Output: median normalised (row-wise and column-wise) RFU whole data frame
# 
# MIDNORM = function(RawM,caliCase){ ###caliCase control which sample type to be applied MidNorm
#   
#   PlateIdUni = levels(factor(RawM$PlateId))
#   
#   Platelist = list()
#   
#   for (plateCounter in 1:length(PlateIdUni)){
#     
#     PlateIdSg = which(RawM$PlateId == PlateIdUni[plateCounter])
#     
#     RawMS = RawM[PlateIdSg,] ### single plate
#     
#     idHyb = which(grepl("HybControlElution",colnames(RawMS))==TRUE)
#     idNonHyb = which(!grepl("HybControlElution",colnames(RawMS))==TRUE)
#     RawMS1 = RawMS[,idNonHyb] ###RawM1 to track "Ratio of Normalization Median to Sample Value"
#     DatStartId = which(colnames(RawMS1) == "CLI") +1
#     DatStartIdP = which(colnames(RawMS) == "CLI") +1
#     
#     sampTypePre = levels(factor(RawMS$SampleType))
#     if (caliCase == 1){sampType = "Calibrator"}
#     else{sampType = sampTypePre[which(sampTypePre!="Calibrator")]}    
#     
#     for (sampTypeCounter in 1:length(sampType)){
#       
#       idSamp = which(RawMS1$SampleType == sampType[[sampTypeCounter]])
#       
#       datZoneP = RawMS1[,DatStartId:ncol(RawMS1)] ### datazone excludes "HybControlElution"
#       
#       datZone = RawMS1[idSamp,DatStartId:ncol(RawMS1)]  ###dataZone only include interested sample type
#       
#       if(length(idSamp)==1) {SampTypeRFU = matrix(datZone,nrow=1)
#       SampTypeMedian = SampTypeRFU
#       SampMedianNorm = as.matrix(SampTypeMedian/SampTypeRFU)}
#       
#       else {SampTypeRFU = datZone
#       SampTypeMedian = apply(SampTypeRFU,2,median)
#       SampMedianNorm = t(apply(SampTypeRFU,1,function(x){SampTypeMedian/x}))}
#       
#       Dilute = Dilution[which(Dilution!="0")] ###SampMedianNorm & Dilute: ncol(SampMedianNorm)=length(Dilute)
#       uniqDilute = levels(factor(Dilute))
#       
#       for (idDilute in (1:length(uniqDilute))){
#         
#         DataDiluteID = which(Dilute==uniqDilute[idDilute])
#         
#         MedianDiluteSingle = matrix(apply(SampMedianNorm[,DataDiluteID],1,median),ncol=1)
#         
#         DataDiluteNormT = apply(datZone[,DataDiluteID],2,function(x){MedianDiluteSingle*x})
#         
#         if (idDilute==1){DataDiluteNorm = DataDiluteNormT}
#         else{DataDiluteNorm = cbind(DataDiluteNorm,DataDiluteNormT)}
#       }
#       
#       if (sampTypeCounter==1) {RawMStemp = cbind(RawMS[idSamp,1:(DatStartIdP-1)],RawMS[idSamp,idHyb],DataDiluteNorm)}
#       else{
#         RawMStempT = cbind(RawMS[idSamp,1:(DatStartIdP-1)],RawMS[idSamp,idHyb],DataDiluteNorm)
#         RawMStemp = rbind(RawMStemp,RawMStempT)
#       }
#     }
#     
#     Platelist[[plateCounter]] = RawMStemp  
#   }
#   
#   MySomaTemp = getMySoma(Platelist)  
#   
#   ###up to date MySomaTemp has different indexing from Raw. We make them consistent
#   rowOrderName = rownames(RawM)
#   rowOrder = vector(mode="numeric",length=nrow(MySomaTemp))
#   for (j in 1:nrow(MySomaTemp)){
#     rowOrder[j] = which(rownames(MySomaTemp) %in% rowOrderName[j])
#   }
#   
#   colOrderName = colnames(RawM)
#   colOrder = vector(mode="numeric",length=ncol(MySomaTemp))
#   for (k in 1:ncol(MySomaTemp)){
#     colOrder[k] = which(colnames(MySomaTemp) %in% colOrderName[k])
#   }
#   
#   MySoma = MySomaTemp[rowOrder,colOrder]
#   return(MySoma)
# }

MIDNORM = function(RawM){ ###caliCase control which sample type to be applied MidNorm
  
  PlateIdUni = levels(factor(RawM$PlateId))
  
  Platelist = list()
  
  for (plateCounter in 1:length(PlateIdUni)){
    
    PlateIdSg = which(RawM$PlateId == PlateIdUni[plateCounter])
    
    RawMS = RawM[PlateIdSg,] ### single plate
    
    idHyb = which(grepl("HybControlElution",colnames(RawMS))==TRUE)
    idNonHyb = which(!grepl("HybControlElution",colnames(RawMS))==TRUE)
    RawMS1 = RawMS[,idNonHyb] ###RawM1 to track "Ratio of Normalization Median to Sample Value"
    DatStartId = which(colnames(RawMS1) == "CLI") +1
    DatStartIdP = which(colnames(RawMS) == "CLI") +1
    
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
    DatStartId = which(colnames(RawMS1) == "CLI") +1
    DatStartIdP = which(colnames(RawMS) == "CLI") +1
    
    sampType = "Calibrator"
    
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
    DatStartId = which(colnames(RawMS1) == "CLI") +1
    DatStartIdP = which(colnames(RawMS) == "CLI") +1
    
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


### initialise required metadata for normalisation methods
initQCnorm <- function(inputfile1,inputfile2){
  
  RawM <- read.adat(inputfile1)
  RawM$SampleType[grep("OA",RawM$SampleId)]="OApoolSample" ### facilitate median norm on QC
  RawM$SampleType[grep("INJ",RawM$SampleId)]="INJpoolSample"
  target = colnames(RawM)[25:ncol(RawM)]  
  ### read in Col^MetaTable
  ###use inputfile1, there is no PlateScale_Reference information, so we use inputfile2 considering all the row and column indices are consistent
  
  con1 = file(inputfile2, "r")
  
  while(TRUE) {
    sline = readLines(con1, n=1)
    
    if(length(sline) == 0){
      print("Meta data not intact")
      break}
    
    slineL = strsplit(sline,'\t')[[1]]
    
    if(length(slineL) > 29){
      
      if(slineL[30]=="Type"){Type <<- slineL[31:length(slineL)]}
      if(slineL[30]=="Dilution"){Dilution <<- slineL[31:length(slineL)]}
      if(slineL[30]=="PlateScale_Reference"){
        PlateScale_Reference <<- as.numeric(slineL[31:length(slineL)])
        break}
      ### 30 is the beginning field of the ^COL_DATA;column names of meta data = SeqId +Target
    }
  }
  close(con1)
  
  return(RawM)
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
  
  DatStartId1 <- which(colnames(SomaM)=="NormScale_0_5")+1
  DatStartId2 <- which(colnames(MySoma)=="CLI")+1
  
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
  
  corShip = vector(mode = "numeric", length=length(colOrderName)) ### correlation between my calculation and Adat
  for(dd in 1:length(corShip)){
    corShip[dd] = cor(MySomaDone[dd+DatStartId2-1],SomaM[dd+DatStartId1-1])
  }
  
  plot(corShip)
  return(corShip)
}

### divide RawM into two disease groups, analyse individually.clinicType:"OA"/"Injury"
ExtractClinicG = function(RawM){
  metadata_xls <- read_excel("/Users/ydeng/Documents/QCstepOA/STEpUP_QCData_Tranche1.xlsx")
  temp1 <- data.frame(as.matrix(metadata_xls)[-1,1:17])
  names(temp1) <- as.matrix(metadata_xls)[1,1:17]
  metadata <- temp1
  PlateInfor = metadata$`Plate number`
  CohortInfor = metadata$`Cohort name`
  GroupInfor = metadata$Group
  bloodStainInfor = metadata$`Grading of SF bloodstaining prior to centrifugation  (if known)`
  bloodStainInfor[bloodStainInfor=="-"]<-NA
  sampleAgeInfor <- 2020 - as.numeric(metadata$`Date of biological sampling`)
  
  selPatientID = matrix(NA,ncol=1,nrow=nrow(metadata))
  PlateTRawM = matrix(NA,ncol=1,nrow=nrow(RawM))
  colnames(PlateTRawM) = "PlateID"
  GroupTRawM = matrix(NA,ncol=1,nrow=nrow(RawM))
  colnames(GroupTRawM) = "diseaseGroup"
  CohortTRawM = matrix(NA,ncol=1,nrow=nrow(RawM))
  colnames(CohortTRawM) = "Corhort"
  bloodStainTRawM = matrix(NA,ncol=1,nrow=nrow(RawM))
  colnames(bloodStainTRawM) = "bloodStain"
  sampleAgeTRawM = matrix(NA,ncol=1,nrow=nrow(RawM))
  colnames(sampleAgeTRawM) = "sampleAge"
  
  STEPupName <- gsub("V1-F","F-V1",metadata$`STEpUP Sample Identification Number (SIN)`)
  
  RawM$SampleId[which(RawM$SampleId == "STEP1409F-V1-HT1")] = "STEP1409-F-V1-HT1"
  
  for (spCounter in 1:length(selPatientID)){
    selPatientID = which(RawM$SampleId == STEPupName[spCounter])
    PlateTRawM[selPatientID]=PlateInfor[spCounter]
    GroupTRawM[selPatientID]=GroupInfor[spCounter]
    CohortTRawM[selPatientID] = bloodStainInfor[spCounter]
    bloodStainTRawM[selPatientID] = bloodStainInfor[spCounter]
    sampleAgeTRawM[selPatientID] = sampleAgeInfor[spCounter]
  }
  
  DatStartId <- which(colnames(RawM)=="CLI")+1
  RawM$SampleType[which(GroupTRawM=="OA")]="OASample"
  RawM$SampleType[which(GroupTRawM=="Injury")]="INJSample"
  
  ### up to now integrate disease group type into RawM, RawM$SampleType including Pool type, disease sample type
  
  # Total protein distribution for OA, INJ individually
  # histOA = apply(RawM[RawM$SampleType=="OASample",DatStartId:ncol(RawM)],1,sum)
  # hist(histOA,breaks=100,main=paste("Total protein distribution of OA"),cex.main=1)
  # histINJ = apply(RawM[RawM$SampleType=="INJSample",DatStartId:ncol(RawM)],1,sum)
  # hist(histINJ,breaks=100,main=paste("Total protein distribution of Injury"),cex.main=1)
  # 
  RawMList = list(mode="vector",length=2)
  RawMList[[1]] = RawM
  RawMList[[2]] = cbind(PlateTRawM,GroupTRawM,CohortTRawM,bloodStainTRawM,sampleAgeTRawM)
  rownames(RawMList[[2]]) <- rownames(RawM)
  
  return(RawMList) 
}


####asssessment
##1: Total protein check. Input normalised RFUs dataframe, output plots

TotalProCheck <- function(exprDat_norm,bloodStain,metadata_reord){
  totalProtein_norm <- apply(exprDat_norm,1,sum)
  par(mfrow=c(1,3))
  hist(totalProtein_norm,breaks=100,xlab="Total Protein",main="total protein distribution")
  plot(as.numeric(bloodStain),totalProtein_norm,xlab="Blood staining grade",ylab="Total Protein")
  abline(lm(totalProtein_norm ~ as.numeric(bloodStain)),col="red",lwd=2)
  text(1,max(totalProtein_norm)*0.9,paste("p =",formatC(cor.test(totalProtein_norm, as.numeric(bloodStain))$p.val,format="e",digits=3)),pos=4,col="red")
  text(1,max(totalProtein_norm)*0.8,paste("cor =",formatC(cor.test(totalProtein_norm, as.numeric(bloodStain))$estimate,format="e",digits=3)),pos=4,col="red")
  boxplot(totalProtein_norm ~ metadata_reord$Group,main=paste("p =",wilcox.test(log(totalProtein_norm) ~ metadata_reord$Group)$p.val),xlab="disease group",ylab="total protein within group")
  myplot <- recordPlot()
  return(myplot)
}

##2: Checks against calibrators. Input selected RFUs calib_norm, and clinicalType ("OA" or "INJ")
# pool %CV
CalibratorCheck <- function(calib_norm,clinicType,calibIDs){
  temp1 <- apply(calib_norm[grep(clinicType,calibIDs),],2,sd)/apply(calib_norm[grep(clinicType,calibIDs),],2,mean)
  #plot(100*quantile(temp1,seq(0,1,length.out=100)),seq(0,1,length.out=100),type="l",xlim=c(0,50),xlab=paste("%CV",clinicType),ylab="Cumulative total",main="Norm")
  plot(100*quantile(temp1,seq(0,1,length.out=100)),seq(0,1,length.out=100),type="l",xlim=c(0,30),xlab=paste("%CV",clinicType),ylab="Cumulative total",main="Norm")
  abline(h=0.8,lty=2)
  abline(v=100*quantile(temp1,0.8),lty=2)
  text(100*quantile(temp1,0.8)-4,0.2,quantile(temp1,0.8),pos=4,col="red")
  ptPart <- recordPlot()
  return(ptPart)
}

CalibratorCheckPlot <- function(calib_norm,calibIDs){
  par(mfrow=c(2,1))
  ptOA <- CalibratorCheck(calib_norm,"OA",calibIDs)
  ptINJ <- CalibratorCheck(calib_norm,"INJ",calibIDs)
  ptOA
  ptINJ
  myplot <- recordPlot()
  return(myplot)
}

##CV breakdowns:
CVbreak <- function(calib_norm,clinicType,calibIDs,calibPlates){
  if(clinicType == "OA"){suffix = "/29"}
  else{suffix = "/25"}
  
  acrossPlates <- calibIDs %in% paste0(clinicType," POOL-HT-",1:5,suffix)
  withinPlates <- calibIDs %in% paste0(clinicType," POOL-HT-",c(1,6),suffix)
  freezeThaw <-  calibIDs %in% paste0(clinicType," POOL-HT-",c(6,25),suffix)
  
  temp1 <- (apply(calib_norm[acrossPlates,],2,sd)/apply(calib_norm[acrossPlates,],2,mean))
  temp2 <- (apply(calib_norm[withinPlates,],2,sd)/apply(calib_norm[withinPlates,],2,mean))
  temp3 <- (apply(calib_norm[freezeThaw,],2,sd)/apply(calib_norm[freezeThaw,],2,mean))
  
  plot(100*quantile(temp1,seq(0,1,length.out=100)),seq(0,1,length.out=100),type="l",xlim=c(0,25),xlab=paste("%CV",clinicType),ylab="Cumulative total",main="Norm",lwd=2)
  lines(100*quantile(temp2,seq(0,1,length.out=100)),seq(0,1,length.out=100),type="l",lwd=2,col="red")
  lines(100*quantile(temp3,seq(0,1,length.out=100)),seq(0,1,length.out=100),type="l",lwd=2,col="blue")
  
  lines(c(0,quantile(temp1,0.8))*100,c(0.8,0.8),lty=2)
  lines(c(quantile(temp1,0.8),quantile(temp1,0.8))*100,c(0.8,0),lty=2)
  lines(c(0,quantile(temp2,0.8))*100,c(0.8,0.8),lty=2,col="red")
  lines(c(quantile(temp2,0.8),quantile(temp2,0.8))*100,c(0.8,0),lty=2,col="red")
  lines(c(0,quantile(temp3,0.8))*100,c(0.8,0.8),lty=2,col="blue")
  lines(c(quantile(temp3,0.8),quantile(temp3,0.8))*100,c(0.8,0),lty=2,col="blue")
  
  legend(15,0.4,c("Across","Within","Freeze"),lwd=2,col=c("black","red","blue"),xpd=T,cex=0.75)
  
  lines(c(0,quantile(temp1,0.8))*100,c(0.8,0.8),lty=2)
  lines(c(quantile(temp1,0.8),quantile(temp1,0.8))*100,c(0.8,0),lty=2)
  lines(c(0,quantile(temp2,0.8))*100,c(0.8,0.8),lty=2,col="red")
  lines(c(quantile(temp2,0.8),quantile(temp2,0.8))*100,c(0.8,0),lty=2,col="red")
  lines(c(0,quantile(temp3,0.8))*100,c(0.8,0.8),lty=2,col="blue")
  lines(c(quantile(temp3,0.8),quantile(temp3,0.8))*100,c(0.8,0),lty=2,col="blue")
  
  ##each of the five plates against the average of the rest
  plot(NA,xlim=c(0,50),ylim=c(0,1),xlab=paste("%CV",clinicType),ylab="Cumulative total",main="Norm,per plate against the others")
  printList=vector(mode="list",length=length(acrossPlates))
  for (i in 1:5){
    temp4 <- abs(calib_norm[which(acrossPlates)[i],] - apply(calib_norm[which(acrossPlates)[-i],],2,mean))/apply(calib_norm[acrossPlates,],2,mean)
    lines(100*quantile(temp4,seq(0,1,length.out=100)),seq(0,1,length.out=100),type="l",lwd=2,col=i)
    printList[[i]]=(c(calibPlates[which(acrossPlates)[i]],mean(temp4)))
  }
  legend(35,0.5,col=1:5,legend=calibPlates[acrossPlates],lwd=2,xpd=T,cex=0.75)
  return()
}

CVbreakPlot <- function(calib_norm,calibIDs,calibPlates){
  par(mfrow=c(2,2))
  CVbreak(calib_norm,"OA",calibIDs,calibPlates)
  CVbreak(calib_norm,"INJ",calibIDs,calibPlates)
  myplot<-recordPlot()
  return(myplot)
}

#Check 3: PCA
PCAglob <- function(MetaRaw,exprDat_norm,whichTplot){
  pc_norm <- prcomp(log10(as.matrix(exprDat_norm)),scale = TRUE)
  
  colnames(MetaRaw) <- c("PlateID","Disease.Group","Corhort","Blood.Stain","Sample.Age")
  PlotDat = data.frame(cbind(pc_norm$x[,1:10],MetaRaw))
  PlotDat[,1:10] = apply(PlotDat[,1:10],2,as.numeric)
  
  fullList = c("PlateID","Disease Group","Corhort","Blood Stain","Sample Age","Variation Explained by Principle Components")
  idTplot = which(fullList == whichTplot)
  
  if(idTplot==6){
    eig.val <- get_eigenvalue(pc_norm)
    egV = eig.val$eigenvalue[1:10] ###eigenvalues
    PCv = eig.val$variance.percent[1:10]     ### corresponding PC variance percentage
    PCcv = eig.val$cumulative.variance.percent[1:10]   ### corresponding cumulative PCv
    PCcount = seq(1:10)
    
    dat = as.data.frame(cbind(PCcount,egV,PCv,PCcv))
    colors <- c("PCv"="cyan","PCcv"="darkcyan")
    ptPCA <- ggplot(dat,aes(x=PCcount)) + geom_point(aes(y=PCv,color="PCv")) + geom_line(aes(y=PCv, color="PCv")) +
      geom_point(aes(y=PCcv,color="PCcv")) + geom_line(aes(y=PCcv, color="PCcv"),linetype="dashed") +
      ggtitle("Variance Explained by PCs") + labs(x = "Principal Component", y ="Proportion of variance explained by PCs(%)", color = "") +
      scale_color_manual(values = colors, labels=c("Cumulative Proportion","Proportion/per PC")) +
      theme(legend.position="bottom",plot.title=element_text(size=12,hjust=0.5),legend.text=element_text(size = 8), axis.title=element_text(size = 8)) +
      scale_x_continuous(breaks=seq(1,10,1))
  }
  else{
    confounderC=idTplot
    ptPCA<- ggpairs(PlotDat, columns=1:10, aes(color= PlotDat[,(confounderC+10)]), 
                    diag=list(continuous=wrap("densityDiag",alpha=0.4)),
                    lower=list(continuous = wrap("points",alpha=0.9,size=0.1)),
                    upper = list(continuous = "blank"),
                    legend = c(1,1)) + labs(fill = colnames(MetaRaw)[confounderC])
  }
  
  return(ptPCA)
}


#Check 4: Techical confounders. Input normANOV(pc_norm$x[,1:5], or exprDat_norm)
#NOTE: number of freeze-thaws looks quite confused.
ConfouderCheck <- function(normANOV,plateID,metadata_reord,bloodStain,sampleAge){
  ps_norm <- data.frame(matrix(NA,nrow=ncol(normANOV),ncol=5))
  names(ps_norm) <- c("Plate","Cohort","Group","Bloodstain","SampleAge")
  
  for (i in 1:ncol(normANOV)){
    tempMd <-  anova(lm(normANOV[,i] ~ as.factor(plateID)))
    ps_norm$Plate[i] <- tempMd$Pr[1]
    tempMd <-  anova(lm(normANOV[,i] ~ as.factor(metadata_reord$`Cohort name`)))
    ps_norm$Cohort[i] <- tempMd$Pr[1]
    tempMd <-  anova(lm(normANOV[,i] ~ as.factor(metadata_reord$Group)))
    ps_norm$Group[i] <- tempMd$Pr[1]
    tempMd <-  anova(lm(normANOV[,i] ~ as.factor(bloodStain)))
    ps_norm$Bloodstain[i] <- tempMd$Pr[1]
    tempMd <-  anova(lm(normANOV[,i] ~ sampleAge))
    ps_norm$SampleAge[i] <- tempMd$Pr[1]
  }
  
  ConfounderTable = ps_norm      ###confounders for per protein
  
  return(ConfounderTable)
}

###confounder per PC / per protein; variation explained 
confounderPlot2 = function(ConfounderTableX){ 
  toPlot = ConfounderTableX
  confounder = c("Plate","Cohort","Group","Bloodstain","SampleAge")
  
  par(mfrow=c(2,5))
  for (ccounter in 1:length(confounder)){
    temp <- p.adjust(toPlot[,ccounter],method="bonferroni")
    plot(ecdf(temp),xlab="padj value",ylab="Cumulative total",main=paste(confounder[ccounter]," effect per protein"))
    abline(v=0.05,lty=2)
    abline(h=ecdf(temp)(0.05),lty=2)
    abline(h=0.8,lty=2)
    abline(v=quantile(temp,0.8),lty=2)
    text(0.05,0.4,ecdf(temp)(0.05),pos=4,col="red")
  }
  for (ccounter in 1:length(confounder)){
    temp <- toPlot[,ccounter]
    qqplot(-log10(seq(0,1,length.out=length(temp))),-log10(temp),xlab="Expected -log10(p)",ylab="Observed -log10(p)",main=paste(confounder[ccounter]," QQplot"),cex.main=0.8,cex.lab=0.8)
    abline(0,1)
  }
  
  myplot <- recordPlot()
  return(myplot)
}


#Check 5: External validation
#check Ben (this is all OA) / Historic data (this is all knee injury)
getpars <- function(sandwich_master,exprDat_norm,par1,par2,metadata_reord) {
  
  temp1 <- sandwich_master[as.character(metadata_reord$`STEpUP Participant Identification Number (PIN)`),]
  
  if (sum(!(is.na(temp1[,par1] + exprDat_norm[,par2]))) == 0) return(c(par1,par2,NA,NA,0))
  md <- cor.test(temp1[,par1],exprDat_norm[,par2])
  c(par1,par2,unlist(md[c("estimate","p.value")]),2+md$parameter)
  
}

### toTest1,toTest2 are biomarker lists from sandwich file and adat file individually
ExtVal1 <- function(sandwich_master,exprDat_norm,toTest1,toTest2,metadata_reord){
  
  CorData_norm = data.frame(matrix(NA,nrow=length(toTest1),ncol=5))
  
  for (compCounter in 1:length(toTest1)){
    CorData_norm[compCounter,] <- getpars(sandwich_master,exprDat_norm,toTest1[compCounter],toTest2[compCounter],metadata_reord)
  }
  colnames(CorData_norm) <- c("ImmunoName","SomaName","cor","Pvalue","N")
  return(CorData_norm) 
}


getpars2 <- function(temp1,par1,par2,keep) {
  if (sum(!(is.na(temp1[keep,par1] + exprDat_norm[keep,par2]))) == 0) return(c(par1,par2,NA,NA,0))
  md <- cor.test(temp1[keep,par1],exprDat_norm[keep,par2])
  c(par1,par2,unlist(md[c("estimate","p.value")]),2+md$parameter)
}

###titleType here, for the convenience of generate title of plots, should be "sandwith_Ben" or "sandwith_Historic"
ExtVal2 <- function(sandwich_master,titleType,exprDat_norm,metadata_reord,plateID,toTest1,toTest2,conTxt){
  temp1 <- sandwich_master[as.character(metadata_reord$`STEpUP Participant Identification Number (PIN)`),]
  boxplot(sapply(unique(plateID),function(p) as.numeric(sapply(1:length(toTest1),function(i) getpars2(temp1,toTest1[i],toTest2[i],plateID == p,conTxt)["estimate.cor"]))),ylab="Correlation coefficient",main=paste(titleType," external check"))
  
  compCounterPl=1
  CorData_plate =list()
  for (p in unique(plateID)){
    x1 <- apply(temp1[plateID == p,toTest1],2,function(x) (x - mean(x,na.rm=T))/sd(x,na.rm=T))
    x2 <- apply(exprDat_norm[plateID == p,toTest2],2,function(x) (x - mean(x,na.rm=T))/sd(x,na.rm=T))
    if (sum(!is.na(x1 + x2)) == 0) next
    CorData_plate[[compCounterPl]] = c(p,cor(as.numeric(x1),as.numeric(x2),use="complete.obs"))
    compCounterPl = compCounterPl + 1
  }
  CorData_plateDone = combineList(CorData_plate)
  return(CorData_plateDone)
}

#check the HT treated data (5 OA and 5 injury)
getpars3 <- function(par1,par2,exprDat_norm) {
  temp1 <- sandwich_master_HAse_UT[as.character(metadata_reord$`STEpUP Participant Identification Number (PIN)`),]
  temp2 <- sandwich_master_HAse_HT[as.character(metadata_reord$`STEpUP Participant Identification Number (PIN)`),]
  temp3 <- sandwich_master_HAse_HTF[as.character(metadata_reord$`STEpUP Participant Identification Number (PIN)`),]
  
  md1 <- cor.test(temp1[,par1],exprDat_norm[,par2])
  md2 <- cor.test(temp2[,par1],exprDat_norm[,par2])
  md3 <- cor.test(temp3[,par1],exprDat_norm[,par2])
  temp <- c(par1,par2,md1$estimate,md1$p.value,2+md1$parameter,md2$estimate,md2$p.value,2+md2$parameter,md3$estimate,md3$p.value,2+md3$parameter)
  names(temp) <- NULL
  temp
}

HTcheck <- function(sandwich_master_HAse,exprDat_norm,toTest1,toTest2){
  
  HACorData_norm <- data.frame(matrix(NA,nrow=5,ncol=11))
  names(HACorData_norm) <- c("SandwichName","SomaName","UT_cor","UT_Pvalue","UT_N","HT_cor","HT_Pvalue","HT_N","HTF_cor","HTF_Pvalue","HTF_N")
  
  for (counter in 1:length(toTest1)){
    HACorData_norm[counter,] <- getpars3(toTest1[counter],toTest2[counter],exprDat_norm)
  }
  
  return(HACorData_norm)
}


#check against predicted R2 from the repeats
R2repeats = function(R2_norm,CorData_norm,clinicalType){
  plot(R2_norm[CorData_norm$SomaName],as.numeric(CorData_norm$cor)^2,xlab="Predicted R2",ylab="Actual R2",xlim=c(0,1.05),ylim=c(0,1.05),main=clinicalType)
  temp <- toupper(substr(CorData_norm$SandwichName,1,nchar(CorData_norm$SandwichName) - 2))
  text(R2_norm[CorData_norm$SomaName],as.numeric(CorData_norm$cor)^2,temp,cex=0.7,pos=2)
  abline(0,1)
  return()
}

### extract MetaRaw data
ExtractAdat <- function(MySoma,MetaRawM){
  exprDat_norm <- data.frame(MySoma[,!grepl("HybControlElution|NonBiotin|None",names(MySoma))])
  exprDat_norm <- as.matrix(MySoma[grep("STEP",MySoma$SampleId),-c(1:29)])
  
  MetaRaw <- as.matrix(MetaRawM[grep("STEP",MySoma$SampleId),])
  
  calib_norm <- data.frame(MySoma[,!grepl("HybControlElution|NonBiotin|None",names(MySoma))])
  calib_norm <- as.matrix(MySoma[grep("POOL",MySoma$SampleId),-c(1:29)])
  
  plasma_norm <- data.frame(MySoma[,!grepl("HybControlElution|NonBiotin|None",names(MySoma))])
  plasma_norm <- as.matrix(MySoma[grep("^[0-9]*$",MySoma$SampleId),-c(1:29)])
  
  calibIDs <-  MySoma$SampleId[grep("POOL",MySoma$SampleId)]
  calibPlates <-  MySoma$PlateId[grep("POOL",MySoma$SampleId)]
  
  AdatExtract <- vector(mode="list",length=6)
  AdatExtract[[1]] = exprDat_norm
  AdatExtract[[2]] = calib_norm
  AdatExtract[[3]] = plasma_norm
  AdatExtract[[4]] = calibIDs
  AdatExtract[[5]] = calibPlates
  AdatExtract[[6]] = MetaRaw
  
  return(AdatExtract)
}

##read the metadata
ExtractTranch <- function(inputfile3,MySoma){
  metadata_xls <- read_excel(inputfile3)
  temp1 <- data.frame(as.matrix(metadata_xls)[-1,1:17])
  names(temp1) <- as.matrix(metadata_xls)[1,1:17]
  metadata <- temp1
  rownames(metadata) <- metadata$`STEpUP Sample Identification Number (SIN)`
  
  #extract STEPUP and plate IDs
  stepupID <- MySoma$SampleId[grep("STEP",MySoma$SampleId)]
  stepupID[stepupID == "STEP1409F-V1-HT1"] <- "STEP1409-F-V1-HT1" #fix an apparent typo
  stepupID <- gsub("F-V1","V1-F",stepupID) 
  plateID <- MySoma$PlateId[grep("STEP",MySoma$SampleId)]
  
  #reorder meta-data
  metadata_reord <- metadata[stepupID,]  ### now metadat_record have the consistent "stepupID" indexing with MySoma
  
  #extract relevant measures
  bloodStain <- as.character(metadata_reord$`Grading of SF bloodstaining prior to centrifugation  (if known)`)
  bloodStain[bloodStain == "-"] <- NA
  bloodStain <- factor(bloodStain,levels=0:4)
  sampleAge <- 2020 - as.numeric(as.character(metadata_reord$`Date of biological sampling`))
  patientAge <- as.numeric(as.character(metadata_reord$`Patient age at Baseline SF sample`))
  
  TranchExtract = vector(mode="list",length=5)
  TranchExtract[[1]] = metadata_reord
  TranchExtract[[2]] = bloodStain
  TranchExtract[[3]] = sampleAge
  TranchExtract[[4]] = patientAge
  TranchExtract[[5]] = plateID
  
  return(TranchExtract) 
}


##read the Master sandwich data
ExtractSandwich = function(inputfile4){
  sandwich_master_xls_1 <- read_excel(inputfile4,sheet=1)
  sandwich_master_xls_2 <- read_excel(inputfile4,sheet=2)
  sandwich_master_xls_3 <- read_excel(inputfile4,sheet=3)
  
  #process Ben data, averaging across replicates
  temp1 <- data.frame(sandwich_master_xls_2)
  temp2 <- (temp1[temp1$replicate == 1,-c(1:2)] + temp1[temp1$replicate == 2,-c(1:2)])/2
  sandwich_master_Ben <- data.frame(PIN=temp1[temp1$replicate == 1,1],temp2)
  rownames(sandwich_master_Ben) <- sandwich_master_Ben$PIN
  
  #process Historic data
  sandwich_master_Historic <- data.frame(sandwich_master_xls_3)
  rownames(sandwich_master_Historic) <- sandwich_master_Historic$PIN
  
  SandwichExtract = vector(mode="list",length=2)
  SandwichExtract[[1]] = sandwich_master_Ben
  SandwichExtract[[2]] = sandwich_master_Historic
  
  return(SandwichExtract)
}


###unqualified RFUs below LoD
LoDdetection <- function(RawM){
  
  ###abstract plateId
  PlateIdUni = levels(factor(RawM$PlateId))
  
  DatStartId <- which(colnames(RawM)=="CLI")+1  ###for calculation convenience, extract data zone only
  
  ### limit of detection calculation LoD = mean(Buffer) + 4.9*MAD(Buffer) (median absolute deviation)
  MBuffer = RawM[which(RawM$SampleType == "Buffer"),]
  datBuffer = MBuffer[,DatStartId:ncol(MBuffer)]
  BuffMedian = apply(datBuffer,2,median)
  LoD = BuffMedian + 4.9* t(apply(apply(datBuffer,1,function(x){abs(x - BuffMedian)}),1,median))
  rownames(LoD) = "LoD"
  
  MSample = RawM[grep("Sample",RawM$SampleType),]
  DatSamp = MSample[,DatStartId:ncol(MSample)]
  
  ### which RFU is below LoD? CompM
  CompM = t(apply(DatSamp,1,function(x){x-LoD}))
  
  ### histgram of proportion below LoD per protein
  ProteinRatio = apply(CompM,2,function(x){length(which(x<0))/nrow(CompM)})
  ### histgram of proportion below LoD per sammple
  SampRatio = apply(CompM,1,function(x){length(which(x<0))/ncol(CompM)})
  
  par(mfrow=c(1,2))
  hist(-log(ProteinRatio),main="(-log)RFUs below LoD among proteins",xlab="ratio",breaks=50,cex.main=0.8)
  hist(-log(SampRatio),main="(-log)RFUs below LoD among samples",xlab="ratio",breaks=50,cex.main=0.8)
  
  myplot <- recordPlot()
  
  return(myplot)
}

### KNN batch effect test
KNNtest <- function(exprDat_norm,MetaRaw){
  ### k nearest neighborhood batch effect test
  for (roundCount in 1:1){
    pc_norm <- prcomp(log10(as.matrix(exprDat_norm)),scale = TRUE)
    distPairs = pc_norm$x[,1:10]
    DisM = as.matrix(dist(distPairs, method = "euclidean", diag = TRUE, upper = TRUE))
    
    sampSizeS = seq(40,100,by=20)  ### different percentage of 1%~25% samples for batch effect test (based on the paper)
    kNearS = seq(25,100,by=25)          ### different neiboughood sizes, based on testing, >250 no rejection at all
    positiveRate = matrix(NA,nrow=length(kNearS),ncol=length(sampSizeS))
    batchTest = vector(mode="list",length=5)
    
    
    for (batchCf in 1:5){ ###batchCf: "PlateID" "diseaseGroup" "Corhort" "bloodStain" "sampleAge"   
      
      tblWh = table(MetaRaw[,batchCf])
      EpectedRio = tblWh/sum(tblWh)
      
      for (sampSct in 1:length(sampSizeS)){
        sampSize = sampSizeS[sampSct]
        randSampId = sample(1:435,sampSize) ### random select which samples to be tested
        
        for (neiSize in 1:length(kNearS)){
          
          kNear = kNearS[neiSize]   ###set neighborhood 
          
          idInMetaR = matrix(NA,nrow=kNear-1,ncol=1)
          testSampYN = matrix(NA,nrow=sampSize,ncol=1)
          
          
          for (testP in  1:sampSize){
            
            testSamp = DisM[randSampId[testP],]  ###define test sample
            neighborD = sort(testSamp)[2:kNear]  ### find its nearest kNear neighbors
            neiborOriID = names(neighborD)       ### find its neighbor name, for batch label matching
            
            neiNameList = rownames(MetaRaw)
            
            for (neiCounter in 1:(kNear-1)){ 
              idInMetaR[neiCounter] = which(neiNameList %in% neiborOriID[neiCounter]) ### meta for tested sample
            }
            
            neiDat = MetaRaw[idInMetaR,batchCf]
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
            else{if (batchCf==5) {reportYN = fisher.test(tblTest,simulate.p.value = TRUE)$p.value} ### sample age too many categories
              else{reportYN = fisher.test(tblTest)$p.value}
            }
            
            ### chi square based multinomial test
            # chi.sq.value <- sum((tblTest[1,]/sampSize - EpectedRio)^2/EpectedRio)
            
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
    
    # boxplot(batchTest[[1]],xlab = "tested sample amount",ylab="Percentage of rejection H0 (Neighborhood size 25 ~ 250)",main = "One Round of Plate Batch Effect Test on Raw RFUs",cex.lab=0.7,cex.main=0.8,cex.axis=0.8)
    
    BatchEffectM = sapply(batchTest,mean)
    names(BatchEffectM) = c("PlateID","Disease Group","Corhort","BloodStain","SampleAge")
    BatchEffectM <- as.table(t(as.matrix(BatchEffectM)))
    rownames(BatchEffectM) = "Rejection Rate"
  }
  return(BatchEffectM)
}
