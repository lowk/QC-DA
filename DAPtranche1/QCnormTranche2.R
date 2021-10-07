###07 OCT, somewhere changed compared to github record(which is correct one to generate released data), need further check.
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
#     DatStartId = which(colnames(RawMS1) == "RMA") +1
#     DatStartIdP = which(colnames(RawMS) == "RMA") +1
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


### initialise required metadata for normalisation methods
initQCnorm <- function(inputfile1,inputfile2){
  
  RawM <- read.adat(inputfile1)
  ### read in Col^MetaTable
  ###use inputfile1, there is no PlateScale_Reference information, so we use inputfile2 considering all the row and column indices are consistent
  
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
      if(slineL[getStartTerm]=="medNormRefSMP_ReferenceRFU"){
        PlateScale_Reference <<- as.numeric(slineL[getStartId:length(slineL)])  
        # load("PlateScale_Reference.Rdata")
        break} ### tranche2 data, Platesclae_Reference is labelled "medNormRefSMP_ReferenceRFU"
      
      ### 30 is the beginning field of the ^COL_DATA;column names of meta data = SeqId +Target
    }
  }
  close(con1)
  
  ### pay attention, tranche2 RawM protein data begin from the 26th col
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
  
  return(list(RawM,ColTable))
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
  ###01 JAN 2021 = 44197 as our uptodate
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
  
  # DatStartId <- which(colnames(RawM)=="CRYBB2.10000.28")
  # DiseaseType <- levels(as.factor(GroupTRawM))
  # for (Dcounter in 1:length(DiseaseType)){
  #   RawM$SampleType[which(GroupTRawM==DiseaseType[Dcounter])]=paste(DiseaseType[Dcounter],"Sample",sep="")
  #   
  # }
  
  ### up to now integrate disease group type into RawM, RawM$SampleType including Pool type, disease sample type
  stepupIDTRaw[stepupIDTRaw == "STEP1409F-V1-HT1"] <- "STEP1409-F-V1-HT1"
  
  BioMeta = cbind(RawM$PlateId,stepupIDTRaw,GroupTRawM,CohortTRawM,bloodStainTRawM,sampleAgeTRaw,patientGenderTRaw,patientAgeTRaw)
  colnames(BioMeta)[1:2] = c("PlateID","STEpUPID")
  rownames(BioMeta) = rownames(RawM)
  
  return(BioMeta) 
}