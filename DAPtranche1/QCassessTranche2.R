##1: Total protein check. Input normalised RFUs dataframe, output plots

TotalProCheck <- function(exprDat_norm,BioMeta){
  totalProtein_norm <- apply(exprDat_norm,1,sum)
  hist(totalProtein_norm,breaks=100,xlab="Total Protein",main="Total Protein Distribution",cex.main=2,cex.lab=1.5,cex.axis=1.5)
  plot(as.numeric(BioMeta$bloodStain),totalProtein_norm,xlab="Blood Staining Grade",ylab="Total Protein",cex.lab=1.5,cex.axis=1.5)
  abline(lm(totalProtein_norm ~ as.numeric(BioMeta$bloodStain)),col="red",lwd=2)
  text(1,max(totalProtein_norm)*0.95,paste("p =",formatC(cor.test(totalProtein_norm, as.numeric(BioMeta$bloodStain))$p.val,format="e",digits=3)),pos=4,col="red",cex=2)
  text(1,max(totalProtein_norm)*0.9,paste("cor =",formatC(cor.test(totalProtein_norm, as.numeric(BioMeta$bloodStain))$estimate,format="e",digits=3)),pos=4,col="red",cex=2)
  # boxplot(totalProtein_norm ~ BioMeta$diseaseGroup,main=paste("p =",
  #                                                             signif(wilcox.test(log(totalProtein_norm) ~ BioMeta$diseaseGroup)$p.val,digits=4)),
  #                                                             xlab="Disease Group",ylab="Total Protein Within Group",cex.main=2,cex.lab=1.5,cex.axis=1.5)
  return(totalProtein_norm)
}
# which(totalProtein_norm<500000000)
# text(10000,250,"tranche1",pos=4,col="red",cex=2)
# text(1500000000,100,"tranche2",pos=4,col="red",cex=2)

##2: Checks against calibrators. Input selected RFUs calib_norm, and clinicalType ("OA" or "INJ")
# pool %CV
CalibratorCheck <- function(calib_norm,clinicType){
  calibIDs <-  calib_normM$SampleId  ###calibID corresponding to calib_norm
  temp1 <- apply(calib_norm[grep(clinicType,calibIDs),],2,sd)/apply(calib_norm[grep(clinicType,calibIDs),],2,mean)
  plot(100*quantile(temp1,seq(0,1,length.out=100)),seq(0,1,length.out=100),type="l",xlim=c(0,100),xlab=paste("%CV",clinicType),ylab="Cumulative total",cex.lab=1.5,main="Repeatiblity of Injury sample replicates")
  abline(h=0.8,lty=2)
  abline(v=100*quantile(temp1,0.8),lty=2)
  text(100*quantile(temp1,0.8)-4,0.2,signif(quantile(temp1,0.8),4),pos=4,col="red",cex=1.5)
  
  # ###pool MAD  
  # tempMatrix = calib_norm[grep(clinicType,calibIDs),]
  # temp2 <- apply(apply(tempMatrix,1,function(x){abs(x-apply(tempMatrix,2,median))}),2,median)
  # plot(quantile(temp2,seq(0,1,length.out=100)),seq(0,1,length.out=100),type="l",xlab=paste("MAD",clinicType),ylab="Cumulative total",main="Norm")
  # abline(h=0.8,lty=2)
  # abline(v=quantile(temp2,0.8),lty=2)
  
  return()
}

# pool variance explained
VarExp <- function(calib_normM,calib_norm,clinicType,exprDat_norm){
  calibIDs <- calib_normM$SampleId
  temp <- apply(calib_norm[grep(clinicType,calibIDs),],2,var)/apply(exprDat_norm,2,var)
  temp[temp > 1] <- 1
  R2_norm <- (1 - temp)^2
  plot(quantile(R2_norm,1 - seq(0,1,length.out=100)),seq(0,1,length.out=100),xlim=c(1,0),type="l",xlab=paste("R2",clinicType),ylab="Cumulative total",main="Norm")
  abline(h=0.8,lty=2)
  abline(v=quantile(R2_norm,1 - 0.8),lty=2)
  return(R2_norm)
}


##CV breakdowns:
CVbreak <- function(calib_norm,clinicType){
  if(clinicType == "OA"){suffix = "/29"
  freezeThaw <-  calibIDs %in% paste0("OA POOL-HT-",c(6,25,26),suffix)
  }else{suffix = "/25"
  freezeThaw <-  calibIDs %in% paste0("INJ POOL-HT-",c(6,25),suffix)}
  
  acrossPlates <- calibIDs %in% paste0(clinicType," POOL-HT-",c(1:5,7:13),suffix)
  withinPlates <- calibIDs %in% paste0(clinicType," POOL-HT-",c(1,6),suffix)
  
  UNSP <- grepl("UNSP",calibIDs)
  
  temp1 <- (apply(calib_norm[acrossPlates,],2,sd)/apply(calib_norm[acrossPlates,],2,mean))
  temp2 <- (apply(calib_norm[withinPlates,],2,sd)/apply(calib_norm[withinPlates,],2,mean))
  temp3 <- (apply(calib_norm[freezeThaw,],2,sd)/apply(calib_norm[freezeThaw,],2,mean))
  temp4 <- (apply(calib_norm[UNSP,],2,sd)/apply(calib_norm[UNSP,],2,mean))
  
  plot(100*quantile(temp1,seq(0,1,length.out=100)),seq(0,1,length.out=100),type="l",xlim=c(0,10),xlab=paste("%CV",clinicType,"Group"),ylab="Cumulative total",main="",lwd=2,cex.lab=1.5)
  lines(100*quantile(temp2,seq(0,1,length.out=100)),seq(0,1,length.out=100),type="l",lwd=2,col="red")
  lines(100*quantile(temp3,seq(0,1,length.out=100)),seq(0,1,length.out=100),type="l",lwd=2,col="blue")
  lines(100*quantile(temp4,seq(0,1,length.out=100)),seq(0,1,length.out=100),type="l",lwd=2,col="green")
  
  lines(c(0,quantile(temp1,0.8))*100,c(0.8,0.8),lty=2)
  lines(c(quantile(temp1,0.8),quantile(temp1,0.8))*100,c(0.8,0),lty=2)
  lines(c(0,quantile(temp2,0.8))*100,c(0.8,0.8),lty=2,col="red")
  lines(c(quantile(temp2,0.8),quantile(temp2,0.8))*100,c(0.8,0),lty=2,col="red")
  lines(c(0,quantile(temp3,0.8))*100,c(0.8,0.8),lty=2,col="blue")
  lines(c(quantile(temp3,0.8),quantile(temp3,0.8))*100,c(0.8,0),lty=2,col="blue")
  lines(c(0,quantile(temp4,0.8))*100,c(0.8,0.8),lty=2,col="green")
  lines(c(quantile(temp4,0.8),quantile(temp4,0.8))*100,c(0.8,0),lty=2,col="green")
  
  legend(8,0.4,c("Across","Within","freezeThaw","Unspun"),lwd=2,col=c("black","red","blue","green"),xpd=T,cex=0.75)

  ##each of the five plates against the average of the rest
  plot(NA,xlim=c(0,20),ylim=c(0,1),xlab=paste("%CV",clinicType,"Group"),ylab="Cumulative total",main="Per Plate Against the Others",cex.lab=1.5)
  printList=vector(mode="list",length=length(acrossPlates))
  mycol=c("black","red","blue","green","yellow","purple","cyan","brown","pink","darkgoldenrod1","darkgreen","deeppink")
  for (i in 1:12){
    temp4 <- abs(calib_norm[which(acrossPlates)[i],] - apply(calib_norm[which(acrossPlates)[-i],],2,mean))/apply(calib_norm[acrossPlates,],2,mean)
    lines(100*quantile(temp4,seq(0,1,length.out=100)),seq(0,1,length.out=100),type="l",lwd=2,col=mycol[i])
    printList[[i]]=(c(calibPlates[which(acrossPlates)[i]],mean(temp4)))
  }
  printListDone = combineList(printList)
  legend(15,0.5,col=mycol,legend=calibPlates[acrossPlates],lwd=2,xpd=T,cex=0.75)
  
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

  return(printListDone)
}


#Check 3: PCA
PCAglob <- function(BioMeta,exprDat_norm){
  pc_norm <- prcomp(as.matrix(exprDat_norm),scale = TRUE)
  #myTSNE <-Rtsne(temp, dims = 2, perplexity=30, verbose=TRUE, max_iter = 500)
  
  # for (plotID in 1:ncol(MetaRaw)){
  #   plot(data.frame(temp),col = factor(MetaRaw[,plotID]), pch=19, cex=0.08)
  #   legend("right",inset=c(-0.04,-1),unique(MetaRaw[,plotID]),lwd=2,col=1:length(unique(MetaRaw[,plotID])),xpd=T,cex=0.7,title=colnames(MetaRaw)[plotID])
  # }
  # 
  # for (plotID in 1:ncol(MetaRaw)){
  #   plot(myTSNE$Y, col=factor(MetaRaw[,plotID]), main="tSNE", xlab="tSNE D1", ylab="tSNE D2", pch=19, cex=0.2)
  #   legend("bottomright",unique(MetaRaw[,plotID]),lwd=2,col=1:length(unique(MetaRaw[,plotID])),xpd=T,cex=0.4,title=colnames(MetaRaw)[plotID])
  # }
  # 
  
  PlotDat = data.frame(cbind(pc_norm$x[,1:10],BioMeta))
  
  for (confounderC in 1:6){
     
    ggpairs(PlotDat, columns=1:10, aes(color= as.factor(PlotDat[,confounderC+10])),
            diag=list(continuous=wrap("densityDiag",alpha=0.4)),
            lower=list(continuous = wrap("points",alpha=0.9,size=0.1)),
            upper = list(continuous = "blank"),
            legend = c(1,1)) + labs(fill = colnames(BioMeta)[confounderC])
    
    PlotDatFil = PlotDat[!is.na(PlotDat$diseaseGroup),]
    PlotDatFil = PlotDat[!is.na(PlotDat$Corhort),]
    PlotDatFil = PlotDat[!is.na(PlotDat$bloodStain)&PlotDat$bloodStain!="Not known"&PlotDat$bloodStain!="No"&PlotDat$bloodStain!="Mild",]
    
    ggpairs(PlotDatFil, columns=1:10, aes(color= as.factor(PlotDatFil[,confounderC+10])),
            diag=list(continuous=wrap("densityDiag",alpha=0.4)),
            lower=list(continuous = wrap("points",alpha=0.9,size=0.1)),
            upper = list(continuous = "blank"),
            legend = c(1,1)) + labs(fill = colnames(BioMeta)[confounderC])
    
    PlotDatFil = PlotDat[!is.na(PlotDat$sampleAge),]
    ggpairs(PlotDatFil, columns=1:10, aes(color= as.factor(floor(as.numeric(PlotDatFil[,confounderC+10])/365))),
            diag=list(continuous=wrap("densityDiag",alpha=0.4)),
            lower=list(continuous = wrap("points",alpha=0.9,size=0.1)),
            upper = list(continuous = "blank"),
            legend = c(1,1)) + labs(fill = colnames(BioMeta)[confounderC])
    
  }
  
   ###extract most important PCs for PPT presentation: "PlateID" PC3 VS PC8; "diseaseGroup" PC2 VS PC6;  "Corhort" PC2 VS PC6; "bloodStain": PC6 PC9; "sampleAge" PC2 PC6 
  ggplot(PlotDat,aes(x=PC5,y=PC6,color=tranche)) + geom_point(size=3) + labs(color="tranche") + theme(legend.title=element_text(size=15),legend.text=element_text(size=16),axis.title=element_text(size=13,face="bold"))
  
  PlotDat = PlotDat[!is.na(PlotDat$PlateID),]
  ggplot(PlotDat,aes(x=PC5,y=PC9,color=PlateID)) + geom_point(size=3) + labs(color="PlateID") + theme(legend.title=element_text(size=15),legend.text=element_text(size=16),axis.title=element_text(size=13,face="bold"))
  
  PlotDat = PlotDat[!is.na(PlotDat$diseaseGroup),]
  ggplot(PlotDat,aes(x=PC2,y=PC7,color=diseaseGroup))+geom_point(size=3) + labs(color="Dissease Group") + theme(legend.title=element_text(size=18),legend.text=element_text(size=18),axis.title=element_text(size=18,face="bold"))
  
  PlotDat = PlotDat[!is.na(PlotDat$Corhort),]
  ggplot(PlotDat,aes(x=PC2,y=PC8,color=Corhort))+geom_point(size=3) + labs(color="Corhort")+ theme(legend.title=element_text(size=18),legend.text=element_text(size=18),axis.title=element_text(size=18,face="bold"))
  
  PlotDat = PlotDat[!is.na(PlotDat$bloodStain)&PlotDat$bloodStain!="Not known"&PlotDat$bloodStain!="No"&PlotDat$bloodStain!="Mild",]
  ggplot(PlotDat,aes(x=PC1,y=PC7,color=as.numeric(bloodStain)))+geom_point(size=3) + labs(color="Blood Stain") + theme(legend.title=element_text(size=18),legend.text=element_text(size=18),axis.title=element_text(size=18,face="bold"))
  
  PlotDat = PlotDat[!is.na(PlotDat$sampleAge),]
  ggplot(PlotDat,aes(x=PC2,y=PC8,color=as.numeric(sampleAge)/365))+geom_point(size=3) +  scale_color_gradient(breaks=seq(0,20,5),limits=c(0,20)) + labs(color="Sample Age") + theme(legend.title=element_text(size=15),legend.text=element_text(size=16),axis.title=element_text(size=13,face="bold"))
  
  
  ### PCA outliers based on top 10 PCs
  myPoints <- pc_norm$x[,1:10]
  centroid <- colMeans(myPoints)
  distanceToCenter = matrix(NA,ncol=1,nrow=nrow(myPoints))
  for (mK in 1: nrow(myPoints)){
    distanceToCenter[mK] <- sqrt(sum((myPoints[mK,]-centroid)^2))
  }
  colnames(distanceToCenter) <- "Distance To Center"
  hist(distanceToCenter, main="Histogtam of Sample Distances to Center Based on Top 10 PCs",xlab="Distance Between Samples and Center")
  length(which(distanceToCenter>(mean(distanceToCenter) + 2*sd(distanceToCenter))))
  
  distanceToCenter2 <- data.frame(distanceToCenter)
  ggplot(distanceToCenter2) + geom_histogram(aes(x=Distance.To.Center,y=..density..),color="lightblue",binwidth=1) + 
    geom_density(aes(x=Distance.To.Center),alpha=.2, fill="#FF6666") + labs(title="Histogtam of Sample Distances to Center Based on Top 10 PCs") +
    xlab("Distance Between Samples and Center") + ylab("Density") + theme(plot.title=element_text(size=15,hjust=0.5))
  
  eig.val <- get_eigenvalue(pc_norm) 
  egV = eig.val$eigenvalue[1:10] ###eigenvalues 
  PCv = eig.val$variance.percent[1:10]     ### corresponding PC variance percentage 
  PCcv = eig.val$cumulative.variance.percent[1:10]   ### corresponding cumulative PCv
  PCcount = seq(1:10)
  
  dat = as.data.frame(cbind(PCcount,egV,PCv,PCcv))
  p1 = ggplot(data=dat, aes(x=PCcount,y=egV)) + geom_point(color="cyan") + geom_line(color="cyan") +
    xlab("Principal Component") + ylab("Eigenvalue") + ggtitle("Scree Plot")
  
  colors <- c("PCv"="cyan","PCcv"="darkcyan")
  p2 =ggplot(dat,aes(x=PCcount)) + geom_point(aes(y=PCv,color="PCv")) + geom_line(aes(y=PCv, color="PCv")) +
    geom_point(aes(y=PCcv,color="PCcv")) + geom_line(aes(y=PCcv, color="PCcv"),linetype="dashed") +
    ggtitle("Variance Explained by PCs") + labs(x = "Principal Component", y ="Proportion of variance explained by PCs(%)", color = "") +
    scale_color_manual(values = colors, labels=c("Cumulative Proportion","Proportion/per PC")) + 
    theme(legend.position="bottom",plot.title=element_text(size=18,hjust=0.5),legend.text=element_text(size = 15), axis.title=element_text(size = 13)) + 
    scale_x_continuous(breaks=seq(1,10,1))
  
  listPCA = list()
  listPCA[[1]] = pc_norm ###first 10 PCs
  listPCA[[2]] = p1  ### plot of Eigenvalue against PCs
  listPCA[[3]] = p2  ### plot of variations explained per PC
  
  return(listPCA)
}


#Check 4: Techical confounders. Input normANOV(pc_norm$x[,1:5], or exprDat_norm)
#NOTE: number of freeze-thaws looks quite confused.
ConfouderCheck <- function(totalProtein_norm,normANOV,BioMeta){
  TPA1 = anova(lm(totalProtein_norm ~ as.factor(BioMeta$PlateID)))$Pr[1]
  TPA2 = anova(lm(totalProtein_norm ~ as.factor(BioMeta$Corhort)))$Pr[1]
  # TPA3 = anova(lm(totalProtein_norm ~ as.factor(BioMeta$diseaseGroup)))$Pr[1]
  TPA4 = anova(lm(totalProtein_norm ~ as.factor(BioMeta$bloodStain)))$Pr[1]
  TPA5 = anova(lm(totalProtein_norm ~ BioMeta$sampleAge))$Pr[1]
  TPA6 = anova(lm(totalProtein_norm ~ BioMeta$tranche))$Pr[1]
  # TPA = matrix(c(TPA1,TPA2,TPA3,TPA4,TPA5,TPA6),nrow=1)
  # colnames(TPA) = c("Plate","Cohort","Group","Bloodstain","SampleAge","tranche")
  
  TPA = matrix(c(TPA1,TPA2,TPA4,TPA5,TPA6),nrow=1)
  colnames(TPA) = c("Plate","Cohort","Bloodstain","SampleAge","tranche")
  
  ps_norm <- data.frame(matrix(NA,nrow=ncol(normANOV),ncol=length(TPA)))
  # names(ps_norm) <- c("Plate","Cohort","Group","Bloodstain","SampleAge","tranche")
  names(ps_norm) <- c("Plate","Cohort","Bloodstain","SampleAge","tranche")
  
  varExp_norm<- ps_norm
  
  for (i in 1:ncol(normANOV)){
    tempMd <-  anova(lm(normANOV[,i] ~ as.factor(BioMeta$PlateID)))
    ps_norm$Plate[i] <- tempMd$Pr[1]
    varExp_norm$Plate[i] <- tempMd$Sum[1]/sum(tempMd$Sum)
    
    tempMd <-  anova(lm(normANOV[,i] ~ as.factor(BioMeta$Corhort)))
    ps_norm$Cohort[i] <- tempMd$Pr[1]
    varExp_norm$Cohort[i] <- tempMd$Sum[1]/sum(tempMd$Sum)
    
    # tempMd <-  anova(lm(normANOV[,i] ~ as.factor(BioMeta$diseaseGroup)))
    # ps_norm$Group[i] <- tempMd$Pr[1]
    # varExp_norm$Group[i] <- tempMd$Sum[1]/sum(tempMd$Sum)
    
    tempMd <-  anova(lm(normANOV[,i] ~ as.factor(BioMeta$bloodStain)))
    ps_norm$Bloodstain[i] <- tempMd$Pr[1]
    varExp_norm$Bloodstain[i] <- tempMd$Sum[1]/sum(tempMd$Sum)
    
    tempMd <-  anova(lm(normANOV[,i] ~ BioMeta$sampleAge))
    ps_norm$SampleAge[i] <- tempMd$Pr[1]
    varExp_norm$SampleAge[i] <- tempMd$Sum[1]/sum(tempMd$Sum)
    
    tempMd <-  anova(lm(normANOV[,i] ~ BioMeta$tranche))
    ps_norm$tranche[i] <- tempMd$Pr[1]
    varExp_norm$tranche[i] <- tempMd$Sum[1]/sum(tempMd$Sum)
  }
  
  ConfounderTable = list()
  ConfounderTable[[1]]=TPA          ###confounders for total protein 
  ConfounderTable[[2]]=ps_norm      ###confounders for first 5 principal components
  ConfounderTable[[3]]=varExp_norm  ### variances explained based on per confounder
  
  return(ConfounderTable)
}


###confounder total proteins
confounderPlot1 = function(ConfounderTable){ 
  toPlot = ConfounderTable[[1]]
  plot(1:5,p.adjust(toPlot[1,],method="bonferroni"),main="confounders total protein",xaxt="n",xlab="", ylab="padj values")
  abline(h=0.05,lty=2)
  xlabels = colnames(toPlot)
  axis(side = 1, at=seq(1:5),labels=FALSE)
  text(seq(1:5),par("usr")[3]-0.002,labels=xlabels, srt=45, adj=1, xpd=TRUE)
  return()
}

###confounder per PC / per protein; variation explained 
confounderPlot2 = function(ConfounderTableX,onWhat){ 
  if (onWhat==1){titleAdd = "PerProtein"}
  else{titleAdd = "PCs"}
  
  for (TLcounter in 2:2){
    toPlot = ConfounderTableX[[TLcounter]]
    # confounder = c("Plate","Cohort","Group","Bloodstain","SampleAge","Tranche")
    confounder = c("Plate","Cohort","Bloodstain","SampleAge","Tranche")
    # titleType = c(titleAdd,"ExplainedVariance")
    
    # for (ccounter in 1:length(confounder)){
    #   temp <- p.adjust(toPlot[,ccounter],method="bonferroni")
    #   plot(ecdf(temp),xlab="padj value",ylab="Cumulative total",main=paste(confounder[ccounter]," ",titleType[TLcounter-1]))
    #   abline(v=0.05,lty=2)
    #   abline(h=ecdf(temp)(0.05),lty=2)
    #   abline(h=0.8,lty=2)
    #   abline(v=quantile(temp,0.8),lty=2)
    #   text(0.05,0.4,ecdf(temp)(0.05),pos=4,col="red")
    # }
    for (ccounter in 1:length(confounder)){
      temp <- toPlot[,ccounter]
      qqplot(-log10(seq(0,1,length.out=length(temp))),-log10(temp),ylim = c(0, 100),xlab="Expected p Value",ylab="Observed p Value",main=paste(confounder[ccounter]," QQplot"),cex.main=1.5,cex.lab=1.5)
      abline(0,1)
    }
  }
  return()
}


#Check 5: External validation

#check Ben (this is all OA) / Historic data (this is all knee injury)
getpars <- function(sandwich_master,exprDat_norm,par1,par2,conTxt) {
  
  temp1 <- sandwich_master[as.character(metadata_reord$`STEpUP Participant Identification Number (PIN)`),]
  
  if (sum(!(is.na(temp1[,par1] + exprDat_norm[,par2]))) == 0) return(c(par1,par2,NA,NA,0))
  md <- cor.test(temp1[,par1],exprDat_norm[,par2])
  c(par1,par2,unlist(md[c("estimate","p.value")]),2+md$parameter)
  
}

### toTest1,toTest2 are biomarker lists from sandwich file and adat file individually
ExtVal1 <- function(sandwich_master,exprDat_norm,toTest1,toTest2,conTxt){
  
  CorData_norm = data.frame(matrix(NA,nrow=length(toTest1),ncol=5))
  names(CorData_norm) <- c("SandwichName","SomaName","cor","Pvalue","N")
  
  for (compCounter in 1:length(toTest1)){
    CorData_norm[compCounter,] <- getpars(sandwich_master,exprDat_norm,toTest1[compCounter],toTest2[compCounter],conTxt)
  }
  return(CorData_norm) 
}


getpars2 <- function(temp1,par1,par2,keep,conTxt) {
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
  plot(R2_norm[CorData_norm$SomaName],as.numeric(CorData_norm$cor)^2,xlab="Predicted R2",ylab="Actual R2",xlim=c(0,1.05),ylim=c(0,1.05),main=paste(clinicalType,"Group"),cex.main=2,cex.lab=1.5)
  temp <- toupper(substr(CorData_norm$SandwichName,1,nchar(CorData_norm$SandwichName) - 2))
  text(R2_norm[CorData_norm$SomaName],as.numeric(CorData_norm$cor)^2,temp,cex=0.9,pos=2)
  abline(0,1)
  return()
}

##Check 7: Blood staining markers
BloodMarker <- function(exprDat_norm){
  BMarker = c("HBB.17137.160","CAT.3488.64","PRDX1.3855.56","CA1.4969.2")
  CommanName = c("Hemoglobin (HBB)","Catalase (CAT)","Peroxiredoxin (PRDX1)","Carbonic anhydrase 1 (CA1)")
  BioMeta <- data.frame(BioMeta)
  BMarkerCor = vector(mode="list",length=length(BMarker))
  par(mfrow=c(1,3))
  selBloodId = which(!is.na(BioMeta$bloodStain)&BioMeta$bloodStain!="Not known"&BioMeta$bloodStain!="No"&BioMeta$bloodStain!="Mild")
  exprDat_normSel = exprDat_norm[selBloodId,]
  BioMetaSel = BioMeta[selBloodId,]
    for (counter in 1:length(BMarker)){
    CorStr1 = cor.test(exprDat_norm[selBloodId,BMarker[counter]],as.numeric(blood))
    boxplot(exprDat_norm[selBloodId,BMarker[counter]] ~ BioMeta$bloodStain[selBloodId],xlab = "blood stain",ylab=paste(CommanName[counter]," Quantity"),main=paste("Correlation Coefficient =", signif(CorStr1$estimate,digit=4)),cex.main=1.5,cex.lab=1.5)
    CorStrOA = cor.test(exprDat_normSel[BioMetaSel$diseaseGroup=="OA",BMarker[counter]],as.numeric(BioMetaSel$bloodStain[BioMetaSel$diseaseGroup=="OA"]))
    boxplot(exprDat_normSel[BioMetaSel$diseaseGroup=="OA",BMarker[counter]] ~ BioMetaSel$bloodStain[BioMetaSel$diseaseGroup=="OA"],xlab = "blood stain", ylab=paste(CommanName[counter]," Quantity"), main=paste("OA group, Correlation Coefficient =", signif(CorStrOA$estimate,digit=4)),cex.main=1.5,cex.lab=1.5)
    CorStrINJ =cor.test(exprDat_normSel[BioMetaSel$diseaseGroup=="Injury",BMarker[counter]],as.numeric(BioMetaSel$bloodStain[BioMetaSel$diseaseGroup=="Injury"]))
    boxplot(exprDat_normSel[BioMetaSel$diseaseGroup=="Injury",BMarker[counter]] ~ BioMetaSel$bloodStain[BioMetaSel$diseaseGroup=="Injury"],xlab = "blood stain", ylab=paste(CommanName[counter]," Quantity"),main=paste("Injury group, Correlation Coefficient =", signif(CorStrINJ$estimate,digit=4)),cex.main=1.5,cex.lab=1.5)
    BMarkerCor[[counter]] = c(CorStr1$estimate,CorSt1r$p.value,CorStrOA$estimate,CorStOAr$p.value,CorStrINJ$estimate,CorSt1INJ$p.value)
  }
  BMarkerCorC = combineList(BMarkerCor)
  BMarkerCorDone = cbind(as.matrix(BMarker),as.matrix(CommanName),BMarkerCorC)
  colnames(BMarkerCorDone) = c("BloodMarker","CommanName","CorTest_Cor","CorTest_Cor_Ttest_p")
  
  return(BMarkerCorDone)
}

##Check 8: Sex markers

data_summary <- function(x) { ### to generate mean/sd for box/violin plot
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

getSepcForSens <- function(s,y,x){
  temp <- sapply(s,function(s1) mean(y[!x] <= quantile(y[x],s1)))
  temp[s == 0] <- 0
  temp[s == 1] <- 1
  return(temp)
}

SexCheck1 <- function(exprDat_norm,metadata_reord){
  
  Bmarker = c("KLK3.8468.19","CGA.LHB.2953.31","CGA.FSHB.3032.11","CGA.CGB.4914.10")
  CommonName = c("log PSA (KLK3)","log luteinizing hormone (CGA/LHB)","log FSH (GCA/LHB)","log HCG (CGA/CGB)")

  idSel = which((!is.na(BioMeta$patientAge)) & (!is.na(BioMeta$patientGender)))
  tempAge = BioMeta$patientAge[idSel]

  for (replaceCounter in 1:length(tempAge)){
    if(tempAge[replaceCounter]<=30){tempAge[replaceCounter]="AgeStage1"
    next}
    if(tempAge[replaceCounter]>30 && tempAge[replaceCounter]<=50) {tempAge[replaceCounter]="AgeStage2"
    next}
    if(tempAge[replaceCounter]>50 && tempAge[replaceCounter]<=70) {tempAge[replaceCounter]="AgeStage3"
    next}
    if(tempAge[replaceCounter]>70) {tempAge[replaceCounter]="AgeStage4"
    next}
  }
  
  for (counter in 1:length(Bmarker)){
    # boxplot(log(exprDat_norm[,Bmarker[counter]]) ~ metadata_reord$`Patient gender`,main=CommonName[counter])
    plotDat = data.frame(cbind(exprDat_norm[idSel,Bmarker[counter]],BioMeta$patientGender[idSel],tempAge))
    colnames(plotDat) = c("SexMarker","Gender","PatientAge")
    
    ### difference between gender
    tResultAll = t.test(exprDat_norm[BioMeta$patientGender=="F",Bmarker[counter]], exprDat_norm[BioMeta$patientGender=="M",Bmarker[counter]], alternative = "two.sided")
    ggplot(data=plotDat,aes(x=Gender,y=log(as.numeric(SexMarker)),color=Gender)) +  geom_violin() + geom_boxplot(width=0.1) + 
      labs(y=CommonName[counter],color="Gender") + ggtitle(paste("t test between genders p = ", signif(tResultAll$p.value,digits=4))) +
      theme(plot.title=element_text(size=25,hjust=0.5),legend.title=element_text(size=20),legend.text=element_text(size = 20), axis.title=element_text(size = 20))
    
    ### difference among age groups
    anovResult = anova(lm(plotDat$SexMarker ~ as.factor(plotDat$PatientAge)))
    ggplot(data=plotDat,aes(x=PatientAge,y=log(as.numeric(SexMarker)),color=PatientAge)) + geom_violin() + geom_boxplot(width=0.1) +
      labs(y=CommonName[counter],color="Patient Age") + ggtitle(paste("Anova test among age groups p = ", signif(anovResult$Pr,digits=4))) + scale_colour_discrete(labels = c("16~30", "31~50","51~70","71~91")) +
      theme(plot.title=element_text(size=25,hjust=0.5),legend.title=element_text(size=20),legend.text=element_text(size = 20), axis.title=element_text(size = 20)) + scale_fill_discrete(labels = c("16~30", "31~50", "51~70","71~91"))
    
    ### difference both age and gender
    tResult=vector(mode="numeric",length=4)
    for(kk in 1:4){
      tResult[kk] = signif(t.test(as.numeric(plotDat[plotDat$PatientAge==paste("AgeStage",kk,sep="") & plotDat$Gender=="F",]$SexMarker), as.numeric(plotDat[plotDat$PatientAge==paste("AgeStage",kk,sep="") & plotDat$Gender=="M",]$SexMarker), alternative = "two.sided")$p.value,4)
    }
    
    ggplot(data=plotDat,aes(x=Gender,y=log(as.numeric(SexMarker)),fill=PatientAge)) + geom_violin(width=1) + stat_summary(fun.data=data_summary,position = position_dodge(1)) +
      labs(y=CommonName[counter],fill="Patient Age: test p value") + ggtitle(paste("t-test of the same age group between genders")) + theme(plot.title = element_text(hjust=0.5)) + 
      scale_fill_discrete(labels = c(paste("16~30:",tResult[1]), paste("31~50:",tResult[2]), paste("51~70:",tResult[3]),paste("71~91:",tResult[4]))) +
      theme(plot.title=element_text(size=25,hjust=0.5),legend.title=element_text(size=20),legend.text=element_text(size = 20), axis.title=element_text(size = 20))
    
  }
  
  
  predictors <- Bmarker
  dat_norm <- data.frame(y=BioMeta[idSel,]$patientGender == "M",exprDat_norm[idSel,predictors])
  names(dat_norm) <- c("y","x1","x2","x3","x4")
  
  ypred_norm_logit <- c()
  ypred_norm_svm <- c()
  pre_norm_list = list()
  
  for (i in 1:length(idSel)){
    ypred_norm_logit[i] <- predict(glm(y ~ x1 + x2 + x3 + x4,data=dat_norm[-i,],family="binomial"),newdata=dat_norm[i,]) 
    ypred_norm_svm[i] <- -attr(predict(svm(y ~ x1 + x2 + x3 + x4,data=dat_norm[-i,],type="C"),newdata=dat_norm[i,],decision.values = T),"decision.values")
  }
  pre_norm_list[[1]] = ypred_norm_logit
  pre_norm_list[[2]] = ypred_norm_svm
  
  return(pre_norm_list)
}

SexCheck2 <- function(pre_norm_list){
  
  ypred_norm_logitP = pre_norm_list[[1]]
  ypred_norm_svmP = pre_norm_list[[2]]
  
  ypred_norm_logit = 1/(1+exp(-ypred_norm_logitP))  ###u=1/(1+exp(-y)), ypred_norm_logitP is the "y", we need u for further mismatch check 
  ypred_norm_svm = 1/(1+exp(-ypred_norm_svmP))
  
  tempAge = as.numeric(BioMeta[idSel,]$patientAge)
  sens <- seq(0,1,length.out=1000)
  keep <- tempAge>0
  
  y <- BioMeta[idSel,]$patientGender[keep] == "M"
  
  spec_norm_logit <- getSepcForSens(sens,ypred_norm_logit[keep],y) ### true positive female ; sens: false positive female
  spec_norm_svm <-  getSepcForSens(sens,ypred_norm_svm[keep],y)
  
  auc_norm_logit <- integrate(getSepcForSens,0,1,y=ypred_norm_logit[keep],x=y,subdivisions=2000)$val
  auc_norm_svm <- integrate(getSepcForSens,0,1,y=ypred_norm_svm[keep],x=y,subdivisions=2000)$val
  
  plot(sens,spec_norm_logit,type="l",xlim=c(0,1),ylim=c(0,1),xlab="1 - Specificity",ylab="Sensitivity",lty=3,lwd=2,main="Sex Biomarker Validation for All Age Groups",col="blue",cex.main=1.5,cex.lab=1.5)
  lines(sens,spec_norm_svm,type="l",lty=1,lwd=2,col="red")
  abline(0,1)
  legend("bottomright",c(paste0("logit AUC=",signif(auc_norm_logit,3)),paste0("svm AUC=",signif(auc_norm_svm,3))),lty=c(3,3,1,1),lwd=2,col=c("black","red"),cex=1.5)
  points(sens[which(spec_norm_logit>0.9)][1],spec_norm_logit[which(spec_norm_logit>0.9)][1],pch=11,col="green",cex=2)
  text(0.08,0.95,"threshold",col = "green",cex=1.5)
  
  # MyThreshhold = sens[which(spec_norm_logit>0.9)[1]]
  MyThreshhold = 0.5
  MID = which(ypred_norm_logit >= MyThreshhold)
  MyPredict = vector(length=length(ypred_norm_logit))
  MyPredict[MID] = 1
  MyPredict[!MID] = 0 
  GroundTruth = BioMeta[idSel,]$patientGender[keep]
  GroundTruth1 <- sub("M",1,GroundTruth)
  GroundTruth2 <- sub("F",0,GroundTruth1)
  GroundTruth <- as.numeric(GroundTruth2)
  mismatchID = which(MyPredict-GroundTruth!=0)
  mismatchID50 = mismatchID[which(mismatchID %in% which(tempAge>50))]
  metadata_reord$`STEpUP Participant Identification Number (PIN)`[mismatchID50] ###gender mismatch ID
  length(mismatchID50)
  ### plot distribution of prediction score for different genders  
  #   sexDataF = data.frame(Gender.In.Record=metadata_reord$`Patient gender`,pred=ypred_norm_logit)
  #   pointsHere=data.frame(x=ypred_norm_logit[mismatchID50],y=rep(0,length(mismatchID50)),TrueGender=metadata_reord$`Patient gender`[mismatchID50],SampleID=metadata_reord$`STEpUP Sample Identification Number (SIN)`[mismatchID50])
  #   
  #   ggplot(data=sexDataF) + geom_histogram(aes(x=pred,y=..density..,color=Gender.In.Record,fill=Gender.In.Record),alpha=0.2,binwidth=0.05) + 
  #     geom_density(aes(x=pred,color=Gender.In.Record)) + ggtitle("Distribution of Gender Prediction Score based on Regression Model: 0->Female, 1->Male") + 
  #     xlab("Prediction Score") + ylab("Density") + theme(plot.title = element_text(hjust=0.5)) +
  #     geom_point(data=pointsHere,aes(x=x,y=y,color=TrueGender))
  #   
  #   FinalGenderMismatchTable = data.frame(PredictionScore=sort(pointsHere$x),GenderInRecord=pointsHere$TrueGender[order(pointsHere$x)],ZSampleId=pointsHere$SampleID[order(pointsHere$x)])
  # ### gender mismatch score distribution plot use another threshold 0.42 where 16 cases reported.
  
  for (replaceCounter in 1:length(tempAge)){
    if(tempAge[replaceCounter]<=30){tempAge[replaceCounter]="AgeStage1"
    next}
    if(tempAge[replaceCounter]>30 && tempAge[replaceCounter]<=50) {tempAge[replaceCounter]="AgeStage2"
    next}
    if(tempAge[replaceCounter]>50 && tempAge[replaceCounter]<=70) {tempAge[replaceCounter]="AgeStage3"
    next}
    if(tempAge[replaceCounter]>70) {tempAge[replaceCounter]="AgeStage4"
    next}
  }
  keepCase = c("AgeStage1","AgeStage2","AgeStage3","AgeStage4")
  ageStage = c("16~30","31~50","51~70","71~91")
  
  for(keepCounter in 1:length(ageStage)){
    keep <- which(tempAge==keepCase[keepCounter])
    y <- metadata_reord$`Patient gender`[keep] == "M"
    spec_norm_logit <- sapply(sens,function(s) getSepcForSens(s,ypred_norm_logit[keep],y))
    spec_norm_svm <- sapply(sens,function(s) getSepcForSens(s,ypred_norm_svm[keep],y))
    
    auc_norm_logit <- integrate(getSepcForSens,0,1,y=ypred_norm_logit[keep],x=y,subdivisions=2000)$val
    auc_norm_svm <- integrate(getSepcForSens,0,1,y=ypred_norm_svm[keep],x=y,subdivisions=2000)$val
    
    plot(sens,spec_norm_logit,type="l",xlim=c(0,1),ylim=c(0,1),xlab="1 - Specificity",ylab="Sensitivity",lty=3,lwd=2,main=paste("Sex Biomarker Prediction Model Validation for Age Group of ",ageStage[keepCounter]),col="blue")
    lines(sens,spec_norm_svm,type="l",lty=1,lwd=2,col="red")
    abline(0,1)
    legend("bottomright",c(paste0("logit AUC=",signif(auc_norm_logit,3)),paste0("svm AUC=",signif(auc_norm_svm,3))),lty=c(3,3,1,1),lwd=2,col=c("black","red"),cex=1.5)
  }
  
  return()
}

ExtractAdat <- function(MySoma,MetaRawM){
  # exprDat_norm <- data.frame(MySoma[,!grepl("HybControlElution|NonBiotin|None",colnames(MySoma))])  ###HybControlElution|NonBiotin|None in total, 276
  # exprDat_norm <- data.frame(MySoma[,!grepl("HybControlElution|Non",colnames(MySoma))])  ###HybControlElution|NonBiotin|None in total, 276
  # exprDat_norm <- as.matrix(exprDat_norm[grep("STEP",MySoma$SampleId),-c(1:24)]) ### up to now exprDat different rows compared to MySoma
  
  # MetaRawClinic <- as.matrix(MetaRawM[grep("STEP",MySoma$SampleId),])
  
  # calib_norm <- data.frame(MySoma[,!grepl("HybControlElution|NonBiotin|None",names(MySoma))])
  # calib_norm <- data.frame(MySoma[,!grepl("HybControlElution|Non",names(MySoma))])
  # calib_norm <- as.matrix(calib_norm[grep("POOL",MySoma$SampleId),-c(1:24)])
  
  # plasma_norm <- data.frame(MySoma[,!grepl("HybControlElution|NonBiotin|None",names(MySoma))])
  plasma_norm <- data.frame(MySoma[,!grepl("HybControlElution|Non",names(MySoma))])
  plasma_norm <- as.matrix(plasma_norm[grep("^[0-9]*$",MySoma$SampleId),-c(1:24)])
  
  calibIDs <-  MySoma$SampleId[grep("POOL",MySoma$SampleId)]  ###calibID corresponding to calib_norm
  calibPlates <-  MySoma$PlateId[grep("POOL",MySoma$SampleId)]
  
  AdatExtract <- vector(mode="list",length=6)
  AdatExtract[[1]] = exprDat_norm
  AdatExtract[[2]] = calib_norm
  AdatExtract[[3]] = plasma_norm
  AdatExtract[[4]] = calibIDs
  AdatExtract[[5]] = calibPlates
  AdatExtract[[6]] = MetaRawClinic ### MetaRawClinic only include clinical sample meta information
  
  return(AdatExtract)
}


##read the metadata
ExtractTranch <- function(inputfile3,MySoma){
  metadata <- read_excel(inputfile3,range="A2:Q438")
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
  sandwich_master_xls_4 <- read_excel(inputfile4,sheet=4,range = "A1:D31", col_names = TRUE)
  
  #process Ben data, averaging across replicates
  temp1 <- data.frame(sandwich_master_xls_2)
  temp2 <- (temp1[temp1$replicate == 1,-c(1:2)] + temp1[temp1$replicate == 2,-c(1:2)])/2
  sandwich_master_Ben <- data.frame(PIN=temp1[temp1$replicate == 1,1],temp2)
  rownames(sandwich_master_Ben) <- sandwich_master_Ben$PIN
  
  #process Historic data
  sandwich_master_Historic <- data.frame(sandwich_master_xls_3)
  rownames(sandwich_master_Historic) <- sandwich_master_Historic$PIN
  
  #process HAse data (splitting by HAse treatment)
  temp3 <- data.frame(sandwich_master_xls_4)
  sandwich_master_HAse_UT <- temp3[temp3$HAse.treatment == "UT",c(1,3,4)]
  sandwich_master_HAse_HT <- temp3[temp3$HAse.treatment == "HT",c(1,3,4)]
  sandwich_master_HAse_HTF <- temp3[temp3$HAse.treatment == "HTF",c(1,3,4)]
  rownames(sandwich_master_HAse_UT) <- sandwich_master_HAse_UT$PIN
  rownames(sandwich_master_HAse_HT) <- sandwich_master_HAse_HT$PIN
  rownames(sandwich_master_HAse_HTF) <- sandwich_master_HAse_HTF$PIN
  
  SandwichExtract = vector(mode="list",length=5)
  SandwichExtract[[1]] = sandwich_master_Ben
  SandwichExtract[[2]] = sandwich_master_Historic
  SandwichExtract[[3]] = sandwich_master_HAse_UT
  SandwichExtract[[4]] = sandwich_master_HAse_HT
  SandwichExtract[[5]] = sandwich_master_HAse_HTF
  
  return(SandwichExtract)
}


###unqualified RFUs below LoD
LoDdetection <- function(RawM,Exp){
  
  
  DatStartId <- which(colnames(RawM)=="CRYBB2.10000.28")  ###for calculation convenience, extract data zone only
  
  ### limit of detection calculation LoD = mean(Buffer) + 4.9*MAD(Buffer) (median absolute deviation)
  MBuffer = RawM[which(RawM$SampleType == "Buffer"),]
  datBuffer = MBuffer[,DatStartId:ncol(MBuffer)]
  BuffMedian = apply(datBuffer,2,median)
  LoD = as.matrix(BuffMedian + 4.9* t(apply(apply(datBuffer,1,function(x){abs(x - BuffMedian)}),1,median)))
  LoD = LoD[!grepl("HybControlElution|Non",colnames(LoD))]
  
  ### which RFU is below LoD? CompM. for combatted Exp, LoD need have log transform
  CompM = t(apply(Exp,1,function(x){x-log(LoD)}))
  
  ### histgram of proportion below LoD per protein
  ProteinRatio = apply(CompM,2,function(x){length(which(x<0))/nrow(CompM)})
  ### histgram of proportion below LoD per sammple
  SampRatio = apply(CompM,1,function(x){length(which(x<0))/ncol(CompM)})
  
  hist(ProteinRatio,main="RFUs below LoD among proteins",xlab="ratio",breaks=30,cex.main=1,cex.lab=1.5)
  hist(SampRatio,main="RFUs below LoD among samples",xlab="ratio",breaks=30,cex.main=1,cex.lab=1.5)
  
  return(CompM)
}
###??? LoB should be used for deleting single RFUs, LoD vs 50% quantile should be used to screen proteins less than 90% power

### combine print into table for print out in txt
combineList <- function(printList){
  for (lscount in 1:length(printList)){
    if(lscount==1) {printListDone=printList[[lscount]]}
    else{printListDoneT = printList[[lscount]]
    printListDone = rbind(printListDone,printListDoneT)}
  }
  return(printListDone)
}

### use pre-saved txt files, to plot sinble graph directly compare two MySoma
plotDirectComp <- function(ps_norm1,ps_norm2){
  qqplot(-log10(seq(0,1,length.out=length(ps_norm1$SampleAge))),-log10(ps_norm2$SampleAge),xlab="Expected -log10(p)",ylab="Observed -log10(p)",main="SampleAge")
  abline(0,1)
  qqplot(-log10(seq(0,1,length.out=length(ps_norm1$SampleAge))),-log10(ps_norm2$Bloodstain),xlab="Expected -log10(p)",ylab="Observed -log10(p)",main="Bloodstain")
  abline(0,1)
  qqplot(-log10(seq(0,1,length.out=length(ps_norm1$SampleAge))),-log10(ps_norm2$Cohort),xlab="Expected -log10(p)",ylab="Observed -log10(p)",main="Cohort")
  abline(0,1)
  qqplot(-log10(seq(0,1,length.out=length(ps_norm1$SampleAge))),-log10(ps_norm2$Plate),xlab="Expected -log10(p)",ylab="Observed -log10(p)",main="Plate")
  abline(0,1)
}

### per protein, confounders, distributions per group, 2nd,3rd,4th moment
# lvcount = levels(factor(MetaRaw[,"PlateID"]))
# 
# skMark = matrix(NA,nrow=length(lvcount),ncol=ncol(exprDat_norm))
# ktsMark = matrix(NA,nrow=length(lvcount),ncol=ncol(exprDat_norm))
# sdMark = matrix(NA,nrow=length(lvcount),ncol=ncol(exprDat_norm))
# 
# for (lvcountER in 1:length(lvcount) ){
#   idPer = which(MetaRaw[,"PlateID"]==lvcount[lvcountER]) 
#   skMark[lvcountER,] = apply(exprDat_norm[idPer,],2,skewness,na.rm=TRUE)
#   colnames(skMark) = colnames(exprDat_norm)
#   ktsMark[lvcountER,] =  apply(exprDat_norm[idPer,],2,kurtosis,na.rm=TRUE)
#   colnames(ktsMark) = colnames(exprDat_norm)
#   sdMark[lvcountER,] = apply(exprDat_norm[idPer,],2,sd,na.rm=TRUE)
#   colnames(sdMark) = colnames(exprDat_norm)
# }

KNNtest <- function(exprDat_norm,MetaRaw){
  ### k nearest neighborhood batch effect test
  for (roundCount in 1:5){
    pc_norm <- prcomp(as.matrix(exprDat_norm),scale = TRUE)
    distPairs = pc_norm$x[,1:10]
    DisM = as.matrix(dist(distPairs, method = "euclidean", diag = TRUE, upper = TRUE))
    
    sampSizeS = seq(10,500,by=100)  ### different percentage of 1%~25% samples for batch effect test (based on the paper)
    kNearS = seq(10,500,by=100)          ### different neiboughood sizes, based on testing, >250 no rejection at all
    positiveRate = matrix(NA,nrow=length(kNearS),ncol=length(sampSizeS))
    batchTest = vector(mode="list",length=6)
    
    
    for (batchCf in 1:6){ ###batchCf: "PlateID" "diseaseGroup" "Corhort" "bloodStain" "sampleAge"   
      
      tblWh = table(MetaRaw[,batchCf])
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
            else{reportYN = fisher.test(tblTest,simulate.p.value = TRUE)$p.value}
            
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
    
    BatchEffectM[[roundCount]] = sapply(batchTest,mean)
    # names(BatchEffectM) = c("tranche","PlateID","Disease Group","Corhort","BloodStain","SampleAge")
    # BatchEffectM <- as.table(t(as.matrix(BatchEffectM)))
    # rownames(BatchEffectM) = "Rejection Rate"
  }
  return(BatchEffectM)
}
