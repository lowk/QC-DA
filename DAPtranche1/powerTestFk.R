library(MASS)
library(NbClust)
library(diceR)
library(mclust)
### local functions: for F(k) calculation

getALFA <- function(thisk,Nd){
  if (thisk == 2) {temp = 1 - 3/(4*Nd)}
  else{temp <- getALFA(thisk-1,Nd) + (1-getALFA(thisk-1, Nd))/6}
  return(temp)
}

FkStatistic <- function(pcDat,thisk){
  Nd = ncol(pcDat)
  kmeanStr <- kmeans(pcDat, thisk, iter.max = 50, nstart = 1) ### here we can extract sample member ship for further comparison with "concensus clustering ect.
  ClustLabel = kmeanStr$cluster
  Sk = kmeanStr$tot.withinss
  if(thisk==1){fs=1}
  else if(FkStatistic(pcDat,thisk-1)[[2]]==0){fs=1}
  else{fs = Sk/(getALFA(thisk,Nd)*FkStatistic(pcDat,thisk-1)[[2]])}
  return(list(fs,Sk,ClustLabel))
}

### define sample size based on our case
# exprDat_normX = exprDat_norm[which(MetaRaw[,"diseaseGroup"]=="Injury"),]
# 
# pcStr <- prcomp(log10(as.matrix(exprDat_normX)),scale = TRUE)
# topPCn <- which(get_eigenvalue(pcStr)$cumulative.variance.percent>80)[1]
# ftSz = topPCn
ftSz = 10
# smpSz = 
smpSz = 754 ### injury 219; OA 754
X=matrix(0,nrow=smpSz,ncol=ftSz)
trueLabel = vector(mode="integer",length=smpSz)

### begin artificial test
### generate mean (mm) and sd(ss)

### for 1 cluster
varianceType=5
vbet.vin1 = vector(mode="numeric",length=varianceType)
TestPower11 = vector(mode="numeric",length=varianceType)
TestPower12 = vector(mode="numeric",length=varianceType)
wrong11= vector(mode="numeric",length=varianceType)
wrong12= vector(mode="numeric",length=varianceType)

maxCounter=5

ss = seq(1,varianceType)

for(k in 1:varianceType){
  
  RecCluster = vector(mode="integer",length=maxCounter)
  RecCluster2 = vector(mode="integer",length=maxCounter)
  
  
  for (testCounter in 1:maxCounter){
    
    for (i in 1:smpSz){
      X[i,]=rnorm(ftSz, mean=2, sd=ss[k])
    }
    vbet.vin1[k]=ss[k]
    
    sig_obj <- sigclust(X, k = 2, nsim = 10, labflag = 0, label = pam_2)
    if(sig_obj==0){RecCluster2[testCounter] = 1}
    else{RecCluster2[testCounter]=2}
    
    ### for each between-cluster to within-cluster variance ratio, repeated samples = maxCounter
    
    fks = vector(mode="numeric",length=6)
    for (myK in 1:6){fks[myK] <-FkStatistic(X,myK)[[1]]} ###myK are the cluster numbers applied to kmeans
    
    if(length(which(fks<0.85))==0){RecCluster[testCounter] = 1
    }else{
      RecCluster[testCounter] = which(0.85-fks == min(abs(0.85-fks[which(0.85-fks>0)])))
    }
  }
  
  TestPower11[k] = length(which(RecCluster==1))/maxCounter
  TestPower12[k] = length(which(RecCluster==1))/maxCounter
  
  wrong11[k] = length(which(RecCluster!=1))/maxCounter
  wrong12[k] = length(which(RecCluster!=1))/maxCounter
}


### for 2 clusters
### Injury: Power >1, mm=30, varianceType=50, alpha=20;  Power=k

varianceType=50
mm=30
alpha=20
maxCounter=10

vbet.vin2 = vector(mode="numeric",length=varianceType)
TestPower2R = vector(mode="numeric",length=varianceType)
TestPower2 = vector(mode="numeric",length=varianceType)
memberAccuracy = vector(mode="numeric",length=varianceType)
sig_obj = vector(mode="numeric",length=varianceType)
memberAccuracyDice = vector(mode="numeric",length=varianceType)
memberAccuracyKD = vector(mode="numeric",length=varianceType)

### generate simulated dataset X. 
for(k in 1:varianceType){
  
  RecCluster = vector(mode="integer",length=maxCounter)
  
  for (testCounter in 1:maxCounter){  
    
    j=sample(1:smpSz,floor(smpSz/2),replace=FALSE) ### setting which sample is 1 or 2
    tt=sample(1:ftSz,floor(ftSz/2),replace=FALSE)  ### setting which feature have mean + tt or mean - tt
    
    for (i in 1:smpSz){
      if(i %in% j){
        for (colI in 1:ftSz){
          if(colI %in% tt){X[i,colI]=rnorm(1, mean=0.1*colI, sd=alpha*k)}
          else{X[i,colI]=rnorm(1, mean=-0.1*colI, sd=alpha*k)}
        }
        trueLabel[i]=1
      }
      else{
        for (colI in 1:ftSz){
          if(colI %in% tt){X[i,colI]=rnorm(1, mean=mm+ alpha*colI, sd=alpha*k)}
          else{X[i,colI]=rnorm(1, mean=mm-0.1*colI, sd=alpha*k)}
        }
        trueLabel[i]=2
      }
    }
    
    cen1 = apply(X[j,],2,mean)
    cen2 = apply(X[-j,],2,mean)
    vbet = sqrt(apply(as.matrix((cen1-cen2)^2),2,sum))/2
    vin = (sum(sqrt(apply((apply(X[j,],1,function(x){(x-cen1)^2})),2,sum))) + 
             sum(sqrt(apply((apply(X[-j,],1,function(x){(x-cen2)^2})),2,sum))))/nrow(X)
    vbet.vin2[k] = vbet/vin
    
    
    fks = vector(mode="numeric",length=6)
    for (myK in 1:6){fks[myK] <-FkStatistic(X,myK)[[1]]}
    
    if(length(which(fks<0.85))==0){RecCluster[testCounter] = 1}
    else{
      RecCluster[testCounter] = which(fks == min(fks))
    }
  }
  
  TestPower2[k] = length(which(RecCluster>1))/maxCounter   ### power >1
  TestPower2R[k] = length(which(RecCluster==2))/maxCounter ### power =k
  memberAccuracy[k] = mclust::adjustedRandIndex(trueLabel,FkStatistic(X,3)[[3]])  ### kmeans membership accuracy
  tempClust <- consensus_cluster(X, nk = 2, p.item = 0.8, reps = 5,algorithms = c("hc","km","hdbscan"))
  # CC <- apply(tempClust, 2:4, impute_knn, data = X, seed = 1)
  # CC_imputed <- impute_missing(CC, X, nk = 2)
  # sig_obj[k] <- sigclust(X, k = 2, nsim = 100, labflag = 0, label = CC_imputed)
  memberDiceR <- apply(tempClust,1,function(x){names(which.max(table(x)))})  
  memberAccuracyDice[k] = mclust::adjustedRandIndex(trueLabel,memberDiceR)   ### concensus membership accuracy
  memberAccuracyKD[k] = mclust::adjustedRandIndex(memberDiceR,FkStatistic(X,3)[[3]]) ### membership consistency between kmeans and concensus
}

plot(sort(vbet.vin2),TestPower2[order(vbet.vin2)],type="b",col="red",main="Power to detect K>1, sample size = Injury group sample size",xlab="Between-cluster to within-cluster variance ratio",ylab="power")

cols=vector()
cols[j] = "blue"
cols[-j]="red"
pairs(X[,1:5],col=cols)

temp = prcomp(X,scale = TRUE)$x
pairs(temp[,1:5],col=cols)


### for 3 clusters. 
### Injury: Power >1, mm=c(0,1,2),varianceType=100,ss=1,alpha=0.01; Power=k

mm=c(0,10,20)
varianceType=50
ss=1
alpha=0.8
maxCounter=10
# pdf("cluster3.pdf")
# par(mfrow=c(1,3))
vbet.vin3 = vector(mode="numeric",length=varianceType)
TestPower3R = vector(mode="numeric",length=varianceType)
TestPower3R2 = vector(mode="numeric",length=varianceType)
TestPower3R3 = vector(mode="numeric",length=varianceType)
TestPower3 = vector(mode="numeric",length=varianceType)
memberAccuracy = vector(mode="numeric",length=varianceType)
memberAccuracy2 = vector(mode="numeric",length=varianceType)
memberAccuracyDice = vector(mode="numeric",length=varianceType)
memberAccuracyKD = vector(mode="numeric",length=varianceType)

j1 = vector(mode="integer")
j2= vector(mode="integer")
j3= vector(mode="integer")
jc1=1
jc2=1
jc3=1

### equal sd, for TestPower3
for(k in 35:varianceType){
  
  RecCluster = vector(mode="integer",length=maxCounter)
  RecCluster2 = vector(mode="integer",length=maxCounter)
  RecCluster3 = vector(mode="integer",length=maxCounter)
  
  for (testCounter in 1:maxCounter){
    
    for (i in 1:smpSz){
      if(i<floor(smpSz/3)){
        tt = sample(1:ftSz,floor(ftSz/2),replace=FALSE)
        for(ttt in 1:ftSz){
          if(ttt %in% tt){X[i,ttt]=rnorm(1, mean=mm[1]+0.1*ttt, sd=ss+alpha*k)}
          else{X[i,ttt]=rnorm(1, mean=mm[1]-0.1*ttt, sd=ss+alpha*k)}
        }
        trueLabel[i]=1
        j1[jc1]=i
        jc1=jc1+1}
      else if(i>2*smpSz/3){
        tt = sample(1:ftSz,floor(ftSz/2),replace=FALSE)
        for(ttt in 1:ftSz){
          if(ttt %in% tt){X[i,ttt]=rnorm(1, mean=mm[2]+0.1*ttt, sd=ss+alpha*k)}
          else{X[i,ttt]=rnorm(1, mean=mm[2]-0.1*ttt, sd=ss+alpha*k)}
        }
        trueLabel[i]=3
        j2[jc2]=i
        jc2=jc2+1}
      else{
        tt = sample(1:ftSz,floor(ftSz/2),replace=FALSE)
        for(ttt in 1:ftSz){
          if(ttt %in% tt){X[i,ttt]=rnorm(1, mean=mm[3]+0.1*ttt, sd=ss+alpha*k)}
          else{X[i,ttt]=rnorm(1, mean=mm[3]-0.1*ttt, sd=ss+alpha*k)}
        }
        trueLabel[i]=2
        j3[jc3]=i
        jc3=jc3+1}
    }
    
    
    cen1 = apply(X[j1,],2,mean)
    cen2 = apply(X[j2,],2,mean)
    cen3 = apply(X[j3,],2,mean)
    vbet = (sqrt(apply(as.matrix((cen1-cen2)^2),2,sum))+sqrt(apply(as.matrix((cen2-cen3)^2),2,sum))+
              sqrt(apply(as.matrix((cen1-cen3)^2),2,sum)))/3
    vin = (sum(sqrt(apply((apply(X[j1,],1,function(x){(x-cen1)^2})),2,sum))) + 
             sum(sqrt(apply((apply(X[j2,],1,function(x){(x-cen2)^2})),2,sum)))+
             sum(sqrt(apply((apply(X[j3,],1,function(x){(x-cen3)^2})),2,sum))))/nrow(X)
    vbet.vin3[k] = vbet/vin
    
    
    fks = vector(mode="numeric",length=6)
    for (myK in 1:6){fks[myK] <-FkStatistic(X,myK)[[1]]}
    
    if(length(which(fks<0.85))==0){
      RecCluster[testCounter] = 1
      RecCluster2[testCounter]=1
    }else{
      RecCluster[testCounter] = which(fks==min(fks))
      tempP <- NbClust(data = X ,distance = "euclidean", min.nc = 2, max.nc = 15, method = "kmeans", index = "all", alphaBeale = 0.1)$Best.partition
      RecCluster2[testCounter] <- length(unique(tempP))
      # tempClust <- consensus_cluster(X, nk = 2:4, p.item = 0.8, reps = 5,algorithms = c("hc","km","hdbscan"))
      # RecCluster3[testCounter]  <- consensus_evaluate(X,tempClust,plot = FALSE)$k
    }
  }
  
  TestPower3[k] = length(which(RecCluster>1))/maxCounter
  TestPower3R[k] = length(which(RecCluster==3))/maxCounter
  TestPower3R2[k] = length(which(RecCluster2==3))/maxCounter
  # TestPower3R3[k] = length(which(RecCluster3==3))/maxCounter
  memberAccuracy[k] = mclust::adjustedRandIndex(trueLabel,FkStatistic(X,3)[[3]])
  memberAccuracy2[k] = mclust::adjustedRandIndex(trueLabel,tempP)
  # tempClust <- consensus_cluster(X, nk = 3, p.item = 0.8, reps = 5,algorithms = c("hc","km","hdbscan"))
  # CC <- apply(tempClust, 2:4, impute_knn, data = X, seed = 1)
  # CC_imputed <- impute_missing(CC, X, nk = 3)
  # sig_obj[k] <- sigclust(X, k = 3, nsim = 100, labflag = 0, label = CC_imputed)
  # memberDiceR <- apply(tempClust,1,function(x){names(which.max(table(x)))})
  # memberAccuracyDice[k] = mclust::adjustedRandIndex(trueLabel,memberDiceR)
  # memberAccuracyKD[k] = mclust::adjustedRandIndex(memberDiceR,FkStatistic(X,3)[[3]])
}
plot(sort(vbet.vin3),TestPower3[order(vbet.vin3)],type="b",lty=2,col="blue",main="Power to detect K>1, sample size = Injury group sample size",xlab="Between-cluster to within-cluster variance ratio",ylab="power")
plot(sort(vbet.vin3),TestPower3R[order(vbet.vin3)],type="b",lty=2,col="blue",main="f(K), Power to detect correct 3 clusters (equal sd wiwthin cluster,sample size = OA sample size",xlab="Between-cluster to within-cluster variance ratio",ylab="power")
plot(sort(vbet.vin3),TestPower3R2[order(vbet.vin3)],type="b",lty=2,col="blue",main="NbClust, Power to detect correct 3 clusters (equal sd wiwthin cluster,sample size = OA sample size",xlab="Between-cluster to within-cluster variance ratio",ylab="power")

# dev.off()

cols=vector()
cols[j1] = "blue"
cols[j2]="red"
cols[j3]="green"
pairs(X[,1:5],col=cols)
temp = prcomp(X,scale = TRUE)$x
pairs(temp[,1:5],col=cols)

plot(sort(vbet.vin2),TestPower2[order(vbet.vin2)],type="b",col="red",main="Power to detect K>1, sample size = OA group sample size",xlab="Between-cluster to within-cluster variance ratio",ylab="power")
lines(sort(vbet.vin3),TestPower3[order(vbet.vin3)],type="b",lty=2,col="blue")
abline(h = 0.8,lty=2)
legend("bottomright", legend=c("ground truth:2 clusters","ground truth:3 clusters"),col=c("red", "blue","green"), lty=1:1, cex=0.8)

plot(sort(vbet.vin2),TestPower2R[order(vbet.vin2)],type="b",col="red",main="Power to detect correct K, sample size = OA group sample size",xlab="Between-cluster to within-cluster variance ratio",ylab="power")
lines(sort(vbet.vin3),TestPower3R[order(vbet.vin3)],type="b",lty=2,col="blue")
lines(sort(vbet.vin3),TestPower3R2[order(vbet.vin3)],type="b",lty=2,col="green")
abline(h = 0.8,lty=2)
legend("bottomright", legend=c("ground truth:2 clusters","ground truth:3 clusters,f(K)","ground truth:3 clusters,NbClust"),col=c("red", "blue","green"), lty=1:1, cex=0.8)

plot(sort(vbet.vin2),memberAccuracy[order(vbet.vin2)],type="b",lty=2,col="red",main="Membership Accuracy \n 2 clusters (equal sd within cluster,sample size = Injury sample size",xlab="Between-cluster to within-cluster variance ratio",ylab="Adjusted Rand Index")
lines(sort(vbet.vin2),memberAccuracyDice[order(vbet.vin2)],type="b",lty=2,col="blue")
lines(sort(vbet.vin2),memberAccuracyKD[order(vbet.vin2)],type="b",lty=2,col="green")
legend("bottomright", legend=c("kmeans vs ground truth","diceR vs ground truth","kmeans vs diceR"),col=c("blue", "red","green"), lty=1:1, cex=0.8)


plot(sort(vbet.vin2),memberAccuracy[order(vbet.vin2)],type="b",lty=2,col="red",main="Membership Accuracy, sample size = OA sample size",xlab="Between-cluster to within-cluster variance ratio",ylab="Adjusted Rand Index",xlim=c(0,1))
lines(sort(vbet.vin3),memberAccuracy[order(vbet.vin3)],type="b",lty=2,col="blue")
abline(h = 0.8,lty=2)
legend("bottomright", legend=c("ground truth:2 clusters","ground truth:3 clusters"),col=c("blue", "red"), lty=1:1, cex=0.8)


plot(sort(vbet.vin3),memberAccuracy[order(vbet.vin3)],type="b",lty=2,col="red",main="Membership Accuracy \n 3 clusters (equal sd within cluster,sample size = Injury sample size",xlab="Between-cluster to within-cluster variance ratio",ylab="Adjusted Rand Index")
lines(sort(vbet.vin3),memberAccuracyDice[order(vbet.vin3)],type="b",lty=2,col="blue")
lines(sort(vbet.vin3),memberAccuracyKD[order(vbet.vin3)],type="b",lty=2,col="green")
legend("bottomright", legend=c("kmeans vs ground truth","diceR vs ground truth","kmeans vs diceR"),col=c("blue", "red","green"), lty=1:1, cex=0.8)

plot(sort(vbet.vin3),TestPower3R2[order(vbet.vin3)],type="l",lty=2,col="red",main="Power to detect correct 3 clusters (equal sd wiwthin cluster,sample size = Injury sample size",xlab="Between-cluster to within-cluster variance ratio",ylab="power")
lines(sort(vbet.vin3),TestPower3R[order(vbet.vin3)],type="l",lty=2,col="blue")
lines(sort(vbet.vin3),TestPower3R3[order(vbet.vin3)],type="l",lty=2,col="green")
legend("bottomright", legend=c("NbClust","f(K)","diceR"),col=c("red", "blue","green"), lty=1:1, cex=0.8)


