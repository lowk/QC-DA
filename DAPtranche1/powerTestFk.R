library(MASS)

getALFA <- function(thisk,Nd){
  if (thisk == 2) {temp = 1 - 3/(4*Nd)}
  else{temp <- getALFA(thisk-1,Nd) + (1-getALFA(thisk-1, Nd))/6}
  return(temp)
}

FkStatistic <- function(pcDat,thisk){
  Nd = ncol(pcDat)
  kmeanStr <- kmeans(pcDat, thisk, iter.max = 50, nstart = 1) 
  Sk = kmeanStr$tot.withinss
  if(thisk==1){fs=1}
  else if(FkStatistic(pcDat,thisk-1)[[2]]==0){fs=1}
  else{fs = Sk/(getALFA(thisk,Nd)*FkStatistic(pcDat,thisk-1)[[2]])}
  return(list(fs,Sk))
}

pcStr <- prcomp(log10(as.matrix(exprDat_normX)),scale = TRUE)
topPCn <- which(get_eigenvalue(pcStr)$cumulative.variance.percent>80)[1]
ftSz = topPCn
smpSz = nrow(exprDat_normX)

### begin artificial test
### generate mean (mm) and sd(ss)
mm = vector(mode="numeric",length=10) 
ss = vector(mode="numeric",length=10)
mm[1]=0
ss[1]=1
for (i in 2:10){
  mm[i]=mm[1]+0.2*(i-1)
  ss[i]=ss[1]-0.1*(i-1)
}

### for 1 cluster
vbet.vin1 = vector(mode="numeric",length=9)
TestPower1 = vector(mode="numeric",length=9)
wrong1= vector(mode="numeric",length=9)

X=matrix(NA,nrow=smpSz,ncol=topPCn)
for(k in 1:10){
  for (i in 1:smpSz){
    X[i,]=rnorm(ftSz, mean=0, sd=ss[k])
  }
  vbet.vin1[k]=ss[k]

maxCounter=10
RecCluster = vector(mode="integer",length=maxCounter)
for (testCounter in 1:maxCounter){
  fks = vector(mode="numeric",length=6)
  for (myK in 1:6){fks[myK] <-FkStatistic(X,myK)[[1]]}
  
  if(length(which(fks<0.85))==0){RecCluster[testCounter] = 1
  }else{
    RecCluster[testCounter] = which(fks<0.85)[1]
  }
}

TestPower1[k] = length(which(RecCluster==1))/maxCounter
wrong1[k] = length(which(RecCluster!=1))/maxCounter
}


### for 2 clusters
vbet.vin2 = vector(mode="numeric",length=9)
TestPower2 = vector(mode="numeric",length=9)
wrong2 = vector(mode="numeric",length=9)

X=matrix(NA,nrow=smpSz,ncol=topPCn)
for(k in 1:9)
{j=sample(1:smpSz,floor(smpSz/2),replace=FALSE)
for (i in 1:smpSz){
  if(i %in% j){X[i,]=rnorm(ftSz, mean=mm[1], sd=ss[1])}
  else(X[i,]=rnorm(ftSz, mean=mm[k+1], sd=ss[k+1]))
}

cen1 = apply(X[j,],2,mean)
cen2 = apply(X[-j,],2,mean)
vbet = sqrt(apply(as.matrix((cen1-cen2)^2),2,sum))/2
vin = (sum(sqrt(apply((apply(X[j,],1,function(x){(x-cen1)^2})),2,sum))) + 
         sum(sqrt(apply((apply(X[-j,],1,function(x){(x-cen2)^2})),2,sum))))/nrow(X)
vbet.vin2[k] = vbet/vin

maxCounter=10
RecCluster = vector(mode="integer",length=maxCounter)
for (testCounter in 1:maxCounter){
  fks = vector(mode="numeric",length=6)
  for (myK in 1:6){fks[myK] <-FkStatistic(X,myK)[[1]]}
  
  if(length(which(fks<0.85))==0){RecCluster[testCounter] = 1
  }else{
    RecCluster[testCounter] = which(fks<0.85)[1]
  }
}

TestPower2[k] = length(which(RecCluster==2))/maxCounter
wrong2[k] = length(which(RecCluster==1))/maxCounter
}
pairs(X[,1:3])

### for 3 clusters
mm = vector(mode="numeric",length=10) 
ss = vector(mode="numeric",length=10)
mm[1]=0
ss[1]=1
for (i in 2:10){
  mm[i]=mm[1]+0.3*(i-1)
  ss[i]=ss[1]-0.1*(i-1)
}


vbet.vin3 = vector(mode="numeric",length=9)
TestPower3 = vector(mode="numeric",length=9)
wrong3 = vector(mode="numeric",length=9)
j1 = vector(mode="integer")
j2= vector(mode="integer")
j3= vector(mode="integer")
jc1=1
jc2=1
jc3=1
for(k in 1:7)
{for (i in 1:smpSz){
  if(i<floor(smpSz/3)){X[i,]=rnorm(ftSz, mean=mm[1]+k*0.2, sd=ss[9])
  j1[jc1]=i
  jc1=jc1+1}
  else if(i>2*smpSz/3){X[i,]=rnorm(ftSz, mean=mm[2]+k*0.2, sd=ss[9])
  j2[jc2]=i
  jc2=jc2+1}
  else{X[i,]=rnorm(ftSz, mean=mm[3]+k*0.2, sd=ss[9])
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

maxCounter=10
RecCluster = vector(mode="integer",length=maxCounter)

for (testCounter in 1:maxCounter){
  fks = vector(mode="numeric",length=6)
  for (myK in 1:6){fks[myK] <-FkStatistic(X,myK)[[1]]}
  
  if(length(which(fks<0.85))==0){RecCluster[testCounter] = 1
  }else{
    RecCluster[testCounter] = which(fks<0.85)[length(which(fks<0.85))]
  }
}

TestPower3[k] = length(which(RecCluster==3))/maxCounter
wrong3[k] = length(which(RecCluster==1))/maxCounter
}


plot(vbet.vin2[1:9],TestPower2[1:9],type="l",col="red",main="Power Analysis (sample size =)",xlab="Between-cluster to within-cluster variance ratio",ylab="power")
lines(vbet.vin3,TestPower3,type="l",col="blue")
lines(vbet.vin1,TestPower1,type="l",col="green")
abline(h = 0.8,lty=2 )
legend("bottomright", legend=c("ground truth: 1 clusters", "ground truth: 2 clusters","ground truth: 3 clusters"),
       col=c("red", "blue","green"), lty=1:1, cex=0.8)

cols=vector()
cols[j] = "blue"
cols[-j]="red"

pairs(pcaClus2$x[,1:3],col=cols)


cols=vector()
cols[j1] = "blue"
cols[j2]="red"
cols[j3]="green"
pairs(pcaClus3$x[,1:3],col=cols)
