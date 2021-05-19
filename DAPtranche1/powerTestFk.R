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


### for 2 clusters
vbet.vin = vector(mode="numeric",length=9)
TestPower = vector(mode="numeric",length=9)


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
vbet.vin[k] = vbet/vin

maxCounter=100
RecCluster = vector(mode="integer",length=maxCounter)
for (testCounter in 1:maxCounter){
  fks = vector(mode="numeric",length=6)
  for (myK in 1:6){fks[myK] <-FkStatistic(X,myK)[[1]]}
  
  if(length(which(fks<0.85))==0){RecCluster[testCounter] = 1
  }else{
    RecCluster[testCounter] = which(fks<0.85)[1]
  }
}

TestPower[k] = length(which(RecCluster==2))/maxCounter
}


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

maxCounter=100
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
}


plot(vbet.vin,TestPower,type="l",col="red",main="Power Analysis (sample size =435)",xlab="Between-cluster to within-cluster variance ratio",ylab="power")
lines(vbet.vin3,TestPower3,type="l",col="blue")
abline(h = 0.8,lty=2 )
legend("bottomright", legend=c("ground truth: 2 clusters", "ground truth: 3 clusters"),
       col=c("red", "blue"), lty=1:1, cex=0.8)

pcaClus2 <- prcomp(X,scale = TRUE)
cols=vector()
cols[j] = "blue"
cols[-j]="red"

pairs(pcaClus2$x[,1:3],col=cols)


pcaClus3 <- prcomp(X,scale = TRUE)
cols=vector()
cols[j1] = "blue"
cols[j2]="red"
cols[j3]="green"
pairs(pcaClus3$x[,1:3],col=cols)

----------------------------------------------------------------------------------------------------------------------------------###different thought

### detect covariance structure of our observed data. Using original data is time consuming, we use eigen modules.
softPower = 6;
adjacency = adjacency(exprDat_norm, power = softPower)
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM
minModuleSize = 50;
tomTree = hclust(as.dist(dissTOM), method = "average");
dynamicMods = cutreeDynamic(dendro = tomTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
dynamicColors = labels2colors(dynamicMods)
MEList = moduleEigengenes(exprDat_norm, colors = dynamicColors)
MEs = MEList$eigengenes

### use top PCs
pcStr <- prcomp(log10(as.matrix(exprDat_norm)),scale = TRUE)
topPCn <- which(get_eigenvalue(pcStr)$cumulative.variance.percent>80)[1]
pcDat = pcStr$x[,1:topPCn]

covInput = MEs ### or pcDat
covInput = pcDat
### gaurantee the covariance strucutre. 
temp = cor(covInput)
ftSz = ncol(covInput)
NUM0.3 = floor(length(which(temp<0.3))/ftSz) ### 42% [0,0.3] 
NUM0.7 = floor(length(which(temp>0.7))/ftSz) ### 17% [0.7,1], 41% [0.3,0.7],
NUM0.3_0.7 = ncol(covInput) - NUM0.3 -NUM0.7

###construct covariance matrix for simulated data
SIMcov = matrix(0,nrow=ftSz,ncol=ftSz) ### covariance matrix SIMcov
for (i in 1:ftSz){
  for(j in 1:ftSz){
    if(i<j){
      randP = sample(1:ftSz,1) ### guarantee the propotion of 0.3,0.7,...
      if(randP<=NUM0.3){randV = runif(1,min=0,max=0.3)}
      else if(randP > NUM0.3 & randP <= NUM0.7) {randV = runif(1,min=0.31,max=0.7)}
      else{randV = runif(1,min=0.71,max=1)}
      SIMcov[i,j]=randV}
  }
}

for (i in 1:ftSz){
  for(j in 1:ftSz){
    if(i==j){SIMcov[i,j]=1}
    if(i>j){SIMcov[i,j]=SIMcov[j,i]}
  }}
# SIMcov is the empricially matching covariance matrix


### begin artificial test
mytestCounter = 40

RecCluster = vector(mode="integer",length=mytestCounter) ###recommended cluster number by f(k)
Nb = vector(mode="integer",length=mytestCounter) ### recommended cluster number by Nb

for(testCounter in 1:mytestCounter){
  ###construct different cluster numbers
  if(testCounter <= mytestCounter/4){
    MUSIM = rep(0,ftSz)
    SigmaSIM = diag(1,ftSz,ftSz)
  }else if(testCounter >mytestCounter/4 & testCounter <2*mytestCounter/4){
    MUSIM=rep(c(0,1),length.out=ftSz)
    SigmaSIM=nearPD(SIMcov, conv.tol = 1e-7)$mat
  }else if(testCounter >2*mytestCounter/4 & testCounter <3*mytestCounter/4){
    MUSIM=rep(c(0,1,2),length.out=ftSz)
    SigmaSIM=nearPD(SIMcov, conv.tol = 1e-7)$mat
  }else{
    MUSIM=rep(c(0,1,2,3),length.out=ftSz)
    SigmaSIM=nearPD(SIMcov, conv.tol = 1e-7)$mat
  }
  
  X = mvrnorm(nrow(exprDat_norm),mu=MUSIM,Sigma =SigmaSIM)
  # X = mvrnorm(nrow(exprDat_normX),mu=rep(c(0,10,100),each=floor(ncol(exprDat_normX)/3)),Sigma = diag(1,ncol(exprDat_normX),ncol(exprDat_normX)))
  
  NbCluster <- length(unique(NbClust(data = X ,distance = "euclidean", min.nc = 2, max.nc = 15, 
                                     method = "kmeans", index = "all", alphaBeale = 0.1)$Best.partition))
  
  fks = vector(mode="numeric",length=5)
  for (myK in 1:5){fks[myK] <-FkStatistic(X,myK)[[1]]}
  
  if(length(which(fks<0.85))==0){
    RecCluster [testCounter] = 1
    Nb[testCounter] = 1
  }else{
    RecCluster[testCounter] = which(fks<0.85)[length(which(fks<0.85))]
    Nb[testCounter] = NbCluster
  }
}


length(which(RecCluster==1))
length(which(RecCluster==2))
length(which(RecCluster==3))

length(which(Nb==1))
length(which(Nb==2))
length(which(Nb==3))
length(which(Nb==4))
