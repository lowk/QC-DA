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

### detect covariance structure of our observed data. Using original data is time consuming, we use eigen modules.
softPower = 6;
adjacency = adjacency(exprDat_normX, power = softPower)
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

# ### use top PCs
# pcStr <- prcomp(log10(as.matrix(exprDat_norm)),scale = TRUE)
# topPCn <- which(get_eigenvalue(pcStr)$cumulative.variance.percent>80)[1]
# pcDat = pcStr$x[,1:topPCn]

### gaurantee the covariance strucutre. 
temp = cor(MEs)
ftSz = ncol(MEs)
NUM0.3 = floor(length(which(temp<0.3))/ftSz) ### 42% [0,0.3] 
NUM0.7 = floor(length(which(temp>0.7))/ftSz) ### 17% [0.7,1], 41% [0.3,0.7],
NUM0.3_0.7 = ncol(MEs) - NUM0.3 -NUM0.7

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
mytestCounter = 300

RecCluster = vector(mode="integer",length=mytestCounter)

for(testCounter in 1:mytestCounter){
  ###construct different cluster numbers
  if(testCounter <= mytestCounter/3){
    MUSIM = rep(0,ftSz)
    SigmaSIM = diag(1,ftSz,ftSz)
  }else if(testCounter >mytestCounter/3 & testCounter <2*mytestCounter/3){
    MUSIM=rep(c(0,0.3),each=floor(ftSz/2))
    SigmaSIM=nearPD(SIMcov, conv.tol = 1e-7)$mat
  }else{MUSIM=rep(c(0,0.3,0.8),each=floor(ftSz/3))
  SigmaSIM=nearPD(SIMcov, conv.tol = 1e-7)$mat
  }
  
  X = mvrnorm(nrow(exprDat_norm),mu=MUSIM,Sigma = SigmaSIM)
  
  fks = vector(mode="numeric",length=5)
  for (myK in 1:5){fks[myK] <-FkStatistic(X,myK)[[1]]}
  
  if(length(which(fks<0.85))==0){RecCluster [testCounter] = 1
  }else{
    RecCluster[testCounter] = which(fks<0.85)[length(which(fks<0.85))]
  }
}

length(which(RecCluster==1))
length(which(RecCluster==2))
length(which(RecCluster==3))
length(which(RecCluster==4))

