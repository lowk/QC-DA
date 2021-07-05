### subsetting data, adjusted rand index to validat robustness of kmeans performed on our data
library(MASS)
library(NbClust)
library(mclust)
library(diceR)

#f(K) statistics:Nd=10
getALFA <- function(thisk,Nd){
  if (thisk == 2) {temp = 1 - 3/(4*Nd)}
  else{temp <- getALFA(thisk-1,Nd) + (1-getALFA(thisk-1, Nd))/6}
  return(temp)
}

FkStatistic <- function(pcDat,thisk){
  Nd = ncol(pcDat)
  kmeanStr <- kmeans(pcDat, thisk, iter.max = 50, nstart = 1) ###?iter.max should be???
  Sk = kmeanStr$tot.withinss
  if(thisk==1){fs=1}
  else if(FkStatistic(pcDat,thisk-1)[[2]]==0){fs=1}
  else{fs = Sk/(getALFA(thisk,Nd)*FkStatistic(pcDat,thisk-1)[[2]])}
  return(list(fs,Sk))
}

### extract expression data for OA/Injury group, when extract exprDat_norm from MySoma non human proteins have been excluded
exprDat_normX = exprDat_norm[which(MetaRaw[,"diseaseGroup"]=="OA"),]
# exprDat_normX = exprDat_norm

### Analysis 1.1 How many clusters are there?
pcStr <- prcomp(log10(as.matrix(exprDat_normX)),scale = TRUE)
topPCn <- which(get_eigenvalue(pcStr)$cumulative.variance.percent>80)[1]
pcDat = pcStr$x[,1:topPCn]

### calculate fks corresponding to cluster number 1:10
fks = vector(mode="numeric",length=10)
for (myK in 1:10)
{fks[myK] <-FkStatistic(pcDat,myK)[[1]]}

if(!any(fks<0.85)){print("only 1 cluster detected")
}else{
  ### propsed cluster number BestCN
  BestCN <- length(unique(NbClust(data = pcDat ,distance = "euclidean", min.nc = 2, max.nc = 15, 
                                  method = "kmeans", index = "all", alphaBeale = 0.1)$Best.partition))
  
  kmeanStr <- kmeans(pcDat, BestCN, iter.max = 100, nstart = 1)
  
  ClustType <- kmeanStr$cluster
  
  ###10-cross validation
  sampID <- 1:nrow(pcDat)
  similairyVK = vector("integer",length=10)
  intv = floor(nrow(pcDat)/10)
  
  for (validK in 1:10){
    testID <- seq(intv*(validK-1)+1,intv*validK)
    trainID <- which(!(sampID%in%testID))
    TrKCenter <- kmeans(pcDat[trainID,],BestCN, iter.max = 100, nstart = 1)$centers
    
    testClusterType= vector("integer",length=length(testID))
    for (testCounter in 1:length(testID)){
      
      distanceK = vector("integer",length=nrow(TrKCenter))
      ###assign cluster membership by finding the nearest center
      for(centerCounter in 1:nrow(TrKCenter)){
        distanceK[centerCounter] = dist(rbind(pcDat[testID[testCounter],],TrKCenter[centerCounter,]))
      }
      ### testClusetType: assigned by the nearest center from the test data
      testClusterType[testCounter] <- which(distanceK==min(distanceK))
    }
    ### similarity between "test data" and "orignial clustering based on the whole data".
    similairyVK[validK] =  mclust::adjustedRandIndex(ClustType[testID],testClusterType)
  }
  
  similarityAve = mean(similairyVK)
}


### PAC (Proportion of Ambiguous Clusters) based concensus clustering, using hierachical clustering and DBSCAN
CC1 <- diceR::consensus_cluster(pcDat, nk =2, p.item = 0.8, reps = 5,algorithms = c("hdbscan"))
consensus_evaluate(pcDat,CC1)
CC2 <- diceR::consensus_cluster(pcDat, nk =2, p.item = 0.8, reps = 5,algorithms = c("hc"))

sig_obj <- sigclust(pcDat, k = 2, nsim = 10, labflag = 0)

