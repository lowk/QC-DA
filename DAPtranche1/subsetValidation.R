library(MASS)
library(NbClust)
library(mclust)

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
  
  for (validK in 1:10){
    testID <- seq(10*(validK-1)+1,10*validK)
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
