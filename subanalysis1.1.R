### Subanalysis 1.1
library(factoextra)
library(NbClust)
library(M3C)
library(GGally)
library(fpc)
library(WGCNA)
library(sparcl)
library(diceR)
source("callMe1.R")

### read in csv data from S drive
myPath <- "/Volumes/Shared/STEpUP OA/STEPUP_DAG_releases/STEPUP_DAG_rel001_final/"
QCFrame <- read.csv(paste0(release_prefix,'clinical/discovery_QApheno_1.csv'))
ClinicFrame1 <- read.csv(paste0(release_prefix,'clinical/discovery_DAPpheno_1.csv'))
ClinicFrame2 <- read.csv(paste0(release_prefix,'clinical/discovery_DAPpheno_2.csv'))

### adjust the ClinicFrame1, ClinicFrame2, QCFrame,consistent with the row order of SomaFrame
### consider whether to include soma adata row names into the csv file
QCFrame <- adjustFrame(QCFrame,SomaFrame)
ClinicFrame1 <- adjustFrame(ClinicFrame1,SomaFrame)
ClinicFrame2 <- adjustFrame(ClinicFrame2,SomaFrame)
CombinedFrameTemp <- cbind(QCFrame,ClinicFrame1,ClinicFrame2,SomaFrame[which(SomaFrame$SampleId %in% QCFrame$sf_iknee_sample_id_number),])
CombinedFrame <- CombinedFrameTemp[,unique(colnames(CombinedFrameTemp))]
# startDatId = which(colnames(CombinedFrame)=="CRYBB2.10000.28") ### column where the RFU starts; col 1:startDatId-1 is the meta table


# library(pvclust) ### consider which package best display dendrogram

### Perform all the analysis in OA injury groups separately, unless otherwise specified.

# Calculate PCA
### Tom may not split the disease group like this, so just use the CombinedFrame
exprDat = CombinedFrame[,which(colnames(CombinedFrame)=="CRYBB2.10000.28"):ncol(CombinedFrame)]

pc_norm <- prcomp(exprDat,scale = TRUE)

# Select top PCs
### get eigen values
eig.val <- factoextra::get_eigenvalue(pc_norm) 
topPC = which(eig.val$cumulative.variance.percent>80)[1]
pcDat = pc_norm$x[,1:topPC]

### set which label would like to be viewd in the PCA plot
myLabel = CombinedFrame$sf_iknee_qc_group
PlotDat = data.frame(cbind(pcDat,myLabel))

numLabel = c(0,1,2,3)
wordLabel = c("OA","Disease","RA","Control")
for (ic in 1:length(numLabel)){
  myLabel = gsub(numLabel[ic],wordLabel[ic],myLabel)  
}
PlotDat$myLabel = myLabel

GGally::ggpairs(PlotDat, columns=1:3, aes(color= as.factor(myLabel)),
                diag=list(continuous=wrap("densityDiag",alpha=0.4)),
                lower=list(continuous = wrap("points",alpha=0.9,size=0.1)),
                upper = list(continuous = "blank"),
                legend = c(1,1)) + labs(fill = "Disease Group")

myUmap = umap(pcDat)
df <- data.frame(x = myUmap$layout[,1],
                 y = myUmap$layout[,2],
                 WhichBatch = factor(PlotDat$myLabel))
colnames(df) = c("D1","D2",colnames("Disease Group"))
ggplot(df, aes(x=D1, y=D2, color = df[,3])) + geom_point() + labs(title="Umap labelled by disease group",color="Disease Group") +
  theme(plot.title=element_text(size=14,hjust=0.5),legend.text=element_text(size = 10), axis.title=element_text(size = 12)) 



clinicType=0 ### select disease group

exprDatFrame = CombinedFrame[CombinedFrame$sf_iknee_qc_group==clinicType,]
exprDat = exprDatFrame[,which(colnames(CombinedFrame)=="CRYBB2.10000.28"):ncol(exprDatFrame)]
### exprDat_norm: normalised and combatted proteomics of a particular disease group.

# Calculate f(K) statistics, silhouette score, gap statistic and elbow methods, visualise them against cluster numbers.

### fks statistic using self written function FkStatisitc
fks = vector(mode="numeric",length=10)
for (myK in 1:10)
{fks[myK] <-FkStatistic(pcDat,myK)[[1]]}
idK = which(fks<0.85)

plot(1:10,fks,type="p",main="define cluster number",xlab="cluster number",ylab="f(K) statistic")
lines(1:10,fks,type="l",col="blue")
abline(h=0.85,lty=2)     
grid()
points(idK,fks[idK],pch=23,col="red")

### gap statistic 
set.seed(123)
factoextra::fviz_nbclust(pcDat, kmeans, nstart = 25,  method = "gap_stat", nboot = 100)+
  labs(subtitle = "Gap statistic method")

# Elbow method
factoextra::fviz_nbclust(pcDat, kmeans, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Elbow method")

# If exist f(K)<0.85 apply NbClust to define the optimal cluster number.
# Cluster top PCs using kmeans (with 10 random starts, i.e. nstart=10 argument)
if(length(idK)>0){
  EndoLabelNbClust <- NbClust::NbClust(data = pcDat ,distance = "euclidean", min.nc = 2, max.nc = 15, 
                                       method = "kmeans", index = "all", alphaBeale = 0.1)$Best.partition
  thisk = length(unique(EndoLabelNbClust))
  kmeanStr <- kmeans(pcDat, thisk, iter.max = 50, nstart = 10) 
}else{print("no significant underlying clusters are found")}

# PCA and UMAP to visualize clustering.
PlotDat = data.frame(cbind(pcDat,kmeanStr$cluster))
colnames(PlotDat)[topPC+1]="EndoLabel"
GGally::ggpairs(PlotDat, columns=1:3, aes(color= as.factor(EndoLabel)),
                diag=list(continuous=wrap("densityDiag",alpha=0.4)),
                lower=list(continuous = wrap("points",alpha=0.9,size=0.1)),
                upper = list(continuous = "blank"),
                legend = c(1,1)) + labs(fill = "EndoLabel")

M3C::umap(t(pcDat),labels=as.factor(kmeanStr$cluster),controlscale=TRUE,scale=3)

# Run different clustering algorithms (Hierarchical clustering, DBSCAN,Hard competitive learning, Neural Gas algorithm)
### hierarchical cluatering 
hier_dist <- dist(pcDat, method = 'euclidean')
hier_clust <- hclust(hier_dist, method = 'average')
plot(hier_clust)

### DBSCAN
db <- fpc::dbscan(pcDat, eps = 1, MinPts = 100)
unique(db$cluster)
# Plot DBSCAN results
plot(db, pcDat, main = "DBSCAN", frame = FALSE)

# Run with different parameter initialisations for each clustering algorithm.

# Generate n new sampled datasets by sampling with replacement from the original dataset, cluster each of these sampled datasets using kmeans with 10 random starts. 
sampleID <- sample(1:nrow(pcDat), nrow(pcDat), replace = TRUE, prob = NULL)
pcDatSamp <- pcDat[sampleID,]
kmeanSamp <- kmeans(pcDatSamp, thisk, iter.max = 50, nstart = 10) 
M3C::umap(t(pcDatSamp),labels=as.factor(kmeanSamp$cluster),controlscale=TRUE,scale=3)

# Clustering on feature space composed by eigenproteins deduced from WGCNA
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(exprDat_normStore, powerVector = powers, verbose = 5)
# Plot the results:
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power(better than hard thresholding)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

###select softpower based on connectivity and scale independence
softPower = 6;
adjacency = adjacency(exprDat_normStore, power = softPower)
# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM
# Call the hierarchical clustering function
tomTree = hclust(as.dist(dissTOM), method = "average")

# # Plot the resulting clustering tree (dendrogram)
# plot(tomTree, xlab="", sub="", main = "Protein expression level clustering on TOM-based dissimilarity",labels = FALSE, hang = 0.04)
minModuleSize = 50;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = tomTree, distM = dissTOM,
                            deepSplit = 4, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
MEList = moduleEigengenes(exprDat_normStore, colors = dynamicColors)
MEs = MEList$eigengenes ###MEs: eigen protein space

fks = vector(mode="numeric",length=10)
for (myK in 1:10)
{fks[myK] <-FkStatistic(MEs,myK)[[1]]}
idK = which(fks<0.85)

if(length(idK)>0){
  EndoLabelNbClust <- NbClust::NbClust(data = MEs ,distance = "euclidean", min.nc = 2, max.nc = 15, 
                                       method = "kmeans", index = "all", alphaBeale = 0.1)$Best.partition
  thisk = length(unique(EndoLabelNbClust))
  kmeanStr <- kmeans(MEs, thisk, iter.max = 50, nstart = 10) 
}else{print("no significant underlying clusters are found")}

M3C::umap(t(MEs),labels=as.factor(kmeanStr$cluster),controlscale=TRUE,scale=3)

# Sparse clustering on original protein measurement using sparscl
perm.out <- HierarchicalSparseCluster.permute(pcDat, wbounds=c(1.5,2:6),nperms=5) ###original expStore out of memory limit
sparsehc <- HierarchicalSparseCluster(dists=perm.out$dists,wbound=perm.out$bestw, method="complete")
par(mfrow=c(1,2))
plot(sparsehc) ### add color leaf labels later

km.perm <- KMeansSparseCluster.permute(pcDat,K=thisk,wbounds=seq(3,7,len=15),nperms=5)
km.out <- KMeansSparseCluster(pcDat,K=thisk,wbounds=km.perm$bestw)
pairs(pcDat[,1:3],col=km.out[[1]]$Cs)

# Perform k means clustering based on normalised tranche1, normalised tranche2, simply scaled both tranches.  

# Calculate adjusted Rand Index between original clustering and each robustness test data resources. 
tree2 <- cutree(hier_clust, k = thisk)
clusterSim <- mclust::adjustedRandIndex(kmeanStr$cluster,tree2)

### PAC (Proportion of Ambiguous Clusters) based consensus clustering, using 
cc1 <- dice(pcDat, nk = thisk, reps = 5, algorithms = c("km"),cons.funs = c("kmodes", "majority"))
BestCN2 <- cc1$indices$pac$k[which(cc1$indices$pac$KM==min(cc1$indices$pac$KM))]

cc2 <- dice(pcDat, nk = thisk, reps = 5, algorithms = c("hc"),cons.funs = c("kmodes", "majority"))
BestCN3 <- cc2$indices$pac$k[which(cc2$indices$pac$HC_Euclidean==min(cc2$indices$pac$HC_Euclidean))]
memHC <- as.numeric(apply(cc2$indices$trim$E.new[[1]],1,function(x){names(sort(summary(as.factor(x)), decreasing=T))[1]}))
similairyKmHc <- mclust::adjustedRandIndex(memHC,kmeanStr$cluster)

memDB <- fpc::dbscan(pcDat,eps = 10, MinPts = 50)

# Perform kmeans clustering (nstart=10) on combined datasets of disease groups.

# Compare clustering results among combined datasets with individual disease group..
