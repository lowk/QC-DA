### identify new pathways
library(WGCNA)
library(GENIE3)
library(igraph)
library(RCy3)
library(MCL)
library(dcGOR) 

exprDat_normX = exprDat_norm[which(MetaRaw[,"diseaseGroup"]=="OA"),]

powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(exprDat_normX, powerVector = powers, verbose = 5)
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

softPower = 6;
# adjacency = adjacency(sigExpTable, power = softPower)
adjacency = adjacency(exprDat_normX, power = softPower)
# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM
# Call the hierarchical clustering function
tomTree = hclust(as.dist(dissTOM), method = "average")

# Plot the resulting clustering tree (dendrogram)
plot(tomTree, xlab="", sub="", main = "Protein expression level clustering on TOM-based dissimilarity",labels = FALSE, hang = 0.04)
minModuleSize = 50;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = tomTree, distM = dissTOM,
                            deepSplit = 4, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
plotDendroAndColors(tomTree, dynamicColors, "Coexpression modules",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Protein dendrogram and module colors (OA group)")

# Merging of modules whose expression profiles are very similar
# Calculate eigengenes
MEList = moduleEigengenes(exprDat_normX, colors = dynamicColors)
MEs = MEList$eigengenes

### EntrezGeneID of proteins: membership in modules. clusterList for further clusterProfiler
clusterListID=list()
clusterName = levels(as.factor(dynamicColors))
for (clusterCounter in 1:length(clusterName)){
  clusterListID[[clusterCounter]] = which(dynamicColors==clusterName[clusterCounter])
}
names(clusterListID)= strsplit(paste("X",seq(1:length(MEs)),collapse=",",sep=""),",")[[1]]### such naming for compareCluster

### k master regulators will be investigated within per module
k=5
masterRegName=matrix(NA,ncol=k,nrow=length(clusterListID))
for (clusterCounter in 1:length(clusterListID)){
  mem = exprDat_normX[,clusterListID[[clusterCounter]]]
  corM = abs(cor(mem,mem))
  igraphMatrixM <- graph_from_adjacency_matrix(corM, mode = "undirected", weighted = TRUE,diag = FALSE, add.colnames = NULL, add.rownames = NA)
  degNODE <- degree(igraphMatrixM, mode="all") 
  betweenNODE<- betweenness(igraphMatrixM)
  closeNODE <-closeness(igraphMatrixM)
  
  vertex.name=names(betweenNODE[order(betweenNODE,decreasing=TRUE)[1:k]])
  # plot(igraphMatrixM, vertex.label=vertex.name, vertex.label.font=0.5, vertex.label.color="black",
  #      vertex.label.cex=.7, vertex.size=betweenNODE/100, vertex.color=names(table(mergedColors)[1]),edge.color="gray85")
  
  masterRegName[clusterCounter,] = vertex.name
}

masterRexp = exprDat_normX[,as.vector(masterRegName)]
weightMat <- GENIE3(cov(masterRexp,masterRexp), treeMethod="ET", K=7, nTrees=10000)
netRaw <-getLinkList(weightMat,threshold=0.05)
adjM <- get.adjacency(graph.edgelist(as.matrix(netRaw[,1:2]), directed=TRUE))
vertex.display = ColTable[which(ColTable[,"Protein Name"] %in% colnames(adjM) ),"EntrezGeneSymbol"]
gra <- graph_from_adjacency_matrix(adjM, mode = "directed", weighted = TRUE,diag = FALSE, add.colnames = NULL, add.rownames = NA)
plot(gra, vertex.label=vertex.display,vertex.size=1, vertex.color="red",vertex.label.font=0.5, vertex.label.cex=.5,edge.arrow.size=0.3,arrow.width=1,vertex.label.dist=1.5,vertex.label.color="black",rescale = TRUE)

### nothing special in CytoScape. STRING has its limitation
masterRegUniProID = ColTable[which(ColTable[,1] %in% as.vector(masterRegName)),2]
string_interaction_cmd <- paste('string protein query taxonID=9606 limit=150 cutoff=0.9 query="',masterRegUniProID,'"',sep="")
response <- commandsGET(string_interaction_cmd)

### MCL on proposed pathway
MCLcluster <- mcl(x = adjM, addLoops = TRUE, allow1 = FALSE)

MCLmem = rownames(adjM)[which(MCLcluster$Cluster==3)]

masterRegUniProID2 = ColTable[which(ColTable[,1] %in% MCLmem),2]

plot(gra, vertex.label=vertex.display, vertex.color=MCLcluster$Cluste, vertex.size=3, vertex.label.font=0.5, vertex.label.cex=.5,edge.arrow.size=0.3,arrow.width=1,vertex.label.dist=1.5,vertex.label.color="black",rescale = TRUE)

dcDAGdomainSim(gra, domains = NULL, method.domain = c("BM.average","BM.max","BM.complete", "average", "max"), method.term = c("Resnik", "Lin", "Schlicker", "Jiang", "Pesquita"), force = TRUE, fast = TRUE,
               parallel = TRUE, multicores = NULL, verbose = TRUE)


g <- dcRDataLoader('onto.GOMF')
Anno <- dcRDataLoader('SCOP.sf2GOMF')
dag <- dcDAGannotate(g, annotations=Anno, path.mode="shortest_paths",verbose=TRUE)
alldomains <- unique(unlist(nInfo(dag)$annotations))
domains <- sample(alldomains,10)
dnetwork <- dcDAGdomainSim(g=dag, domains=domains,
                           method.domain="BM.average", method.term="Resnik", parallel=FALSE,
                           verbose=TRUE)
dnetwork
data <- data.frame(aSeeds=c(1,0,1,0,1), bSeeds=c(0,0,1,0,1))
rownames(data) <- id(dnetwork)[1:5]
coutput <- dcRWRpipeline(data=data, g=dnetwork, parallel=FALSE)
coutput
