library(NbClust)
library(fpc) # Compute DBSCAN using fpc package
library(corrplot)
library(caret)
library(WGCNA)
options(stringsAsFactors = FALSE)
library(fgsea)
library(msigdbr)
library(QuaternaryProd)
library(RCy3)
library(gProfileR)
library(RColorBrewer)
library(EnrichmentBrowser)
library(RCurl)
library(tidyverse)
library(broom)
library(umap)
library(lme4)
library(enrichplot)
library(clusterProfiler)
library(pathview)
library(DOSE)
library(igraph)
library(ggvenn)
library(MASS)
library(RColorBrewer)
library(gplots)


### extract expression data for OA/Injury group, when extract exprDat_norm from MySoma non human proteins have been excluded
exprDat_normX = exprDat_norm[which(MetaRaw[,"diseaseGroup"]=="Injury"),]
# exprDat_normX = exprDat_norm

### Analysis 1.1 How many clusters are there?
pcStr <- prcomp(log10(as.matrix(exprDat_normX)),scale = TRUE)
topPCn <- which(get_eigenvalue(pcStr)$cumulative.variance.percent>80)[1]
pcDat = pcStr$x[,1:topPCn]

# for (clusterK in 1:10){
#   kmeanStr <- kmeans(pcDat, clusterK, iter.max = 50, nstart = 1)
#   fviz_cluster(kmeanStr, exprDat_normX, ellipse.type = "norm")
# }

# Elbow method
fviz_nbclust(pcDat, kmeans, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Elbow method")

# Silhouette method
fviz_nbclust(pcDat, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")

# Gap statistic
# nboot = 50 to keep the function speedy. 
# recommended value: nboot= 500 for your analysis.
# Use verbose = FALSE to hide computing progression.
set.seed(123)
fviz_nbclust(pcDat, kmeans, nstart = 25,  method = "gap_stat", nboot = 100)+
  labs(subtitle = "Gap statistic method")

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

fks = vector(mode="numeric",length=10)
for (myK in 1:10)
{fks[myK] <-FkStatistic(pcDat,myK)[[1]]}
idK = which(fks<0.85)

plot(1:10,fks,type="p",main="define cluster number",xlab="cluster number",ylab="f(K) statistic")
lines(1:10,fks,type="l",col="blue")
abline(h=0.85,lty=2)     
grid()
points(idK,fks[idK],pch=23,col="red")

### DBSCAN
db <- fpc::dbscan(pcDat, eps = 0.15, MinPts = 5)
# Plot DBSCAN results
plot(db, pcDat, main = "DBSCAN", frame = FALSE)

fviz_cluster(db, pcDat, stand = FALSE, frame = FALSE, geom = "point")

### NbClust: 
EndoLabelNbClust <- NbClust(data = pcDat ,distance = "euclidean", min.nc = 2, max.nc = 15, 
                            method = "kmeans", index = "all", alphaBeale = 0.1)$Best.partition

kmeanStr <- kmeans(pcDat, 2, iter.max = 50, nstart = 1)
EndoLabel <- replace(kmeanStr$cluster,kmeanStr$cluster==2,0)


### Analysis 1.2 Identifying characteristic proteins of clusters
fit_glm <- glm(EndoLabel ~ exprDat_normX, family = binomial(link = "logit"),control=list(maxit = 50))
###summary(fit_glm) shows all predictors non significant. 
###Possible reason: strong multicolinearity; Case of complete separation 

# corExpr <- cor(exprDat_normX) ### correlation matrix of protein expression level
# length(which(corExpr>0.8))
# car::vif(fit_glm) ### multicolinearity

sigPrInd = vector(mode="integer")
sigPrPl=vector(mode="numeric")
sigPrPt=vector(mode="numeric")
sigPrOR=vector(mode="numeric")
sigPrLFC=vector(mode="numeric")

###apply logistic regression & t test
for (i in 1:ncol(exprDat_normX)){
  fit_glm <- glm(EndoLabel ~ exprDat_normX[,i], family = binomial(link = "logit"),control=list(maxit=50))
  temp <- summary(fit_glm)[["coefficients"]]
  Prz1 <- temp[nrow(temp),ncol(temp)]
  Prz2 <-t.test(exprDat_normX[,i][EndoLabel==1],exprDat_normX[,i][EndoLabel==0])$p.value
  # not right now, check how to compare warning message
  # if (warnings()=="glm.fit: fitted probabilities numerically 0 or 1 occurred"){completeProId[i]=i}
  sigPrPl[i]=Prz1
  sigPrPt[i]=Prz2
  sigPrOR[i]=exp(temp[nrow(temp),1])
  sigPrLFC[i]=log2(mean(exprDat_normX[which(EndoLabel==1),i])/mean(exprDat_normX[which(EndoLabel==0),i]))
}
sigPrP <- p.adjust(sigPrPl,method="bonferroni")
#sigPrPt <- p.adjust(sigPrPt,method="BH")
sigPrInd = which(sigPrP<0.05)

### all protein
# allPr = cbind(colnames(exprDat_normX),sigPrP,sigPrOR,sigPrLFC)
# colnames(allPr) <- c("Protein Signature","p values of z test", "odds ratio","log2 fold change")
# 
# allProGeneId = vector(mode="integer",length=nrow(sigPr))
# for (allProCount in 1:nrow(allPr)){
#   allProGeneId[allProCount] = which(ColTable[,"Protein Name"]==allPr[,"Protein Signature"][allProCount])
# }
# 
# allPrTable = data.frame(cbind(allPr,ColTable[allProGeneId,]))
# allPrTable <- allPrTable[which(allPrTable$EntrezGeneID!=""),] 
# allExpTableP <- exprDat_normX[,which(allPrTable$EntrezGeneID!="")]### colnames as RawM
# allExpTable <- allExpTableP ### creat a copy, with "database friendly col name", sigExpTableP could be for reference or check if needed
# colnames(allExpTable) <- allPrTable[,"EntrezGeneSymbol"]
# 
# uniAll = unique(allPrTable[,"EntrezGeneID"])
# uniAllID = vector(mode="integer",length=length(uniAll)) ### write such transformation into a function, too many places in the code require such operation
# for (uniAllCounter in 1:length(uniAll)){
#   uniAllID[uniAllCounter] = which(allPrTable[,"EntrezGeneID"]==uniAll[uniAllCounter])[1]
# }

### core table for enrichment analysis: in clude all the proteins.
# allPrTable = allPrTable[uniAllID,]
# allExpTable = allExpTable[,uniAllID]

### protein signature for endotype 1, Corresponding p value and odd ratio
sigPr = cbind(colnames(exprDat_normX)[sigPrInd],sigPrP[sigPrInd],sigPrOR[sigPrInd],sigPrLFC[sigPrInd])
colnames(sigPr) <- c("Protein Signature","p values of z test", "odds ratio","log2 fold change")
sigExp = exprDat_normX[,sigPrInd]

sigProGeneId = vector(mode="integer",length=nrow(sigPr))
for (sigProCount in 1:nrow(sigPr)){
  sigProGeneId[sigProCount] = which(ColTable[,"Protein Name"]==sigPr[,"Protein Signature"][sigProCount])
}

sigPrTable = data.frame(cbind(sigPr,ColTable[sigProGeneId,]))
sigPrTable <- sigPrTable[which(sigPrTable$EntrezGeneID!=""),] 
sigExpTableP <- sigExp[,which(sigPrTable$EntrezGeneID!="")]### colnames as RawM
sigExpTable <- sigExpTableP ### creat a copy, with "database friendly col name", sigExpTableP could be for reference or check if needed
colnames(sigExpTable) <- sigPrTable[,"EntrezGeneSymbol"]

### enrichment analysis "duplicate gene names, fgsea may produce unexpected results", we select unique "EntrezGeneSymbol"
uniPro = unique(sigPrTable[,"EntrezGeneID"])
uniID = vector(mode="integer",length=length(uniPro)) ### write such transformation into a function, too many places in the code require such operation
for (uniCounter in 1:length(uniPro)){
  uniID[uniCounter] = which(sigPrTable[,"EntrezGeneID"]==uniPro[uniCounter])[1]
}

###Core talbe here: protein signature table (sigPrTable) & sig protein expression level table(sigExpTable)
sigPrTable = sigPrTable[uniID,]
sigExpTable = sigExpTable[,uniID]

### module analysis
### simplest 1-step network construction and module detection function
# net = blockwiseModules(sigExpTable, power = 6,
#                        TOMType = "unsigned", minModuleSize = 30, maxBlockSize=6000,
#                        reassignThreshold = 0, mergeCutHeight = 0.25,
#                        numericLabels = TRUE, pamRespectsDendro = FALSE,
#                        saveTOMs = TRUE, verbose = 3)
# 
# consMEs = net$MEs; ### eigen genes
# moduleLabels = net$colors;
# consTree = net$dendrograms[[1]];
# plotDendroAndColors(consTree, moduleColors,
#                     "Module colors",
#                     dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05,
#                     main = "Consensus protein dendrogram and module colors based on euclidean distance")
### 1-step should be enough for our analysis

### tried regression on eigen genes: some egien genes are perfect enough for the endotype clustering
endo_glm <- glm(EndoLabel ~ as.matrix(consMEs), family = binomial(link = "logit"),control=list(maxit=50))
### no significance can be detected when all eigns are involved. and "glm.fit: fitted probabilities numerically 0 or 1 occurred" reported

sigEigenPrInd = vector(mode="integer",length=11)
for (i in 1:11){
  endo_glm <- glm(EndoLabel ~ as.matrix(consMEs[,i]), family = binomial(link = "logit"),control=list(maxit=50))
  temp <- summary(endo_glm)[["coefficients"]]
  Prz <- temp[nrow(temp),ncol(temp)]
  # if (is.null(warnings())=="FALSE"){completeProId[i]=i}
  if (Prz < 0.05/ncol(consMEs)) {sigEigenPrInd[i]=i}
}
###number 8,9, egien have complete separation
sigEigenPr = consMEs[,c(3,5,6,7,10,11)]
endo_glm2 <- glm(EndoLabel ~ as.matrix(sigEigenPr), family = binomial(link = "logit"),control=list(maxit=200))

corrplot(cor(consMEs),order="hclust")
testCor = cor(consMEs)
which(testCor>0.7&testCor<0.99)
testCor[testCor>0.7&testCor<0.99]

### more sophisticated way to construct the coexpression network and module analysis
# Choose a set of soft-thresholding powers
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
adjacency = adjacency(sigExpTable, power = softPower);
# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM
# Call the hierarchical clustering function
tomTree = hclust(as.dist(dissTOM), method = "average");

# Plot the resulting clustering tree (dendrogram)
plot(tomTree, xlab="", sub="", main = "Protein expression level clustering on TOM-based dissimilarity",labels = FALSE, hang = 0.04)
minModuleSize = 50;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = tomTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
plotDendroAndColors(tomTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Protein dendrogram and module colors (injury group)")

# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
# plotTOM = dissTOM^7;
# # Set diagonal to NA for a nicer plot
# diag(plotTOM) = NA;
# Call the plot function
# TOMplot(plotTOM, tomTree, dynamicColors, main = "Network heatmap plot, all genes")

### Merging of modules whose expression profiles are very similar
# Calculate eigengenes
MEList = moduleEigengenes(sigExpTable, colors = dynamicColors)
MEs = MEList$eigengenes

### test correlation between coexpression model and endotype
SigPMeLogit = vector(mode="numeric",length=ncol(MEs))
SigPTtest = vector(mode="numeric",length=ncol(MEs))
for (i in 1:ncol(MEs)){
  fit_glm <- glm(EndoLabel ~ MEs[,i], family = binomial(link = "logit"),control=list(maxit=50))
  tempLogit <- summary(fit_glm)[["coefficients"]]
  PmeLogit = temp[nrow(temp),ncol(temp)]
  PmeTtest = t.test(MEs[,i][EndoLabel==1],MEs[,i][EndoLabel==0])$p.value
  SigPMeLogit[i]=signif(PmeLogit,digits=3)
  SigPTtest[i] = signif(PmeTtest,digits=3)
}
SigPMeLogitAdj = p.adjust(SigPMeLogit,method="bonferroni")
SigPTtestAdj = p.adjust(SigPTtest,method="bonferroni")
ModulName = names(table(dynamicColors))
MeEndo = rbind(SigPMeLogit,SigPMeLogitAdj,SigPTtest,SigPTtestAdj)
colnames(MeEndo) = ModulName
rownames(MeEndo) = c("p(logit)","padj(logit)","p(ttest)","padj(ttest)")
MeEndo

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
plot(METree, main = "Clustering of module proteins",
     xlab = "", sub = "")

MEDissThres = 0.4 #0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(sigExpTable, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors
table(mergedColors)
# Convert numeric lables into colors
mergedMEs = merge$newMEs
plotDendroAndColors(tomTree, mergedColors,
                    "coexpression module",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

### add sigPrTable col of correlation model label (here as color names) 
sigPrTable <- cbind(sigPrTable,mergedColors)

### EntrezGeneID of proteins: membership in modules. clusterList for further clusterProfiler
clusterList = list()
clusterListID=list()
clusterName = levels(as.factor(mergedColors))
for (clusterCounter in 1:length(clusterName)){
  clusterList[[clusterCounter]] = sigPrTable[which(mergedColors==clusterName[clusterCounter]),"EntrezGeneID"]
  clusterListID[[clusterCounter]] = which(mergedColors==clusterName[clusterCounter])
}
names(clusterList)=c("X1","X2","X3","X4","X5","X6","X7","X8") ### such naming for compareCluster
#names(clusterList)=c("X1","X2","X3") ### such naming for compareCluster
#names(clusterList)=c("X1","X2","X3","X4","X5")
## Tom plot is time consuming (as stated in package tutorial). For reproducibility, we set the random seed, and size =300
nSelect = 300
set.seed(10);
select = sample(ncol(sigExpTable), size = nSelect);
selectTOM = dissTOM[select, select];
# There’s no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = mergedColors[select];

# Taking the dissimilarity to a power, say 7, makes the plot more informative by effectively changing
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes (300) Injury group")
### in my opinion, TOMplot for selected genes says nothing!

### when all the protein involved, too computational expensive in igraph. We visualise in different modules
colnames(dissTOM) = sigPrTable[,"EntrezGeneSymbol"]
rownames(dissTOM) = sigPrTable[,"EntrezGeneSymbol"]

###MEDiss eigen protein distance.TOMnet selected proteins.
igraphMatrix <- graph_from_adjacency_matrix(MEDiss, mode = "undirected", weighted = TRUE,diag = FALSE, add.colnames = NULL, add.rownames = NA)

###node importance scores: degree, betweenness, closeness.
degNODE <- degree(igraphMatrix, mode="all") 
betweenNODE<- betweenness(igraphMatrix)
closeNODE <-closeness(igraphMatrix)

l <- layout.circle(igraphMatrix)
plot(igraphMatrix,layout=l,vertex.label=colnames(MEDiss),vertex.label.font=2, vertex.label.color=rgb(0.1,0.7,0.8,0.5),
     vertex.label.cex=.7, vertex.size=betweenNODE*6, edge.color="gray85")

### top n membership for each MEs
topN=20
corMembership = abs(cor(sigExpTable,mergedMEs))
MEgroup = colnames(corMembership)
memberME = list()
dissMember = matrix(0,nrow=nrow(sigExpTable),ncol=topN*length(MEgroup))
for (group in 0:(length(MEgroup)-1)){
  memberME[[1+group]] = sort(corMembership[,1+group],decreasing=TRUE)[1:topN]
  dissMember[,(1+group*topN):(topN+group*topN)] = sigExpTable[,names(memberME[[1+group]])]
}

memberAdj = abs(cor(dissMember,dissMember))  ###Tom matrix not working well here, so we directly use correlation matrix

igraphMatrixMEM <- graph_from_adjacency_matrix(memberAdj, mode = "undirected", weighted = TRUE,diag = FALSE, add.colnames = NULL, add.rownames = NA)
degNODE <- degree(igraphMatrixMEM, mode="all") 
betweenNODE<- betweenness(igraphMatrixMEM)
closeNODE <-closeness(igraphMatrixMEM)
vertex.name=vector()
for(vertexCounter in 1:length(MEgroup)){
  temp.name=names(memberME[[vertexCounter]])
  vertex.name=c(vertex.name,temp.name)
}
vertex.color = rep(sub("ME","",MEgroup),each=5)
plot(igraphMatrixMEM, vertex.label=vertex.name, vertex.label.font=0.5, vertex.label.color="black",
     vertex.label.cex=.7, vertex.size=betweenNODE/30, vertex.color=vertex.color,edge.color="gray85")

### Analysis 1.3: Bioinformatic characterisation of clusters

###pathway enrichment, all using R package.
rankPro <- as.numeric(sigPrTable[,"log2.fold.change"])
names(rankPro) <- sigPrTable[,"EntrezGeneID"]
rankPro2 <- as.numeric(sigPrTable[,"log2.fold.change"]) ### rownames using UniProt, convenient for KEGG identifier transgormation
names(rankPro2) <- sigPrTable[,"UniPro.ID"]
ranks <- sort(rankPro,decreasing=TRUE)
ranks2 <- sort(rankPro2,decreasing=TRUE)
rankName <- names(ranks)
rankName2 <-names(ranks2)

### enrichplot package
### barplot
enrich1 <- DOSE::enrichDGN(rankName,pvalueCutoff = 0.05,pAdjustMethod = "BH",universe,minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,readable = FALSE)
barplot(enrich1, showCategory=10)
enrich2 <- clusterProfiler::gseKEGG(geneList= ranks2,organism= "hsa", nPerm= 10000,minGSSize= 3,maxGSSize=2000,pvalueCutoff = 0.05,pAdjustMethod = "BH",keyType = "uniprot")
#enrich3 <- clusterProfiler::gseGO(geneList=ranks, ont ="ALL", keyType = "ENTREZID",exponent = 1,nPerm= 10000,minGSSize = 10,maxGSSize = 500,pvalueCutoff = 0.05, verbose = TRUE, OrgDb = 'org.Hs.eg.db',pAdjustMethod = "BH",by="fgsea")

###dotplot
p1 <- dotplot(enrich1, showCategory=20) + ggtitle("dotplot for DisGeNET")
p2 <- dotplot(enrich2, showCategory=20) + ggtitle("dotplot for KEGG")
plot_grid(p1, p2, ncol=2)

enrich11 <- setReadable(enrich1, 'org.Hs.eg.db', 'ENTREZID')
p3 <- cnetplot(enrich11, foldChange=ranks)
p33 <- cnetplot(enrich11, categorySize="pvalue", foldChange=ranks)
p333 <- cnetplot(enrich11, foldChange=ranks, circular = TRUE, colorEdge = TRUE)
cowplot::plot_grid(p3, p33, p333, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))

p4 <- cnetplot(enrich2, foldChange=ranks)
p44 <- cnetplot(enrich2, categorySize="pvalue", foldChange=ranks2)
p444 <- cnetplot(enrich2, foldChange=ranks2, circular = TRUE, colorEdge = TRUE)
cowplot::plot_grid(p4, p44, p444, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))

emapplot(enrich1, showCategory = 20)
emapplot(enrich2, showCategory = 20)
emapplot(enrich2,legend_n=2) 
emapplot(enrich2,pie="count", pie_scale=1.5, layout="kk")

ridgeplot(enrich2) ###ridge plot for GSEA result

### WGCNA modules based pathway 
xx <- compareCluster(clusterList, fun = "enrichKEGG", organism="hsa", pvalueCutoff=0.05)
emapplot(xx)
# emapplot(xx,legend_n=2) 
# emapplot(xx,pie="count")
# emapplot(xx,pie="count", pie_scale=1.5, layout="kk")
# cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])

### pubmed  
terms <- enrich2$Description[1:8]
pmcplot(terms, 2010:2020, proportion=FALSE)

### if have time, I hope to develop spreading power score for nodes based on hsa05022 pathway network
hsa04530<- pathview(gene.data=ranks2,pathway.id = "hsa04530",species="hsa",limit=list(gene=max(abs(ranks2)), cpd=1))

### (1) PPI network then incorporate with our expression level in CytoScape
### interested in modPro in co-express network. 
modPro = names(unlist(memberME)) ### top 5 in each module
modID = which(sigPrTable[,"EntrezGeneSymbol"] %in% modPro)
netProDat <- sigPrTable[modID,c("EntrezGeneSymbol","EntrezGeneID","log2.fold.change","p.values.of.z.test","mergedColors")]
string_interaction_cmd <- paste('string protein query taxonID=9606 limit=150 cutoff=0.9 query="',paste(netProDat$EntrezGeneSymbol, collapse=","),'"',sep="")
response <- commandsGET(string_interaction_cmd)
loadTableData(netProDat[,c("EntrezGeneSymbol","log2.fold.change","mergedColors")],table.key.column = "display name",data.key.column = "EntrezGeneSymbol")  #default data.frame key is row.names

### amplify betweenness, visualise again in Cytoscape
current_nodetable_colnames <- getTableColumnNames(table="node", network =  "OA1")
current_nodetable <- getTableColumns('node', current_nodetable_colnames,network = "OA1")
current_nodetable2 = as.data.frame(cbind(current_nodetable[,"display name"],current_nodetable[,"BetweennessCentrality"]*3000))
current_nodetable2[,2] <- as.numeric(current_nodetable2[,2])
colnames(current_nodetable2) = c("EntrezGeneSymbol2","BetweennessCentrality")
loadTableData(current_nodetable2,table.key.column = "display name",data.key.column = "EntrezGeneSymbol2")  #default data.frame key is row.names

new.nodes = selectNodes(nodes=as.list(modPro),by.col = "display name",preserve.current.selection = TRUE,network = "OA1")
selectNodes(nodes=new.nodes)
commandsPOST('diffusion diffuse') 
createSubnetwork('selected',subnetwork.name = 'OA2')

###incorporate our correlation modules into network, module color from WGCNA

### then in cytoscape "Style-> fill color -> column = merged colors + Mapping type = "Passthrough";" 

### (2) co-expression network and generic enrichment map in CytoScape
sig_cor <- memberAdj  ### correlation matrix use top membershipi proteins correlation module analysis
colnames(sig_cor) <- vertex.name
rownames(sig_cor) <- vertex.name
sig_cor[row(sig_cor) == col(sig_cor)] <- 0 ### diagonal set from 1 to 0 elimnate self correlation
sig_cor[which(sig_cor<0.8)] <- 0 ### set hard threshold 0.9
sig_cor <- sig_cor[which(rowSums(sig_cor)!= 0),which(colSums(sig_cor) !=0)] ### remove all 0 protein from correlation matrix

#write out the correlation file
correlation_filename <- file.path(getwd(),"CorexpressionNetwork","cor_matrix.txt")
write.table(sig_cor,file = correlation_filename, col.names  = TRUE, row.names = FALSE, sep = "\t", quote=FALSE)

amat_url <- "aMatReader/v1/import"
amat_params = list(files = list(correlation_filename),delimiter = "TAB",undirected = FALSE,ignoreZeros = TRUE,interactionName = "correlated with",rowNames = FALSE)

### display network based on memberAdj (top proteins regarding to membership)
response <- cyrestPOST(operation = amat_url,body = amat_params,base.url ="http://localhost:1234")

current_network_id <- response$data["suid"]

current_nodetable_colnames <- getTableColumnNames(table="node",  network =  current_network_id)

layoutNetwork('cose',network = as.numeric(current_network_id))

### combine our information to nodetable
loadTableData(netProDat,table.key.column = "name",data.key.column = "EntrezGeneSymbol")  #default data.frame key is row.names

# pathway enrichment analysis in CytoScape. function returns a data frame in the generic EM file format. Work well, still can be good alternative.
# runGprofiler <- function(genes,current_organism = "hsapiens", 
#                          significant_only = F, set_size_max = 200, 
#                          set_size_min = 3, filter_gs_size_min = 5 , exclude_iea = F){
#   
#   gprofiler_results <- gprofiler(genes ,
#                                  significant=significant_only,ordered_query = F,
#                                  exclude_iea=exclude_iea,max_set_size = set_size_max,
#                                  min_set_size = set_size_min,
#                                  correction_method = "fdr",
#                                  organism = current_organism,
#                                  src_filter = c("GO:BP","REAC"))
#   
#   #filter results
#   gprofiler_results <- gprofiler_results[which(gprofiler_results[,'term.size'] >= 3
#                                                & gprofiler_results[,'overlap.size'] >= filter_gs_size_min ),]
#   
#   # gProfileR returns corrected p-values only.  Set p-value to corrected p-value
#   if(dim(gprofiler_results)[1] > 0){
#     em_results <- cbind(gprofiler_results[,
#                                           c("term.id","term.name","p.value","p.value")], 1,
#                         gprofiler_results[,"intersection"])
#     colnames(em_results) <- c("Name","Description", "pvalue","qvalue","phenotype","genes")
#     
#     return(em_results)
#   } else {
#     return("no gprofiler results for supplied query")
#   }
# }
# 
# ###Create an enrichment map with the returned g:Profiler results.
# current_node_table <- getTableColumns(table= "node",network = as.numeric(current_network_id))
# em_results <- runGprofiler(current_node_table$name)
# em_results_filename <-file.path(getwd(),"CorexpressionNetwork",paste("gprofiler_cluster_enr_results.txt",sep="_"))
# write.table(em_results,em_results_filename,col.name=TRUE,sep="\t",row.names=FALSE,quote=FALSE)
# 
# #write out the g:Profiler results
# em_command = paste('enrichmentmap build analysisType="generic" ', 
#                    'pvalue=',"0.05", 'qvalue=',"0.05",
#                    'similaritycutoff=',"0.25",
#                    'coeffecients=',"JACCARD",
#                    'enrichmentsDataset1=',em_results_filename ,
#                    sep=" ")
# 
# #ereturn the suid of newly created network.
# em_network_suid <- commandsRun(em_command)
# renameNetwork("Enrichmentmap", network=as.numeric(em_network_suid))

###sublocation enrichment

### retrive sublocation from UniProt using UniProt.ws package
# availableUniprotSpecies(pattern="Homo sapiens") ### find taxID for human
# up <- UniProt.ws(taxId=9606) 
# 
# SubLocation=vector(mode="list",length=nrow(ColTable))###ColTable[,2] is UniProtId
# for (proCounter in 1:nrow(ColTable)) try({
#   UniProtId = ColTable[,2][proCounter]
#   SubLocation[[proCounter]]<- select(up,UniProtId,"SUBCELLULAR-LOCATIONS", "UNIPROTKB")
# })
### saveRDS(SubLocation, file = "SubLocation.rds")
### slow, 2-3 seconds for one UniProtId; but a good programmatically way to retrieve sublocation 

### downloaded csv from https://www.proteinatlas.org/about/download
SubLocation2 <- read.csv("subcellular_location.tsv",sep="\t")

### match the gene symbol between our sigPrTable and reference subcellular_location.tsv
LocId = vector(mode="integer",length=nrow(sigPrTable))
delSigID=vector(mode="integer")
delk=1
for (geneNameId in 1:nrow(sigPrTable)){
  if(!any(SubLocation2[,"Gene.name"]== sigPrTable[,"EntrezGeneSymbol"][geneNameId])){LocId[geneNameId]=""
  delSigID[delk]=geneNameId
  delk=delk+1}
  else{LocId[geneNameId] = which(SubLocation2[,"Gene.name"]== sigPrTable[,"EntrezGeneSymbol"][geneNameId])[1]
  }
}
SubLocationOrder = SubLocation2[which(LocId!=""),"Main.location"]
SubLocationDat = cbind(sigPrTable[-delSigID,],as.matrix(SubLocationOrder))
colnames(SubLocationDat)[ncol(SubLocationDat)] = "Subcellular location"
SubLocationTypeR = levels(as.factor(SubLocationDat[,"Subcellular location"]))

### extract all the sublocation types
k=1
SubLocationType=vector()
for (locationCouter in 1:length(SubLocationTypeR)){
  templocal = SubLocationTypeR[locationCouter]
  if(!grepl(";",templocal)){SubLocationType[k]=templocal
  k=k+1}
  else{templocal2 <- strsplit(templocal,";")[[1]]
  for(splitC in 1:length(templocal2)){SubLocationType[k+splitC-1]=templocal2[splitC]}
  k=k+length(templocal2)
  }
}

SubLocationTypeA = names(table(SubLocationType)) ### all the subtypes included in csv file
slices <- table(SubLocationType)
lbls <- names(table(SubLocationType))
pie(slices, labels = lbls, main="Pie Chart of Sub-cellular locations") ###composition of sublocation
###construct sublocation genesets: GeneSetSublocation
GeneSetSublocation=vector(mode="list",length=length(SubLocationTypeA))
for (subCounter in 1:length(SubLocationTypeA)){
  GeneSetSublocation[[subCounter]] = SubLocationDat[grep(SubLocationTypeA[subCounter],SubLocationOrder),"EntrezGeneID"]
} 

Exosomoe = read.table("/Users/ydeng/Documents/QCstepOA/CorexpressionNetwork/exosome.tsv",header=TRUE,sep="\t")$GeneSymbol
ExosomoeCase = SubLocationDat[SubLocationDat[,"EntrezGeneSymbol"]%in%Exosomoe,"EntrezGeneID"]
GeneSetSublocation[[length(SubLocationTypeA)+1]] = ExosomoeCase
names(GeneSetSublocation) = c(SubLocationTypeA,"Exosomoe")

rankPro3 <- as.numeric(SubLocationDat[,"p.values.of.z.test"])
names(rankPro3) <- SubLocationDat[,"EntrezGeneID"]
ranks3 <- sort(rankPro3,decreasing=TRUE)

### subcellular location enrichment test
fgseaLoca <- fgsea(pathways = GeneSetSublocation, stats=ranks3, scoreType = "pos",minSize=15,maxSize=500)
plotGseaTable(GeneSetSublocation, ranks3, fgseaLoca, gseaParam = 0.3)### need to change the plot title "Pathway"
plotEnrichment(GeneSetSublocation[["Exosomoe"]],ranks3) + labs(title="Exosome")

### pathway analysis using MSigDB database
all_gene_sets <- gmtPathways("/Users/ydeng/Documents/QCstepOA/CorexpressionNetwork/msigdb.v7.4.entrez.gmt")
fgsea <- fgsea(pathways = all_gene_sets, stats=ranks3, scoreType = "pos", eps = 0.0,minSize=15, maxSize=500)

topPathwaysUp <- fgsea[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgsea[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(all_gene_sets[topPathways], ranks3, fgsea,gseaParam=0.1)

### display only independent pathways 
collapsedPathways <- collapsePathways(fgsea[order(pval)][1:10],all_gene_sets, ranks3)
mainPathways <- fgsea[pathway %in% collapsedPathways$mainPathways][order(-NES), pathway]
plotGseaTable(all_gene_sets[mainPathways], ranks3, fgsea, gseaParam = 0.1)

### cell/tissue type enrichment 
cell.tissue <- read.csv("/Users/ydeng/Documents/QCstepOA/CorexpressionNetwork/normal_tissue.tsv",sep="\t")

### matching sigPrTable to cellType/tissueType
LocId = vector(mode="integer",length=nrow(sigPrTable))
delSigID=vector(mode="integer")
delk=1
for (geneNameId in 1:nrow(sigPrTable)){
  if(!any(cell.tissue[,"Gene.name"]== sigPrTable[,"EntrezGeneSymbol"][geneNameId])){LocId[geneNameId]=""
  delSigID[delk]=geneNameId
  delk=delk+1}
  else{LocId[geneNameId] = which(cell.tissue[,"Gene.name"]== sigPrTable[,"EntrezGeneSymbol"][geneNameId])[1]
  }
}
cellTypeOrder = cell.tissue[which(LocId!=""),"Cell.type"]
cellTypeDat = cbind(sigPrTable[-delSigID,],as.matrix(cellTypeOrder))
colnames(cellTypeDat)[ncol(cellTypeDat)] = "Cell Type"

tissueTypeOrder = cell.tissue[which(LocId!=""),"Tissue"]
tissueTypeDat = cbind(sigPrTable[-delSigID,],as.matrix(tissueTypeOrder))
colnames(tissueTypeDat)[ncol(tissueTypeDat)] = "Tissue Type"

###cellTypeDat and tissueTypeDat: dataframe with information for further enrichment test
cellTypeF <- cellTypeDat[,"Cell Type"]
slices <- table(cellTypeF)
lbls <- names(table(cellTypeF))
pie(slices, labels = lbls, main="Pie Chart of cell types",radius = 1, cex = 0.3) ###composition of cell types

### extract all the cell types in the tsv: cellType
cellType = levels(as.factor(cellTypeF))

###construct cellType genesets: CellSet
CellSet=vector(mode="list",length=length(cellType))
for (subCounter in 1:length(cellType)){
  CellSet[[subCounter]] = cellTypeDat[grep(cellType[subCounter],cellTypeOrder),"EntrezGeneID"]
} 
names(CellSet) = cellType 

### prepare corresponding p values
rankPro4 <- as.numeric(cellTypeDat[,"p.values.of.z.test"])
names(rankPro4) <- cellTypeDat[,"EntrezGeneID"]
ranks4 <- sort(rankPro4,decreasing=TRUE)

### Cell Type enrichment test
fgseaCell <- fgsea(pathways = CellSet, stats=ranks4, scoreType = "pos",minSize=15,maxSize=500)
topPathwaysUp <- fgseaCell[ES > 0][head(order(pval), n=30), pathway]
topPathwaysDown <- fgseaCell[ES < 0][head(order(pval), n=30), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(CellSet[topPathways], ranks4, fgseaCell, gseaParam=0.3)
#plotGseaTable(CellSet, ranks4, fgseaCell, gseaParam = 0.3)### need to change the plot title "Pathway"

### extract all the tissue types
tissueTypeF <- tissueTypeDat[,"Tissue Type"]
deleID = which(tissueTypeF=="N/A")
tissueTypeF <- tissueTypeDat[-deleID,"Tissue Type"]
tissueType = unique(tissueTypeF)
slices <- table(tissueType)
lbls <- names(table(tissueType))
pie(slices, labels = lbls, main="Pie Chart of tissue types") ###composition of cell types

###construct tissue type genesets: TissueSet
TissueSet=vector(mode="list",length=length(tissueType))
for (subCounter in 1:length(tissueType)){
  TissueSet[[subCounter]] = tissueTypeDat[grep(tissueType[subCounter],tissueTypeOrder),"EntrezGeneID"]
} 
names(TissueSet) = tissueType 

### prepare corresponding p values
rankPro5 <- as.numeric(tissueTypeDat[,"p.values.of.z.test"])
names(rankPro5) <- tissueTypeDat[,"EntrezGeneID"]
ranks5 <- sort(rankPro5,decreasing=TRUE)

### Tissue Type enrichment test
fgseaTissue <- fgsea(pathways = TissueSet, stats=ranks5, scoreType = "pos",minSize=15,maxSize=500)
plotGseaTable(TissueSet, ranks5, fgseaTissue, gseaParam = 0.2)### need to change the plot title "Pathway"
### repeated code, write into function when have time.


### upstream analysis
upstreamDat = sigPrTable[,c("EntrezGeneID","log2.fold.change","p.values.of.z.test")]
colnames(upstreamDat) <- c("entrez", "fc", "pvalue") ### name required by the package, package require unique gene ID
upstreamDat$fc = as.numeric(upstreamDat$fc)
upstreamDat$pvalue = as.numeric(upstreamDat$pvalue)

quaternary_resultsINJ <- RunCRE_HSAStringDB(upstreamDat, method = "Quaternary",
                                         fc.thresh = log2(1.3), pval.thresh = 0.05/nrow(upstreamDat),
                                         only.significant.pvalues = TRUE,
                                         significance.level = 0.05,
                                         epsilon = 1e-16, progressBar = FALSE,
                                         relations = NULL, entities = NULL)
goodPid = which(quaternary_resultsINJ[,"pvalue"]!=-1 & quaternary_resultsINJ[,"pvalue"]<0.05/nrow(upstreamDat) & quaternary_resultsINJ[,"symbol"]!="No-Symbol")
print(paste(length(goodPid),"significant upstream regulators are found."))

###Top 10 regulators
quaternary_resultsINJ[goodPid, c("uid","symbol","regulation","pvalue")][1:5,]
x= list("UP Regulation" = which(quaternary_resultsINJ[goodPid,"regulation"]=="up"),"Down Regulation" = which(quaternary_resultsINJ[goodPid,"regulation"]=="down"))
ggvenn(x,fill_color = c("#0073C2FF", "#CD534CFF"),stroke_size = 0.5, set_name_size = 4)

###Analysis 1.4 Clinical characteristics of endotypes
### read in clinical data
stepUpID <- as.matrix(read_excel("/Users/ydeng/Documents/QCstepOA/clinic/STEpUPOA_DATA_2021-04-21_0820.xlsx",range="A1:A141", col_names = TRUE))
age_sex <- read_excel("/Users/ydeng/Documents/QCstepOA/clinic/STEpUPOA_DATA_2021-04-21_0820.xlsx",range="J1:K141", col_names = TRUE)
kl_grade_worst <- read_excel("/Users/ydeng/Documents/QCstepOA/clinic/STEpUPOA_DATA_2021-04-21_0820.xlsx",range="IF1:IF141", col_names = TRUE)
womac_pain_score <- read_excel("/Users/ydeng/Documents/QCstepOA/clinic/STEpUPOA_DATA_2021-04-21_0820.xlsx",range="QQ1:QQ141", col_names = TRUE)
bmi <- read_excel("/Users/ydeng/Documents/QCstepOA/clinic/STEpUPOA_DATA_2021-04-21_0820.xlsx",range="O1:O141", col_names = TRUE)
cohort <- read_excel("/Users/ydeng/Documents/QCstepOA/clinic/STEpUPOA_DATA_2021-04-21_0820.xlsx",range="D1:D141", col_names = TRUE)
clinic <- cbind(stepUpID,cohort,age_sex,bmi,kl_grade_worst,womac_pain_score)

stepUpID <- as.matrix(read_excel("/Users/ydeng/Documents/QCstepOA/clinic/STEpUPOA_OMBMM_MENTOR_DATA_2021-05-20_1313.xlsx",range="A1:A259", col_names = TRUE))
age_sex <- read_excel("/Users/ydeng/Documents/QCstepOA/clinic/STEpUPOA_OMBMM_MENTOR_DATA_2021-05-20_1313.xlsx",range="J1:K259", col_names = TRUE)
kl_grade_worst <- read_excel("/Users/ydeng/Documents/QCstepOA/clinic/STEpUPOA_OMBMM_MENTOR_DATA_2021-05-20_1313.xlsx",range="IF1:IF259", col_names = TRUE)
womac_pain_score <- read_excel("/Users/ydeng/Documents/QCstepOA/clinic/STEpUPOA_OMBMM_MENTOR_DATA_2021-05-20_1313.xlsx",range="QQ1:QQ259", col_names = TRUE)
bmi <- read_excel("/Users/ydeng/Documents/QCstepOA/clinic/STEpUPOA_OMBMM_MENTOR_DATA_2021-05-20_1313.xlsx",range="O1:O259", col_names = TRUE)
cohort <- read_excel("/Users/ydeng/Documents/QCstepOA/clinic/STEpUPOA_OMBMM_MENTOR_DATA_2021-05-20_1313.xlsx",range="D1:D259", col_names = TRUE)
clinic <- cbind(stepUpID,cohort,age_sex,bmi,kl_grade_worst,womac_pain_score)
clinic <- cbind(stepUpID,cohort,age_sex,bmi)

### generate cliniDat: table combined clinic information and proteomics 
matchStepID = vector(mode="integer",length=length(stepUpID))
for (StepIDcounter in 1:length(stepUpID)){
  matchStepID[StepIDcounter] = grep(stepUpID[StepIDcounter],as.matrix(metadata_reord[,"STEpUP Participant Identification Number (PIN)"]))[1]
}
matchRowName = rownames(MetaRaw)[matchStepID]

matchSigPrID = vector(mode="integer",length=length(matchRowName))
for (matchSigCounter in 1:length(matchSigPrID)){
  matchSigPrID[matchSigCounter] = which(rownames(sigExpTable) == matchRowName[matchSigCounter])[1]
}
clinicDatP <- cbind(clinic,EndoLabel[matchSigPrID])
clinicDat <- clinicDatP[!is.na(clinicDatP$womac_pain_score)&!is.na(clinicDatP$kl_grade_worst)&!is.na(clinicDatP$bmi_imported),]
clinicDat <- clinicDatP[!is.na(clinicDatP$bmi_imported),]

colnames(clinicDat)[ncol(clinicDat)] = "EndoLabel"
clinicDat$sex <- sub("f",0,clinicDat$sex)
clinicDat$sex <- sub("m",1,clinicDat$sex)
clinicDat[,-1] <- apply(clinicDat[,-1],2,as.numeric)
clinicDat$KLOA <- clinicDat$kl_grade_worst
clinicDat$KLOA <- apply(as.matrix(clinicDat$KLOA,nrow=1),1,function(x){if(x>2){x=1}
  else{x=0}})  ### add column wether KLOA

### apply Generalised Linear Mixed Model
# clinic_glm <- glmer(clinicDat$EndoLabel~ cbind(clinicDat$age,clinicDat$womac_pain_score,clinicDat$kl_grade_worst,clinicDat$bmi_imported)  + (1 | clinicDat$sex) ,family = binomial)
# summary(clinic_glm)

sex.endo = table(clinicDat$sex,clinicDat$EndoLabel)
sex.endoTest <- chisq.test(sex.endo) 
sex.endoTest.p.value = sex.endoTest$p.value

KLOA.endo = table(clinicDat$KLOA,clinicDat$EndoLabel)
KLOA.endoTest <- chisq.test(KLOA.endo) 
KLOA.endoTest.p.value = KLOA.endoTest$p.value

age.endoTest <- glm(clinicDat$EndoLabel ~ clinicDat$age, family = binomial(link = "logit"))
age.endoTest.p.value <- summary(age.endoTest)[["coefficients"]][2,4]

bim.endoTest <- glm(clinicDat$EndoLabel ~ clinicDat$bmi_imported, family = binomial(link = "logit"))
bim.endoTest.p.value <-summary(bim.endoTest)[["coefficients"]][2,4]

pain.endoTest <- glm(clinicDat$EndoLabel ~ clinicDat$womac_pain_score, family = binomial(link = "logit"))
pain.endoTest.p.value <- summary(pain.endoTest)[["coefficients"]][2,4]

normKL = (clinicDat$kl_grade_worst - mean(clinicDat$kl_grade_worst))/sd(clinicDat$kl_grade_worst)
kl.endoTest <- glm(clinicDat$EndoLabel ~ normKL, family = binomial(link = "logit"))
kl.endoTest.p.value <- summary(kl.endoTest)[["coefficients"]][2,4]

pAll=c(sex.endoTest.p.value, KLOA.endoTest.p.value, age.endoTest.p.value, bim.endoTest.p.value, pain.endoTest.p.value, kl.endoTest.p.value)
#pAll=c(sex.endoTest.p.value,age.endoTest.p.value, bim.endoTest.p.value)
pAll = p.adjust(pAll,method="bonferroni")
endoTest = matrix(rep(pAll,each = 2),nrow=2)
colnames(endoTest) = c("sex","KLOA","age","BIM","pain","KL")
#colnames(endoTest) = c("sex","age","BIM")
rownames(endoTest) = c("OA","Injury")

my_palette = colorRampPalette(brewer.pal(8, "Blues"))(ncol(endoTest)*nrow(endoTest))
heatmap.2(endoTest, Colv = NA,col =my_palette,main="p values of clinical features for testing \n association with endotypes test \n among OA and Injury",
          margins =c(9,12),trace="none",cexRow=2,keysize=1.5,rowsep=c(1,1))
###RowSideColors=c("green","red")

exp = exprDat_normX[matchSigPrID,][!is.na(clinicDatP$womac_pain_score)&!is.na(clinicDatP$kl_grade_worst)&!is.na(clinicDatP$bmi_imported),sigPrInd]
pcTemp <- prcomp(log10(as.matrix(exp)),scale = TRUE)

### find the first topPCn explain the 80% variance
topPCn <- which(get_eigenvalue(pcTemp)$cumulative.variance.percent>80)[1]
pcDatTemp = pcTemp$x[,1:topPCn]
clinicEndo = cbind(pcDatTemp,clinicDat)

PCcount = seq(1:(topPCn+5))
ACVA = get_eigenvalue(pcTemp)$cumulative.variance.percent[1:(topPCn+5)]
ScreeDat = as.data.frame(cbind(PCcount,ACVA))
ggplot(ScreeDat,aes(x=PCcount)) + geom_point(y=ACVA,color="darkcyan") + geom_line(aes(y=ACVA, color="cyan"), linetype="dashed") +
  ggtitle("Variance Explained by PCs") + labs(x = "Principal Component", y ="Proportion of variance explained by PCs(%)", color = "") +
  theme(legend.position = "none",plot.title=element_text(size=18,hjust=0.5), axis.title=element_text(size = 13)) + 
  scale_x_continuous(breaks=seq(1,10,1))

ggpairs(clinicEndo, columns=1:topPCn, aes(color=as.vector(clinicEndo$clinicEndo)),
        diag=list(continuous=wrap("densityDiag",alpha=0.4)),
        lower=list(continuous = wrap("points",alpha=0.9,size=0.1)),
        upper = list(continuous = "blank"))

pairs(clinicEndo[,1:topPCn],
      col =c("green","red")[as.factor(clinicEndo$EndoLabel)],
      pch = c(5,8)[as.factor(clinicEndo$sex)], cex=0.8,                     
      main = "KL defined OA (shape) overlay with endotype(color) on PCA")

chisq.test(table(clinicEndo[,"KLOA"],clinicEndo[,"sex"]))$p.value

ggplot(data=clinicEndo,aes(x=as.factor(EndoLabel),y=normKL,color=as.factor(EndoLabel))) +  geom_violin() + geom_boxplot(width=0.1) + 
  labs(x="Endotype", y="Normalised KL score",color="Endotype") + ggtitle(paste("Normalised KL score distribution within endotypes p =", signif(kl.endoTest.p.value,digits=4))) +
  theme(plot.title=element_text(size=20,hjust=0.5),legend.title=element_text(size=20),legend.text=element_text(size = 20), axis.title=element_text(size = 20))

###overlay clinical characteristics on UMAP

###Analysis 2.1 Biomarker of clinical features
biomarkerFit1 <- lm(pain ~ exprDat_normX + age + sex + cohort)
biomarkerFit2 <- lm(radiographic  ~ exprDat_normX + age + sex + cohort)
###Diagnostic Plots
###Residuals vs Fitted;Normal Q-Q;Scale-Location;Residuals vs Leverage
model.diag.metrics <- augment(model)
par(mfrow = c(2, 2))
plot(model)


