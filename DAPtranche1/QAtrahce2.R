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

exprDat_normX = exprDat_norm[which(MetaRaw[,"diseaseGroup"]=="Injury"),!grepl("HybControlElution|NonBiotin|None",colnames(exprDat_norm))]

### Analysis 1.1 How many clusters are there?
pcStr <- prcomp(log10(as.matrix(exprDat_normX)),scale = TRUE)
pcDat = pcStr$x[,1:10]

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
fviz_nbclust(pcDat, kmeans, nstart = 25,  method = "gap_stat", nboot = 50)+
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
sigPrP=vector(mode="numeric")
sigPrOR=vector(mode="numeric")
sigPrLFC=vector(mode="numeric")
k=1
for (i in 1:ncol(exprDat_normX)){
  fit_glm <- glm(EndoLabel ~ exprDat_normX[,i], family = binomial(link = "logit"),control=list(maxit=50))
  temp <- summary(fit_glm)[["coefficients"]]
  Prz <- temp[nrow(temp),ncol(temp)]
  # not right now, check how to compare warning message
  # if (warnings()=="glm.fit: fitted probabilities numerically 0 or 1 occurred"){completeProId[i]=i}
  if (Prz < 0.05/ncol(exprDat_normX)) {sigPrInd[k]=i
  sigPrP[k]=Prz
  sigPrOR[k]=exp(temp[nrow(temp),1])
  sigPrLFC[k]=log(mean(exprDat_normX[which(EndoLabel==1),i])/mean(exprDat_normX[which(EndoLabel==0),i]))
  k=k+1}
}

### protein signature for endotype 1, Corresponding p value and odd ratio
sigPr = cbind(colnames(exprDat_normX)[sigPrInd],sigPrP,sigPrOR,sigPrLFC,EndoLabel[sigPrInd])
colnames(sigPr) <- c("Protein Signature","p values of z test", "odds ratio","log fold change","Endotype")
sigExp = exprDat_normX[,sigPrInd]

sigProGeneId = vector(mode="integer",length=nrow(sigPr))
for (sigProCount in 1:nrow(sigPr)){
  sigProGeneId[sigProCount] = which(ColTable[,"Protein Name"]==sigPr[,"Protein Signature"][sigProCount])
}

###Core talbe here: protein signature table (sigPrTable) & sig protein expression level table(sigExpTable)
sigPrTable = data.frame(cbind(sigPr,ColTable[sigProGeneId,]))
sigPrTable <- sigPrTable[which(sigPrTable$EntrezGeneID!=""),] 
sigExpTableP <- sigExp[,which(sigPrTable$EntrezGeneID!="")]### colnames as RawM
sigExpTable <- sigExpTableP ### creat a copy, with "database friendly col name"
colnames(sigExpTable) <- sigPrTable[,"EntrezGeneSymbol"]

### module analysis
### simplest 1-step network construction and module detection function
net = blockwiseModules(sigExpTable, power = 6,
                       TOMType = "unsigned", minModuleSize = 30, maxBlockSize=6000,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE, verbose = 3)

consMEs = net$MEs; ### eigen genes
moduleLabels = net$colors;
moduleColors = labels2colors(moduleLabels)
consTree = net$dendrograms[[1]];
plotDendroAndColors(consTree, moduleColors,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Consensus protein dendrogram and module colors based on euclidean distance")
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
minModuleSize = 30;
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
plotTOM = dissTOM^7;
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# Call the plot function
# TOMplot(plotTOM, tomTree, dynamicColors, main = "Network heatmap plot, all genes")

### Merging of modules whose expression profiles are very similar
# Calculate eigengenes
MEList = moduleEigengenes(sigExpTable, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
plot(METree, main = "Clustering of module proteins",
     xlab = "", sub = "")

MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(sigExpTable, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors
mergedMEs = merge$newMEs
plotDendroAndColors(tomTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


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
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes (300) injury group")
### in my opinion, TOMplot for selected genes says nothing!

### add sigPrTable col of correlation model label (here as color names) 
sigPrTable <- cbind(sigPrTable,mergedColors)

### Analysis 1.3: Bioinformatic characterisation of clusters
###pathway enrichment
rankPro <- as.numeric(sigPrTable[,"p.values.of.z.test"])
names(rankPro) <- sigPrTable[,"EntrezGeneID"]
ranks <- sort(rankPro,decreasing=TRUE)

### pathway analysis using MSigDB database
all_gene_sets <- gmtPathways("/Users/ydeng/Documents/QCstepOA/CorexpressionNetwork/msigdb.v7.4.entrez.gmt")
fgsea <- fgsea(pathways = all_gene_sets, stats=ranks, scoreType = "pos", eps = 0.0,minSize=15, maxSize=500)
head(fgsea)

topPathwaysUp <- fgsea[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgsea[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(all_gene_sets[topPathways], ranks, fgsea,gseaParam=0.3)

### display only independent pathways 
collapsedPathways <- collapsePathways(fgsea[order(pval)][1:10],all_gene_sets, ranks)
mainPathways <- fgsea[pathway %in% collapsedPathways$mainPathways][order(-NES), pathway]
plotGseaTable(all_gene_sets[mainPathways], ranks, fgsea, gseaParam = 0.3)

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
for (geneNameId in 1:nrow(sigPrTable)){
  if(!any(SubLocation2[,"Gene.name"]== sigPrTable[,"EntrezGeneSymbol"][geneNameId])){LocId[geneNameId]=""}
  else{LocId[geneNameId] = which(SubLocation2[,"Gene.name"]== sigPrTable[,"EntrezGeneSymbol"][geneNameId])}
}
SubLocationOrder = SubLocation2[which(LocId!=""),"Main.location"]
SubLocationDat = cbind(sigPrTable[which(LocId!=""),],as.matrix(SubLocationOrder))
colnames(SubLocationDat)[7] = "Subcellular location"
SubLocationType = levels(as.factor(SubLocationDat[,"Subcellular location"]))
### consider how to classify multiple location cases

###construct enquiry geneset based on subcellular location
CytosolCase = sigPrTable[grep("Cytosol",SubLocationOrder),"EntrezGeneID"]
NucleoliCase = sigPrTable[grep("Nucleoli",SubLocationOrder),"EntrezGeneID"]
Exosomoe = read.table("/Users/ydeng/Documents/QCstepOA/CorexpressionNetwork/exosome.tsv",header=TRUE,sep="\t")$GeneSymbol
ExosomoeCase = sigPrTable[sigPrTable[,"EntrezGeneSymbol"]%in%Exosomoe,"EntrezGeneID"]
GeneSetSublocation = list("Cytosol"=CytosolCase,"Nucleoli"=NucleoliCase,"Exosome"=ExosomoeCase) 

### subcellular location enrichment test
fgseaLoca <- fgsea(pathways = GeneSetSublocation, stats=ranks, scoreType = "pos",minSize=15,maxSize=500)
###FDR (Benjamini–Hochberg procedure for adjusted padj)
plotGseaTable(GeneSetSublocation, ranks, fgseaLoca, gseaParam = 0.5)### need to change the plot title "Pathway"
plotEnrichment(GeneSetSublocation[["Cytosol"]],ranks) + labs(title="Cytosol")  
plotEnrichment(GeneSetSublocation[["Nucleoli"]],ranks) + labs(title="Nucleoli")
plotEnrichment(GeneSetSublocation[["Exosome"]],ranks) + labs(title="Exosome")

### network analysis in cytoscape, including 
### (1)PPI network incorporate with expression level
### (2)correlation expression network and generic enrichment map

### (1) PPI network then incorporate with our expression level
netProDat <- sigPrTable[,c("EntrezGeneSymbol","EntrezGeneID","log.fold.change","p.values.of.z.test","mergedColors")]
string_interaction_cmd <- paste('string protein query taxonID=9606 limit=150 cutoff=0.9 query="',paste(netProDat$EntrezGeneSymbol[1:50], collapse=","),'"',sep="")
commandsGET(string_interaction_cmd)

###incorporate our correlation modules into network, modul color from WGCNA
loadTableData(netProDat[,c("EntrezGeneSymbol","log.fold.change","mergedColors")],table.key.column = "display name",data.key.column = "EntrezGeneSymbol")  #default data.frame key is row.names
### then in cytoscape "Style-> fill color -> column = merged colors + Mapping type = "Passthrough";" 

### (2) correlation expression network and generic enrichment map
sig_cor <- dissTOM[1:100,1:100] ### correlation matrix use dissTom from correlation module analysis
colnames(sig_cor) <- sigPrTable[,"EntrezGeneSymbol"][1:100]
rownames(sig_cor) <- sigPrTable[,"EntrezGeneSymbol"][1:100]
# sig_cor[row(sig_cor) == col(sig_cor)] <- 0 ### diagonal set from 1 to 0 elimnate self correlation
sig_cor[which(sig_cor<0.90)] <- 0 ### set hard threshold 0.9
sig_cor <- sig_cor[which(rowSums(sig_cor)!= 0),which(colSums(sig_cor) !=0)] ### remove all 0 protein from correlation matrix

#write out the correlation file
correlation_filename <- file.path(getwd(),"CorexpressionNetwork","cor_matrix.txt")
write.table(sig_cor,file = correlation_filename, col.names  = TRUE, row.names = FALSE, sep = "\t", quote=FALSE)

amat_url <- "aMatReader/v1/import"
amat_params = list(files = list(correlation_filename),delimiter = "TAB",undirected = FALSE,ignoreZeros = TRUE,interactionName = "correlated with",rowNames = FALSE)

### display network based on dissTOM
response <- cyrestPOST(operation = amat_url,body = amat_params,base.url ="http://localhost:1234")
current_network_id <- response$data["suid"]

current_nodetable_colnames <- getTableColumnNames(table="node",  network =  current_network_id)

layoutNetwork('cose',network = as.numeric(current_network_id))

### combine our information to nodetable
loadTableData(netProDat,table.key.column = "name",data.key.column = "EntrezGeneSymbol")  #default data.frame key is row.names

# enrichment analysis. function returns a data frame in the generic EM file format.
runGprofiler <- function(genes,current_organism = "hsapiens", 
                         significant_only = F, set_size_max = 200, 
                         set_size_min = 3, filter_gs_size_min = 5 , exclude_iea = F){
  
  gprofiler_results <- gprofiler(genes ,
                                 significant=significant_only,ordered_query = F,
                                 exclude_iea=exclude_iea,max_set_size = set_size_max,
                                 min_set_size = set_size_min,
                                 correction_method = "fdr",
                                 organism = current_organism,
                                 src_filter = c("GO:BP","REAC"))
  
  #filter results
  gprofiler_results <- gprofiler_results[which(gprofiler_results[,'term.size'] >= 3
                                               & gprofiler_results[,'overlap.size'] >= filter_gs_size_min ),]
  
  # gProfileR returns corrected p-values only.  Set p-value to corrected p-value
  if(dim(gprofiler_results)[1] > 0){
    em_results <- cbind(gprofiler_results[,
                                          c("term.id","term.name","p.value","p.value")], 1,
                        gprofiler_results[,"intersection"])
    colnames(em_results) <- c("Name","Description", "pvalue","qvalue","phenotype","genes")
    
    return(em_results)
  } else {
    return("no gprofiler results for supplied query")
  }
}

###Create an enrichment map with the returned g:Profiler results.
current_node_table <- getTableColumns(table= "node",network = as.numeric(current_network_id))
em_results <- runGprofiler(current_node_table$name)
em_results_filename <-file.path(getwd(),"CorexpressionNetwork",paste("gprofiler_cluster_enr_results.txt",sep="_"))
write.table(em_results,em_results_filename,col.name=TRUE,sep="\t",row.names=FALSE,quote=FALSE)

#write out the g:Profiler results
em_command = paste('enrichmentmap build analysisType="generic" ', 
                   'pvalue=',"0.05", 'qvalue=',"0.05",
                   'similaritycutoff=',"0.25",
                   'coeffecients=',"JACCARD",
                   'enrichmentsDataset1=',em_results_filename ,
                   sep=" ")

#ereturn the suid of newly created network.
em_network_suid <- commandsRun(em_command)
renameNetwork("Enrichmentmap", network=as.numeric(em_network_suid))

### upstream analysis
upstreamDat = sigPrTable[,c("EntrezGeneSymbol","log.fold.change","p.values.of.z.test")]
colnames(upstreamDat) <- c("entrez", "fc", "pvalue") ### name required by the package
upstreamDat$fc = as.numeric(upstreamDat$fc) 
upstreamDat$pvalue = as.numeric(upstreamDat$pvalue)
uniqueID=vector()
for (uniqueCouter in 1:length(unique(upstreamDat$entrez))){
  uniqueID[uniqueCouter] = which(upstreamDat$entrez== unique(upstreamDat$entrez)[uniqueCouter])
  
}
upstreamDatF = cbind(upstreamDat[uniqueID,]) ### package require unique gene ID
quaternary_results <- RunCRE_HSAStringDB(upstreamDatF, method = "Quaternary",
                                         fc.thresh = log(1.3), pval.thresh = 0.05,
                                         only.significant.pvalues = TRUE,
                                         significance.level = 0.05,
                                         epsilon = 1e-16, progressBar = FALSE,
                                         relations = NULL, entities = NULL)
quaternary_results[1:4, c("uid","symbol","regulation","pvalue")]

###Analysis 1.4 Clinical characteristics of endotypes
### read in clinical data
stepUpID <- as.matrix(read_excel("/Users/ydeng/Documents/QCstepOA/clinic/STEpUPOA_DATA_2021-04-21_0820.xlsx",range="A1:A141", col_names = TRUE))
age_sex <- read_excel("/Users/ydeng/Documents/QCstepOA/clinic/STEpUPOA_DATA_2021-04-21_0820.xlsx",range="J1:K141", col_names = TRUE)
kl_grade_worst <- read_excel("/Users/ydeng/Documents/QCstepOA/clinic/STEpUPOA_DATA_2021-04-21_0820.xlsx",range="IF1:IF141", col_names = TRUE)
womac_pain_score <- read_excel("/Users/ydeng/Documents/QCstepOA/clinic/STEpUPOA_DATA_2021-04-21_0820.xlsx",range="QQ1:QQ141", col_names = TRUE)
bmi <- read_excel("/Users/ydeng/Documents/QCstepOA/clinic/STEpUPOA_DATA_2021-04-21_0820.xlsx",range="O1:O141", col_names = TRUE)
cohort <- read_excel("/Users/ydeng/Documents/QCstepOA/clinic/STEpUPOA_DATA_2021-04-21_0820.xlsx",range="D1:D141", col_names = TRUE)
clinic <- cbind(stepUpID,cohort,age_sex,bmi,kl_grade_worst,womac_pain_score)
### generate cliniDat: table combined clinic information and proteomics 
matchStepID = vector(mode="integer",length=length(stepUpID))
for (StepIDcounter in 1:length(stepUpID)){
  matchStepID[StepIDcounter] = grep(stepUpID[StepIDcounter],as.matrix(metadata_reord[,"STEpUP Participant Identification Number (PIN)"]))
}
matchRowName = rownames(MetaRaw)[matchStepID]

matchSigPrID = vector(mode="integer",length=length(matchRowName))
for (matchSigCounter in 1:length(matchSigPrID)){
  matchSigPrID[matchSigCounter] = which(rownames(sigExpTable) == matchRowName[matchSigCounter])
}
clinicDatP <- cbind(clinic,EndoLabel[matchSigPrID])
clinicDat <- clinicDatP[!is.na(clinicDatP$womac_pain_score)&!is.na(clinicDatP$kl_grade_worst)&!is.na(clinicDatP$bmi_imported),]
colnames(clinicDat)[ncol(clinicDat)] = "EndoLabel"
clinicDat$sex <- sub("f",0,clinicDat$sex)
clinicDat$sex <- sub("m",1,clinicDat$sex)
clinicDat[,-1] <- apply(clinicDat[,-1],2,as.numeric)
### apply Generalised Linear Mixed Model
clinic_glm <- glmer(clinicDat$EndoLabel~ cbind(clinicDat$age,clinicDat$womac_pain_score,clinicDat$kl_grade_worst,clinicDat$bmi_imported)  + (1 | clinicDat$sex) ,family = binomial)
summary(clinic_glm)

###Residuals vs Fitted;Normal Q-Q;Scale-Location;Residuals vs Leverage
model.diag.metrics <- augment(model)
par(mfrow = c(2, 2))
plot(model)

###overlay clinical characteristics on PCA and UMAP
umapX = umap(exprDat_normX)

###Analysis 2.1 Biomarker of clinical features
biomarkerFit1 <- lm(pain ~ exprDat_normX + age + sex + cohort)
biomarkerFit2 <- lm(radiographic  ~ exprDat_normX + age + sex + cohort)
###Diagnostic Plots


