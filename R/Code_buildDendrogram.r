# This code builds the dendrogram 

# FIRST NEED TO SET THE WORKING DIRECTORY TO THE MAIN FOLDER
# This is the folder where "Start.RData" is located

#-----------------------------------------------------------------------------------------------
print("Read in the data and load relevant libraries")
outputFolder  = "output/"
scriptsFolder = "R/"
inputFolder   = "data/"
dataFolder    = "MTG/"

# Load these libraries
library(beeswarm)
library(WGCNA);
library(edgeR)   
library(feather)
library(dendextend)
library(monocle)
library(ggplot2)
library(dplyr)
options(stringsAsFactors=FALSE)

# Read in the extra scripts
source(paste0(scriptsFolder,"Support_extraFunctions.r"))
source(paste0(scriptsFolder,"Support_violin_functions.R"))
source(paste0(scriptsFolder,"Support_heatmap_functions.R"))

# Read in the data
Expr.dat <- feather(paste(dataFolder,"data.feather",sep=""))   # CPM
anno     <- read_feather(paste(dataFolder,"anno.feather",sep="")) 
exprData <- as.matrix(Expr.dat[,colnames(Expr.dat)[colnames(Expr.dat)!="sample_id"]])
dend     <- NULL


#---------------------------------------------------------------------------------------------
print("Determine top cluster marker genes here for building the tree and displaying in heatmaps and violin plots.")

# Save cluster info variable and excluded clusters
clusterType = anno$class_label 
cl          = anno$cluster_id
names(cl)   = Expr.dat$sample_id
includeClas = c("GABAergic","Glutamatergic","Non-neuronal")
excludeClas = sort(setdiff(clusterType,includeClas))
kpSamp      = !is.element(clusterType,excludeClas) 

# Get median expression per cluster and the proportions
normDat    = log2(exprData+1)
rownames(normDat) = Expr.dat$sample_id
normDat    = t(normDat)
exprThresh = 1
medianExpr = do.call("cbind", tapply(names(cl), cl, function(x) rowMedians(normDat[,x]))) 
propExpr   = do.call("cbind", tapply(names(cl), cl, function(x) rowMeans(normDat[,x]>exprThresh))) 
rownames(medianExpr) <- rownames(propExpr) <- rownames(normDat)  
clusterOrd = as.character(sort(unique(anno$cluster_id)))
clusterKp  = as.character(sort(unique(anno$cluster_id[kpSamp])))
medianExpr = medianExpr[,clusterOrd]
propExpr   = propExpr[,clusterOrd]

# Exclude non-informative genes from cluster names and make sure at least one informative gene is in the shiny
kpGn         = apply(medianExpr[,clusterKp],1,max)>0
shinyFail    = (grepl("\\-",rownames(medianExpr)))|is.element(substr(rownames(medianExpr),1,1),0:9)  # -'s crash the shiny, numbers add an "X"
excludeGenes = sort(rownames(medianExpr)[grepl("LOC",rownames(medianExpr))|grepl("LINC",rownames(medianExpr))|
  grepl("FAM",rownames(medianExpr))|grepl("ORF",rownames(medianExpr))|grepl("KIAA",rownames(medianExpr))|
  grepl("FLJ",rownames(medianExpr))|grepl("DKFZ",rownames(medianExpr))|
  grepl("RPS",rownames(medianExpr))|grepl("RPL",rownames(medianExpr))|shinyFail])  # Also exclude ribosomal genes

# Also exclude mitochondrial genes and genes on sex chromosomes  
sexGenes     = scan(paste0(inputFolder,"sex_genes.txt"),what="character",sep="\n")
mitoGenes    = scan(paste0(inputFolder,"mito_genes.txt"),what="character",sep="\n")
excludeGenes = sort(unique(c(excludeGenes,sexGenes,mitoGenes)))

keepGenes    = setdiff(rownames(medianExpr)[kpGn],excludeGenes)
topGenes     = getTopMarkersByProp(propExpr[keepGenes,clusterKp],1,NULL,medianExpr[keepGenes,clusterKp],fcThresh=0.5,minProp=0.5)
topGenesAll  = getTopMarkersByProp(propExpr[kpGn,clusterKp],1,NULL,medianExpr[kpGn,clusterKp],fcThresh=0.5,minProp=0.5)
print(paste("Fraction of LOCxxxxxx markers:",mean(substr(topGenesAll,1,3)=="LOC")))
# [1] "Fraction of LOCxxxxxx markers: 0.32"


#-----------------------------------------------------------------------------------------------
print("Build clustering trees for nuclei based on the cluster calls")

function(sampleInfo, # Subsetted Samp.dat file for samples from non-outlier clusters
  classNameColumn = "cluster_type_label", # Will also be added to kpColumns (this is inh vs. exc. vs. glia)
  layerNameColumn = "layer_label",    # Where is the layer info stored (NULL if none)
  matchNameColumn = "cellmap_label",  # Where is the comparison info stored (e.g., closest mapping mouse type)
  regionNameColumn = "Region_label",  # Where is the region info stored (NULL if none)
  newColorNameColumn = "cellmap_color",  # Where is the column for selecting new cluster colors (e.g., color of closest mapping mouse type) [NULL keeps current colors]
  otherColumns = NULL, # Additional columns to transfer from Samp.dat (Note: "cluster_id", "cluster_label", and "cluster_color" are required)
  topGenes = NULL,                      #  A vector of the top genes for each cluster, with the vector names as cluster_id values (NULL if none)
  classLevels = c("inh","exc","glia"), # A vector of the levels for broad classes.  Set to NULL for none.
  getLayer = function(x) return(as.numeric(substr(as.character(x),2,2))),  # Function for converting layer text to numeric layer.  MAY NEED TO UPDATE
  sep="_")  # For renaming



# Get cluster info
anno$layer_label = anno$brain_subregion_label
clusterInfoAll = updateAndOrderClusters(as.data.frame(anno),topGenes=topGenes, classLevels=c(includeClas,excludeClas),
  regionNameColumn=NULL,newColorNameColumn="cluster_color", classNameColumn="class_label",matchNameColumn="cluster_label")
sampleInfo  = as.data.frame(anno[kpSamp,])

sampleInfo$old_cluster_label = sampleInfo$cluster_label
sampleInfo$old_cluster_id    = sampleInfo$cluster_id
sampleInfo$old_cluster_color = sampleInfo$cluster_color

clusterInfo = clusterInfoAll[is.element(clusterInfoAll$class_label,includeClas),]

majorGenes  = c("GAD1","SLC17A7","TYROBP","AQP4","PDGFRA","OPALIN","NOSTRIN")  #CSPG4
majorLabels = c("Inh","Exc","Micro","Astro","OPC","Oligo","Endo")
broadGenes  = c("PAX6","LAMP5","VIP","SST","PVALB","LINC00507","RORB","THEMIS","FEZF2","TYROBP","FGFR3","PDGFRA","OPALIN","NOSTRIN")

keepGenes2  = c(keepGenes,majorGenes,broadGenes)
newCluster  = renameClusters(sampleInfo,clusterInfo, propExpr[keepGenes2,clusterKp], medianExpr[keepGenes2,clusterKp], 
    propDiff = 0, propMin=0.4, medianFC = 1, propLayer=0.1,excludeGenes = excludeGenes,
    majorGenes=majorGenes, majorLabels=majorLabels, broadGenes=broadGenes)
clusterInfo$cluster_label = newCluster[clusterInfo[,1]]

# Reorder cluster information
clustSplit = strsplit(clusterInfo$cluster_label," ")
clusterInfo$broadGene = as.character(lapply(clustSplit,function(x) return(x[3])))
clusterInfo$topGene = as.character(lapply(clustSplit,function(x) return(x[length(x)])))

# Read in and update new colors
#newCols = read.csv(paste0(inputFolder,"clusterColors.csv"),row.names=1)
#newCols = newCols[match(clusterInfo$old_cluster_label,rownames(newCols)),]
#clusterInfo$cluster_color = newCols$final_color
#clusterColor <- clusterInfo$cluster_color
#names(clusterColor) <- clusterInfo$cluster_label
#clusterInfo$lrank = newCols$mouse_order
clusterInfo = clusterInfo[order(as.numeric(as.character(clusterInfo$cluster_id))),]
clusterInfo$lrank = 1:dim(clusterInfo)[1]
l.rank      = setNames(as.integer(clusterInfo$lrank), clusterInfo$cluster_label)

#write.csv(cbind(clusterInfo[,1:2],newCluster[clusterInfo[,1]]),"newClusterNames_Final2.csv",row.names=FALSE)

# Subset data to only include expressed genes and desired non-outlier cell types
kpSamp2    <- as.character(clusterInfo$cluster_id)
kpGn2      <- apply(medianExpr[,kpSamp2],1,max)>0
medianExpr2 = medianExpr[kpGn2,kpSamp2]
propExpr2   = propExpr[kpGn2,kpSamp2]
colnames(medianExpr2) <- colnames(propExpr2) <- clusterInfo$cluster_label

# Determine a score for cell type specificity (marker/beta score)
specificityScoreRank <- getBetaScore(propExpr2,FALSE)

# Build and reorder the dendrogram
topNgenes  = 1200
dend <- getDend(medianExpr2[specificityScoreRank<=topNgenes,])
dend <- reorder.dend(dend,l.rank)
dend <- collapse_branch(dend, 0.01)
dend <- dend %>% set("leaves_pch", 19) %>% set("leaves_cex", 2) 
save(dend, file=paste0(outputFolder,"dend.RData"))     
save(dend, file=paste0(dataFolder,"dend.rda"))  

#-----------------------------------------------------------------------------------------------
print("Update the annotation file with the new cluster names, ids, and colors")

clusterInfo = clusterInfo[match(labels(dend),clusterInfo$cluster_label),]
# We read this in from the website, so writing out the same information is unnecessary
#clusterId <- clusterInfo$cluster_id <- 1:length(clusterInfo$cluster_label)
#names(clusterId) <- clusterLabel <- clusterInfo$cluster_label
#names(clusterLabel) <- clusterInfo$old_cluster_label

#sampleInfo$cluster_label = clusterLabel[sampleInfo$cluster_label]
#sampleInfo$cluster_id = clusterId[sampleInfo$cluster_label]
#sampleInfo$cluster_color = clusterColor[sampleInfo$cluster_label]

# Plot the dendrogram
pdf(paste0(outputFolder,"Fig1c_clusterDendrogramFinal.pdf"),height=12,width=20)
main = paste("Based on top",topNgenes,"specificity genes, correlation distance")
label_color2 = clusterInfo$cluster_color[match(dend %>% labels,clusterInfo$cluster_label)]
rankCompleteDendPlot(dend=dend,label_color=label_color2,main=main,node_size=4)
dev.off()

pdf(paste0(outputFolder,"Fig1c_clusterDendrogramFinal_small.pdf"),height=5,width=10)
main = paste("Based on top",topNgenes,"specificity genes, correlation distance")
label_color2 = clusterInfo$cluster_color[match(dend %>% labels,clusterInfo$cluster_label)]
rankCompleteDendPlot(dend=dend,label_color=label_color2,main=main,node_size=0.0001)
dev.off()

#-----------------------------------------------------------------------------------------------
# Not needed since this information is written in a previous code document
#print("Rewrite the updated data files to a new directory")

#newInputFolder = paste0(inputFolder,"humanNew")
#dir.create(newInputFolder)

#exprOut = Expr.dat[match(sampleInfo$sample_id,Expr.dat$sample_id),]
#desc <- read_feather(paste(inputFolder,"human/desc.feather",sep="")) 

#write_feather(exprOut,paste(newInputFolder,"/data.feather",sep=""))   # CPM
#write_feather(sampleInfo,paste(newInputFolder,"/anno.feather",sep="")) 
#write_feather(desc,paste(newInputFolder,"/desc.feather",sep="")) 
#save(dend, file=paste0(newInputFolder,"/dend.rda")) 

save(clusterInfo,file=paste0(dataFolder,"/clusterInfo.rda"))
## WE MIGHT NEED "desc.feather"



