# This code runs meta-neighbor in human

# FIRST NEED TO SET THE WORKING DIRECTORY TO THE MAIN FOLDER
# This is the folder where "Start.RData" is located

print("-----Run MetaNeighbor: https://github.com/maggiecrow/MetaNeighbor")

########################################################################################
print("Load relevant libraries and functions")

outputFolder  = "output/"
scriptsFolder = "R/"
inputFolder   = "data/"
dataFolder    = "MTG/"

library(gplots)
library(RColorBrewer)
library(dplyr);
library(feather)
library(readxl)
library(WGCNA)

source(paste0(scriptsFolder,"Support_2016-11-03-runMetaNeighbor.R"))  # MetaNeighbor code

########################################################################################
print("Read in the data")

anno       <- read_feather(paste(dataFolder,"anno.feather",sep="")) 
Expr.dat   <- feather(paste(dataFolder,"data.feather",sep=""))   # FPKM
Expr.dat   <- Expr.dat[match(anno$sample_id,Expr.dat$sample_id),] # Make sure the expression matches the sample information
datIn      <- as.matrix(Expr.dat[,names(Expr.dat)!="sample_id"])
rownames(datIn) <- Expr.dat$sample_id
datIn      <- t(datIn)

kpSamp     <- is.element(anno$class_label,c("GABAergic","Glutamatergic","Non-neuronal"))
cl3        <- factor(anno$cluster_label,levels = sort(unique(anno$cluster_label)))
names(cl3) <- colnames(datIn)
clustersF  <- droplevels(cl3[kpSamp])
clusts     <- levels(clustersF)
clustTypes <- anno$class_label[kpSamp]

# Reorder clusters to match the tree
load(paste0(dataFolder,"dend.rda"))
ord = match(labels(dend),clusts)
clusts     <- levels(clustersF) <- clusts[ord]
clustCols  <- anno$cluster_color[match(clusts,anno$cluster_label)]
cl3        <- factor(anno$cluster_label,levels = clusts)
names(cl3) <- colnames(datIn)
clustersF  <- droplevels(cl3[kpSamp])
clusters   <- sort(unique(as.character(clustersF)))

# Subset the data to only include genes that are expressed
norm.dat   <- datIn[,kpSamp]
norm.dat   <- norm.dat[rowMeans(norm.dat)>1,]


########################################################################################
print("Read in gene sets from Paul et al 2017 paper to directly compare results with previous data set.")

hgncS3   = read_excel(paste0(inputFolder,"Paul2017_Table S3.xlsx"), skip=3)
customS4 = read_excel(paste0(inputFolder,"Paul2017_Table S4.xlsx"), skip=2) # This is also not used in the current file I sent

genesets = list()
for (i in 1:dim(hgncS3)[1]){
 gs  = as.character(hgncS3[i,1])
 gns = as.character(hgncS3[i,])
 gns = toupper(gns[!is.na(gns)])
 gns = sort(intersect(gns,rownames(norm.dat)))
 if(length(gns)>=5) genesets[[gs]] = gns
}
for (i in 1:dim(customS4)[1]){
 gs  = as.character(customS4[i,1])
 gns = as.character(customS4[i,])
 gns = toupper(gns[!is.na(gns)])
 gns = sort(intersect(gns,rownames(norm.dat)))
 if(length(gns)>=5) genesets[[gs]] = gns
}

kpGene = NULL
for (set in names(genesets))   kpGene = c(kpGene,genesets[[set]])
kpGene = sort(unique(kpGene))


########################################################################################
print("Subset into two groups of up to 20 random cells per group, 10 times (to get error bars on our values).")

subSamp  = 20
numIters = 10 
useMeL <- rand.labL <- list()
for (j in 1:numIters){
 seed = j
 rand.lab = rep(0,length(clustersF))
 names(rand.lab) = clustersF
 for (cli in clusters){
  set.seed(seed)
  seed   = seed+1*j
  kp     = which(clustersF==cli)
  ord    = sample(1:length(kp),min(length(kp),subSamp))
  ord1   = ord[1:round(length(ord)/2)]
  ord2   = ord[(1+round(length(ord)/2)):length(ord)]
  rand.lab[kp[ord1]] = 1
  rand.lab[kp[ord2]] = 2
 }
 useMeL[[j]] = rand.lab>0
 rand.labL[[j]]  = rand.lab
}

########################################################################################
print("Build the variables and run MetaNeighbor (this takes **SEVERAL HOURS** to run)")

clType  = c("GABAergic","Glutamatergic","Non-neuronal")
clTypes = NULL
for (cli in clusters) clTypes = c(clTypes,names(sort(-table(anno$class_label[anno$cluster_label==cli])))[1])
names(clTypes) = clusters

suppressWarnings(dir.create(paste0(outputFolder,"auroc")))
aurocL = list()
for (j in 1:numIters){
 auroc = NULL
 for (cli in clusters) {
  print(paste(j,"-",cli))
  sameClass    = clustTypes == clTypes[cli]
  useMe        = useMeL[[j]]&sameClass
  rand.lab     = rand.labL[[j]]
  dataTmp      = norm.dat[kpGene,useMe]
  randexp.lab  = as.character(rand.lab[useMe])
  celltype.lab = t(t((clustersF[useMe]==cli)-1+1))
  rownames(celltype.lab) = names(clustersF[useMe])
 
  fn = paste0(outputFolder,"auroc/MN_tmp")
  AUROC.scores = suppressWarnings(run_MetaNeighbor(data = dataTmp, experiment_labels = randexp.lab, 
      celltype_labels = celltype.lab, genesets = genesets, file_ext=fn))
  auroc = rbind(auroc,AUROC.scores)
 }
 rownames(auroc) = clusters
 aurocL[[j]] = auroc 
 save(aurocL,file=paste0(outputFolder,"aurocValues.rda"))
}

# Average out the permutation runs (may not make sense to do this here)
auroc <- maxAuroc <- minAuroc <- aurocL[[1]]
for (j in 2:length(aurocL)) {
  auroc = aurocL[[j]] + auroc
  minAuroc = pmin(minAuroc,aurocL[[j]])
  maxAuroc = pmax(maxAuroc,aurocL[[j]])
}
auroc = auroc / length(aurocL)

# Estimate of expected range over 10 iterations 
mean(maxAuroc-minAuroc)
# [1] 0.2261661


########################################################################################
print("Summarize results for different functional categories.")

getAurocSD <- function(aurocIn,kp){
  tmp=NULL
  for (i in 1:length(aurocIn))
	tmp = rbind(tmp,colMeans(aurocIn[[i]][kp,],na.rm=TRUE))
  out = apply(tmp,2,sd,na.rm=TRUE)
  return(out)
}

meanAuroc <- sdAuroc <- NULL
for (cli in clType){
  meanAuroc = cbind(meanAuroc,colMeans(auroc[clTypes==cli,],na.rm=TRUE))
  sdAuroc = cbind(sdAuroc,getAurocSD(aurocL,clTypes==cli))
}
colnames(meanAuroc) <- paste("Mean",clType)
colnames(sdAuroc) <- paste("SD",clType)
numGenes  = NULL
for (nm in rownames(meanAuroc)) numGenes = c(numGenes,length(genesets[[nm]]))
meanAuroc = cbind(meanAuroc,sdAuroc,numGenes)
meanAuroc = meanAuroc[order(-rowMeans(meanAuroc[,1:2])),]
write.csv(meanAuroc,paste0(outputFolder,"meanAuroc_human.csv"))

print("These values will be compared with mouse V1/ALM in the next code document.")
