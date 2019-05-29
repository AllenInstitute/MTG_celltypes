# This code runs meta-neighbor in mouse and compares the result to human

# FIRST NEED TO SET THE WORKING DIRECTORY TO THE MAIN FOLDER
# This is the folder where "Start.RData" is located

print("-----Run MetaNeighbor: https://github.com/maggiecrow/MetaNeighbor")

########################################################################################
print("Load relevant libraries and functions")

outputFolder  = "output/"
scriptsFolder = "R/"
inputFolder   = "data/"
vispFolder    = "VISp/"
almFolder     = "ALM/"

library(gplots)
library(RColorBrewer)
library(dplyr);
library(feather)
library(readxl)
library(WGCNA)

source(paste0(scriptsFolder,"Support_2016-11-03-runMetaNeighbor.R"))  # MetaNeighbor code

########################################################################################
print("Read in the data")

annoV <- read_feather(paste(vispFolder,"anno.feather",sep="")) 
ExprV <- feather(paste(vispFolder,"data.feather",sep=""))   
ExprV <- ExprV[match(annoV$sample_id,ExprV$sample_id),] # Make sure the expression matches the sample information

annoA <- read_feather(paste(almFolder,"anno.feather",sep="")) 
ExprA <- feather(paste(almFolder,"data.feather",sep=""))   
ExprA <- ExprA[match(annoA$sample_id,ExprA$sample_id),] # Make sure the expression matches the sample information

datV  <- as.matrix(ExprV[,colnames(ExprV)!="sample_id"])
datA  <- as.matrix(ExprA[,colnames(ExprA)!="sample_id"])
datM  <- log2(rbind(datV,datA)+1)
datM  <- t(datM)

annoM <- rbind(annoV,annoA)

# Subset to only include relevant clusters
clustersM <- scan(paste0(inputFolder,"mouseClusters.txt"),what="character",sep="\n")
kpSampM   <- is.element(annoM$cluster_label,clustersM)
datM      <- datM[,kpSampM]
annoM     <- annoM[kpSampM,]


#############################
print("Properly format the data")

kpSamp     <- is.element(annoM$class_label,c("GABAergic","Glutamatergic","Non-Neuronal","Endothelial"))
cl3        <- factor(annoM$cluster_label,levels = sort(unique(annoM$cluster_label)))
names(cl3) <- colnames(datM)
clustersF  <- droplevels(cl3[kpSamp])
clusts     <- levels(clustersF)
clustTypes <- annoM$class_label[kpSamp]

# Reorder clusters to match the tree
ord        <- match(clustersM,clusts)
clusts     <- levels(clustersF) <- clusts[ord]
clustCols  <- annoM$cluster_color[match(clusts,annoM$cluster_label)]
cl3        <- factor(annoM$cluster_label,levels = clusts)
names(cl3) <- colnames(datM)
clustersF  <- droplevels(cl3[kpSamp])
clusters   <- sort(unique(as.character(clustersF)))

# Subset the data to only include genes that are expressed
norm.dat   <- datM[,kpSamp]
norm.dat   <- norm.dat[rowMeans(norm.dat)>1,]


########################################################################################
print("Read in gene sets from Paul et al 2017 paper to directly compare results with previous data set.")

hgncS3   = read_excel(paste0(inputFolder,"Paul2017_Table S3.xlsx"), skip=3)
customS4 = read_excel(paste0(inputFolder,"Paul2017_Table S4.xlsx"), skip=2) # This is also not used in the current file I sent

genesets = list()
for (i in 1:dim(hgncS3)[1]){
 gs  = as.character(hgncS3[i,1])
 gns = as.character(hgncS3[i,])
 gns = sort(intersect(gns,rownames(norm.dat)))
 if(length(gns)>=5) genesets[[gs]] = gns
}
for (i in 1:dim(customS4)[1]){
 gs  = as.character(customS4[i,1])
 gns = as.character(customS4[i,])
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
print("Build the variables and run MetaNeighbor (this takes **APROXIMATELY THREE DAYS** to run)")

clType  = c("GABAergic","Glutamatergic","Non-Neuronal","Endothelial")
clTypes = NULL
for (cli in clusters) clTypes = c(clTypes,names(sort(-table(annoM$class_label[annoM$cluster_label==cli])))[1])
names(clTypes) = clusters
clType = c("GABAergic","Glutamatergic","Non-Neuronal")
clTypes[clTypes=="Endothelial"] = "Non-Neuronal"

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
 save(aurocL,file=paste0(outputFolder,"aurocValues.csv"))
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
mean(maxAuroc-minAuroc,na.rm=TRUE)
# [1] 0.2011957


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
write.csv(meanAuroc,paste0(outputFolder,"meanAuroc_mouse.csv"))







########################################################################################
########################################################################################
########################################################################################
########################################################################################


print("Compare mouse AUROC values for interneurons with human values calculated in previous code document.")

## Required libraries
library(pheatmap)
library(grid)
library(ggplot2)
library(ggrepel)

## Load AUROC
auroc.sp <- list()
aurocL <- read.csv(paste0(outputFolder,"meanAuroc_mouse.csv"))
colnames(aurocL)[1] <- "hgnc"
colnames(aurocL)[-1] <- paste0("mouse_", colnames(aurocL)[-1])
auroc.sp[["mouse"]] <- aurocL
aurocL <- read.csv(paste0(outputFolder,"meanAuroc_human.csv"))
colnames(aurocL)[1] <- "hgnc"
colnames(aurocL)[-1] <- paste0("human_", colnames(aurocL)[-1])
auroc.sp[["human"]] <- aurocL

## Merge AUROCs into one data frame
auroc.df <- merge(auroc.sp[["mouse"]], auroc.sp[["human"]], by = "hgnc")
# Fix mislabeled HGNC family
auroc.df$hgnc <- sub("X5..nucleotidases", 
                     "5-hydroxytryptamine receptors, G protein-coupled", 
                     auroc.df$hgnc)
write.csv(auroc.df, file = paste0(outputFolder,"human_mouse_auroc.csv"), row.names = FALSE)


## Plot the correlations between species and classes
pw.cor <- cor(auroc.df[, grep("Mean", colnames(auroc.df))])
ph1 <- pheatmap(pw.cor, clustering_method = "average", show_colnames = FALSE, 
         color = colorRampPalette(brewer.pal(n = 9, "RdPu"))(100))
pdf(paste0(outputFolder,"human_mouse_auroc_heatmap.pdf"), width = 5, height = 3)
grid.newpage()
grid.draw(ph1$gtable)
dev.off()

print(cor(auroc.df$mouse_Mean.GABAergic, auroc.df$human_Mean.GABAergic))
# [1] 0.9157582


## Plot the scatterplot for the paper
plot.lab <- which(abs(auroc.df$human_Mean.GABAergic - auroc.df$mouse_Mean.GABAergic) > 0.2 | 
                    (auroc.df$human_Mean.GABAergic > 0.85 & 
                    auroc.df$mouse_Mean.GABAergic > 0.85 &
                    auroc.df$human_numGenes < 50))
g.inh <- ggplot(auroc.df, aes(x = mouse_Mean.GABAergic, y = human_Mean.GABAergic)) +
  geom_abline(slope = 1, intercept = 0, color = "blue") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey") +
  geom_errorbarh(aes(xmin = mouse_Mean.GABAergic - mouse_SD.GABAergic, 
                     xmax = mouse_Mean.GABAergic + mouse_SD.GABAergic), 
                 color = "grey90") +
  geom_errorbar(aes(ymin = human_Mean.GABAergic - human_SD.GABAergic, 
                    ymax = human_Mean.GABAergic + human_SD.GABAergic),
                color = "grey90") +
  geom_point(aes(size = human_numGenes)) +
  scale_size_continuous(breaks = c(10, 40, 100)) +
  #geom_text_repel(data = auroc.df[plot.lab, ], aes(label = hgnc), size = 2) +  # Current version omits text
  xlab("Mouse classification accuracy (mean AUROC)") +
  ylab("Human classification accuracy (mean AUROC)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 
plot(g.inh)

ggsave(g.inh, filename = paste0(outputFolder,"Fig6a_human_mouse_inh_auroc.pdf"), width = 6.75, height = 5)
