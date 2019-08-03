# This code makes most of the violin plots and group-averaged heatmaps in the manuscript.

# FIRST NEED TO SET THE WORKING DIRECTORY TO THE MAIN FOLDER
# This is the folder where "Start.RData" is located


#-----------------------------------------------------------------------------------------------
print("Read in the data and load relevant libraries")
outputFolder  = "output/"
scriptsFolder = "R/"
dataFolder    = "MTG/"
mouseFolder   = "mouse/"
inputFolder   = "data/"
fdir          = dataFolder

# Load these libraries
library(beeswarm)
library(WGCNA);
library(edgeR)   
library(feather)
library(dendextend)
library(monocle)
library(ggplot2)
library(dplyr)
library(cowplot)
library(matrixStats)
options(stringsAsFactors=FALSE)

# Read in the extra scripts
source(paste0(scriptsFolder,"Support_extraFunctions.r"))
source(paste0(scriptsFolder,"Support_violin_functions.R"))
source(paste0(scriptsFolder,"Support_heatmap_functions.R"))

# Read in the data
Expr.dat <- feather(paste(dataFolder,"data.feather",sep=""))   # CPM
anno <- read_feather(paste(dataFolder,"anno.feather",sep="")) 
exprData <- as.matrix(Expr.dat[,colnames(Expr.dat)[colnames(Expr.dat)!="sample_id"]])
load(paste0(dataFolder,"/dend.rda")) 
load(paste0(dataFolder,"clusterInfo.rda")) 
all_clusters = as.numeric(clusterInfo$cluster_id)
  
#-----------------------------------------------------------------------------------------------
print("Make trimmed mean average heatmap for figure 1.")

broadGenes = c("GAD1","ADARB2","PAX6","LAMP5","VIP","LHX6","SST","PVALB","SLC17A7","LINC00507",
  "RORB","THEMIS","FEZF2","SLC1A3","PDGFRA","FGFR3","OPALIN","NOSTRIN","TYROBP")
  
broad_plot <- group_heatmap_plot(data_source = fdir,
                   genes = broadGenes,
                   group_by = "cluster",
                   clusters = all_clusters,
                   calculation = "trimmed_mean",
                   labelheight = 40,
                   showcounts = F)
ggsave(paste0(outputFolder,"Fig1_broad_markers.pdf"),broad_plot,height = 4, width = 10)


#-----------------------------------------------------------------------------------------------
print("Make trimmed mean average heatmap for figure 2D.")

RORB_Genes = c("SLC17A7","RORB","CNR1","PRSS12","ALCAM","MET","MME","NTNG1","HS3ST4","CUX2",
               "PCP4","GRIN3A","GRIK3","CRHR2","TPBG","POSTN","SMYD1")
  
RORB_plot <- group_heatmap_plot(data_source = fdir,
                   genes = RORB_Genes,
                   group_by = "cluster",
                   clusters = c(50:52,54:59,62),
                   calculation = "trimmed_mean",
                   labelheight = 40,
                   showcounts = F)
ggsave(paste0(outputFolder,"Fig7A_RORB_subtypes.pdf"),RORB_plot,height = 4, width = 6)


#-----------------------------------------------------------------------------------------------
## NOTE: Additional heatmaps are plotted after the violin plots
#-----------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------
print("Recalculate marker genes to see if any LOC genes would be better markers than the named genes.")

# Save cluster info variable and excluded clusters
clusterType = anno$class_label 
cl          = anno$cluster_id
names(cl)   = Expr.dat$sample_id
includeClas = c("inh","exc","glia")
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

kpGn       = rep(TRUE,dim(propExpr)[1])  #betaScore>=minBeta
specGenes  = getTopMarkersByPropNew(propExpr=propExpr[kpGn,], medianExpr=medianExpr[kpGn,], propDiff = 0, propMin=0.4, medianFC = 1)
topIsLoc   = (substr(specGenes,1,3)=="LOC")|(substr(specGenes,1,4)=="LINC")|(substr(specGenes,1,4)=="DKFZ")|     #DKFZp686K1684 = "PAX6-AS1"
   (substr(specGenes,nchar(specGenes)-3,nchar(specGenes))=="-AS1")

# Function for getting the plotting genes (to avoid typing it out a bunch of times)
getPlotGenes <-function(plotGenes,kp){
 for (k in which (kp)){
  plotGenes = c(plotGenes,ifelse(topIsLoc[k],specGenes[k],"none"))
  if(specGenes[k]!=gsub("-","\\.",specGenes[k])) plotGenes=c(plotGenes,gsub("-","\\.",specGenes[k])) # Convert - to . for plotting
  plotGenes = c(plotGenes,clusterInfo$topGene[k])
 }
 return(setdiff(plotGenes,"none"))
}  
   
#-----------------------------------------------------------------------------------------------
print("Make violin plots for figures 2-4.")

# Excitatory
kp = clusterInfo$class_label=="Glutamatergic"
plotGenes = getPlotGenes(c("LAMP5","LINC00507","RORB","THEMIS","FEZF2"),kp)
group_violin_plot(data_source = dataFolder, clusters=all_clusters[kp], genes=plotGenes)
ggsave(paste0(outputFolder,"Fig2_exc_violinPlots_specificGenes_withLocs.pdf"), width = 12, height = 15)

# Inhibitory 
kp = is.element(1:length(clusterInfo$class_label),1:27)   # ADARB2
plotGenes = getPlotGenes(c("ADARB2","PAX6","LAMP5","VIP","LHX6","SST","HTR3A"),kp)
group_violin_plot(data_source = dataFolder, clusters=all_clusters[kp], genes=plotGenes)
ggsave(paste0(outputFolder,"Fig3_inh_ADARB2_violinPlots_specificGenes_withLocs.pdf"), width = 14, height = 15)

kp = is.element(1:length(clusterInfo$class_label),28:45)  # LHX6
plotGenes = getPlotGenes(c("ADARB2","LHX6","SST","PVALB"),kp)
group_violin_plot(data_source = dataFolder, clusters=all_clusters[kp], genes=plotGenes)
ggsave(paste0(outputFolder,"Fig3_inh_LHX6_violinPlots_specificGenes_withLocs.pdf"), width = 11, height = 12)

# Glia      # Expand genes to plot manually...
kp    = clusterInfo$class_label=="Non-neuronal"
opc   = c("PDGFRA","COL20A1","OLIG2","PRRX1")  #CSPG4
astro = c("FGFR3","SLC14A1","DIO2","LINC00982","GFAP","AQP4","MT1F","ID3")  # WIF1
oligo = c("OPALIN","PLP1","MAG","KLK6") # "TMEM235","CD22","KLK6"
endo  = c("NOSTRIN","EBF1","ITIH5","EMCN") # "CLDN5","VIM","PALMD"
micro = c("TYROBP","C3","CX3CR1","CSF1R")  # LST1
plotGenes = c("SLC1A3",opc,astro,oligo,endo,micro)
group_violin_plot(data_source = dataFolder, clusters=all_clusters[kp], genes=plotGenes, max_width = 15)
ggsave(paste0(outputFolder,"Fig4_glia_violinPlots_specificGenes.pdf"), width = 4, height = 10)


#-----------------------------------------------------------------------------------------------
print("Select other violin plots for other figure panels.")

# Pvalb 
plotGenes = c("GAD1","PVALB","UNC5B","NOG","COL15A1")
group_violin_plot(data_source = dataFolder, clusters=all_clusters[39:45], genes=plotGenes)
ggsave(paste0(outputFolder,"FigED10d_PVALB_chandelier_markers.pdf"), width = 7, height = 7)

# ET types 
plotGenes = c("SLC17A7","NPTX1","FAM84B","POU3F1")
group_violin_plot(data_source = dataFolder, clusters=all_clusters[46:69], genes=plotGenes)
ggsave(paste0(outputFolder,"FigED12e_ET_cluster_markers.pdf"), width = 17, height = 6)

# Extended data figure on LOC genes being good markers
plotGenes = c("RORB","CNR1","LINC01164")
group_violin_plot(data_source = dataFolder, clusters=all_clusters[51:54], genes=plotGenes)
ggsave(paste0(outputFolder,"FigED8a_RORB_LOCs.pdf"), width = 6, height = 6)

plotGenes = c("RORB","LOC105376081")
group_violin_plot(data_source = dataFolder, clusters=all_clusters[51:54], genes=plotGenes)
ggsave(paste0(outputFolder,"FigED8c_RORB_LOCs.pdf"), width = 6, height = 6)

plotGenes = c("RORB","CRYM","LOC401134")
group_violin_plot(data_source = dataFolder, clusters=all_clusters[51:53], genes=plotGenes)
ggsave(paste0(outputFolder,"FigED8d_RORB_LOCs.pdf"), width = 6, height = 6)

plotGenes = c("PVALB","LOC102723415")
group_violin_plot(data_source = dataFolder, clusters=all_clusters[41:43], genes=plotGenes)
ggsave(paste0(outputFolder,"FigED8e_PVALB_LOCs.pdf"), width = 6, height = 6)


#-----------------------------------------------------------------------------------------------
## NOTE: JITTER PLOTS
#-----------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------
print("Make layer jitter plots for figures 2-4.")

# Excitatory
anno$layer_id     <- anno$brain_subregion_id
anno$layer_color  <- anno$brain_subregion_color
anno$layer_label  <- anno$brain_subregion_label
exc_layers <- build_layer_plot(anno, dend, 0, cluster_ids = which(clusterInfo$class_label=="Glutamatergic"),textSize=1.8, maxPerCluster=100)  # fillColor =
exc_plot <- plot_grid(exc_layers,align = "v",nrow = 1,rel_widths = 1,rel_heights = 1,labels = "")
ggsave(paste0(outputFolder,"Fig2_exc_jitterLayerPlot.pdf"), width = 8, height = 6)

# Inhibitory 
adarb2_layers <- build_layer_plot(anno,dend,0,cluster_ids = 1:27,textSize=1.8, maxPerCluster=100)
adarb2_plot <- plot_grid(adarb2_layers,align = "v",nrow = 1,rel_widths = 1,rel_heights = 1,labels = "")
ggsave(paste0(outputFolder,"Fig3_inh_ADARB2_jitterLayerPlot.pdf"), width = 8, height = 6)

pvalb_sst_layers <- build_layer_plot(anno,dend,0,cluster_ids = 28:45,textSize=1.8, maxPerCluster=100)
ps_plot <- plot_grid(pvalb_sst_layers,align = "v",nrow = 1,rel_widths = 1,rel_heights = 1,labels = "")
ggsave(paste0(outputFolder,"Fig3_inh_LHX6_jitterLayerPlot.pdf"), width = 8, height = 6)

# Glia
glia_layers <- build_layer_plot(anno, dend, 0, cluster_ids = which(clusterInfo$class_label=="Non-neuronal"),textSize=1.8, maxPerCluster=100)  # fillColor =
glia_plot <- plot_grid(glia_layers,align = "v",nrow = 1,rel_widths = 1,rel_heights = 1,labels = "")
ggsave(paste0(outputFolder,"Fig4_glia_jitterLayerPlot.pdf"), width = 4, height = 6)


#-----------------------------------------------------------------------------------------------
## NOTE: More heatmaps are here along with marker genes for Supplemental table 2
#-----------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------
print("Find the top genes in every cluster and show them all in a heatmap.")

# Read in the levelInfo
levelInfo  = read.csv(paste0(inputFolder,"FinalHumanMTGclusterAnnotation_update.csv"))
levels     = list(all = levelInfo$final.tree.order)
for (cn in c("level1","level2","level3")){
  cols = table(levelInfo[,cn])
  cols = names(cols)[cols>2]
  for (cnn in cols){
    levels[[cnn]] = levelInfo$final.tree.order[levelInfo[,cn]==cnn]
  }
}

#-----------------------------------------------------------------------------------------------
print("Find the top genes in each cluster (overall and for all three levels).")

N = 5
outGenesL <- list()

for (cn in names(levels)){
  kp <- levels[[cn]]
  outGenes <- NULL
  for (i in 1:N) {
    print(paste(cn,i))
    outTmp <- getTopMarkersByPropNew(propExpr = propExpr[,kp], medianExpr = medianExpr[,kp], 
                                          propDiff = 0, propMin=0.4, medianFC = 1, excludeGenes = outGenes)
    outGenes <- c(outGenes,outTmp)
  }
  dim(outGenes) = c(length(outGenes)/N,N)
  outTmp <- apply(outGenes,1,paste,collapse=", ")
  outTmp <- gsub(", none","",outTmp)
  names(outTmp) <- kp
  outGenesL[[cn]] <- outTmp
}

#-----------------------------------------------------------------------------------------------
print("Save these genes to a table and output")

outTable <- levelInfo[,c("final.cluster","final.tree.order","level1","level2","level3")]
outTable$specific_genes = outGenesL[["all"]]
outTable$level1_genes = NA
outTable$level2_genes = NA
outTable$level3_genes = NA
for (cn in c("level1","level2","level3")){
  gn   = paste0(cn,"_genes")
  cols = table(levelInfo[,cn])
  cols = names(cols)[cols>2]
  for (cnn in cols){
	outTable[match(names(outGenesL[[cnn]]),outTable$final.tree.order),gn] = outGenesL[[cnn]]
  }
}
write.csv(outTable,paste0(outputFolder,"Table_S2_markerGenes.csv"))


#-----------------------------------------------------------------------------------------------
print("Plot these genes")

# Class genes defined using manual selection and NS forest.
classRF <- c("GAD1","GAD2","ADARB2","CXCL14","LOC105375415","FBXL7","SV2C","VIP","LHX6","NXPH1","PLCH1","SST","SOX6","PVALB","TAC1","SLC17A7","SATB2","PDE7B","FAM19A2","GRIK3","PDZRN4","SLC1A3")

# The specific genes defined above, separated into inhibitory and other.
isI <- substr(outTable$final.cluster,1,3)=="Inh"
plotGenesI <- plotGenesEN <- classRF
for (i in which(isI)){
 gnTmp      <- strsplit(outTable[i,"specific_genes"],", ")[[1]]
 plotGenesI <- c(plotGenesI, gnTmp)
}
for (i in which(!isI)){
 gnTmp      <- strsplit(outTable[i,"specific_genes"],", ")[[1]]
 plotGenesEN <- c(plotGenesEN, gnTmp)
}

all_clusters <- 1:75
markers_plotI <- group_heatmap_plot(data_source = dataFolder,
                   genes = plotGenesI,
                   group_by = "cluster",
                   clusters = all_clusters,
                   calculation = "median",
                   labelheight = 7,
                   showcounts = F)
markers_plotEN <- group_heatmap_plot(data_source = dataFolder,
                   genes = plotGenesEN,
                   group_by = "cluster",
                   clusters = all_clusters,
                   calculation = "median",
                   labelheight = 8,
                   showcounts = F)
				   
save_plot(paste0(outputFolder,"FigED4_top5_markers_Inh.pdf"),
          markers_plotI,base_height = 17, base_width = 8)
save_plot(paste0(outputFolder,"FigED4_top5_markers_ExNN.pdf"),
          markers_plotEN,base_height = 15, base_width = 8)
#
