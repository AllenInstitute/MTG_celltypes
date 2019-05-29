# This code makes some scatterplots comparing neurosurgical vs. postmortem of matched types.

# FIRST NEED TO SET THE WORKING DIRECTORY TO THE MAIN FOLDER
# This is the folder where "Start.RData" is located


#-----------------------------------------------------------------------------------------------
print("Read in the data and load relevant libraries")
outputFolder  = "output/"
scriptsFolder = "R/"
dataFolder    = "MTG/"

# Load these libraries
library(feather)
options(stringsAsFactors=FALSE)

# Read in the data
anno       <- read_feather(paste(dataFolder,"anno.feather",sep="")) 
Expr.dat   <- feather(paste(dataFolder,"data.feather",sep=""))   # FPKM
Expr.dat   <- Expr.dat[match(anno$sample_id,Expr.dat$sample_id),] # Make sure the expression matches the sample information
datIn      <- as.matrix(Expr.dat[,names(Expr.dat)!="sample_id"])
rownames(datIn) <- Expr.dat$sample_id
datIn      <- t(datIn)
norm.dat   <- log2(datIn+1)

#-----------------------------------------------------------------------------------------------
print("Compare neurosurgical and postmortem data.")

isPM  <- is.element(anno$donor_label,c("H16.24.010","H200.1023","H200.1025","H200.1030"))
homs  <- list(
  SST = paste("Sst",1:5),
  Layer4 = "Exc L4/5 IT",
  Layer5IT = paste("Exc L5/6 IT",1:3)
)

pdf(paste0(outputFolder,"FigED2c_neurosurgicalVsPostmortem.pdf"))
for (nm in c("Layer5IT", "Layer4", "SST")){
  ns <- (!isPM)&is.element(anno$homology_cluster_label,homs[[nm]])
  pm <- isPM&is.element(anno$homology_cluster_label,homs[[nm]])
  nsDat <- rowMeans(norm.dat[,ns])
  pmDat <- rowMeans(norm.dat[,pm])
  
  plot(nsDat,pmDat,pch=19,cex=0.3,main=paste(nm,"- R =",signif(cor(nsDat,pmDat),2)),
       xlab=paste("Neurosurgical -",sum(ns)),ylab=paste("Postmortem -",sum(pm)))
}
dev.off()
