# This code performs the analysis of the mFISH data in the mansucript and plots the results

# FIRST NEED TO SET THE WORKING DIRECTORY TO THE MAIN FOLDER
# This is the folder where "Start.RData" is located


#-----------------------------------------------------------------------------------------------
print("Read in the data and load relevant libraries")

outputFolder  = "output/"
scriptsFolder = "R/"
inputFolder   = "data/"
dataFolder    = "MTG/"

library(dplyr);
library(feather)
library(pdist)
library(WGCNA)
library(gplots)
library(mfishtools)
library(dendextend)
library(tidyverse)
library(gridExtra)
library(viridis)
options(stringsAsFactors=FALSE)

anno     <- read_feather(paste(dataFolder,"anno.feather",sep="")) 
Expr.dat <- feather(paste(dataFolder,"data.feather",sep="")) 
Expr.dat <- Expr.dat[match(anno$sample_id,Expr.dat$sample_id),] # Make sure the expression matches the sample information
datIn    <- as.matrix(Expr.dat[,names(Expr.dat)!="sample_id"])
rownames(datIn) <- Expr.dat$sample_id
datIn    <- t(datIn)

load(paste0(dataFolder,"/dend.rda")) 
load(paste0(dataFolder,"clusterInfo.rda")) 
rownames(clusterInfo) <- clusterInfo$cluster_label

##################################################################

print("Find the median expression in each non-outlier cluster")

kpSamp     <- 1:dim(anno)[1]
cl3        <- factor(anno$cluster_label,levels = sort(unique(anno$cluster_label)))
names(cl3) <- colnames(datIn)
clustersF  <- droplevels(cl3[kpSamp])

normDat    <- log2(datIn[,kpSamp]+1)
medianExpr <- do.call("cbind", tapply(names(clustersF), clustersF, function(x) apply(normDat[,x],1,median))) 
medianExpr <- medianExpr[,levels(clustersF)]
rownames(medianExpr) <- rownames(normDat) 


##################################################################
print("Determine the predicted level of mappability using this gene panel based on FACS data.")

## THESE ARE FILES WITH RESULTS FROM THE THREE EXPERIMENTS
nm = c("smFISH_Human_L2-L3_2018_04_output","smFISH_Human_L2-L3_Replicate2_output","smFISH_Human_L2-L3_Replicate3_output")

fn  = paste0(inputFolder,nm[1],".csv")
fishIn = read.csv(fn,row.names=1)
gns = intersect(rownames(medianExpr),colnames(fishIn))

kpSamp2 = colSums(normDat[gns,]>0)>1
facsCor = corTreeMapping(medianExpr,normDat[gns,kpSamp2])
facsCl  <- rownames(facsCor)[apply(facsCor,2,which.max)]
facsConf = table(factor(as.character(clustersF[kpSamp2]),levels=intersect(labels(dend),clustersF[kpSamp2])),
  factor(facsCl,levels=intersect(labels(dend),facsCl)))
facsFrac <- facsConf / rowSums(facsConf)

# Which leaves have which nodes?
has_any_labels <- function(sub_dend, the_labels) any(labels(sub_dend) %in% the_labels)
node_labels <- NULL
for (lab in labels(dend))
  node_labels <- cbind(node_labels,noded_with_condition(dend, has_any_labels,the_labels=lab))
rownames(node_labels) <- get_nodes_attr(dend,"label")
colnames(node_labels) <- labels(dend)

# Which clusters agree at the node level?
agreeNodes = apply(cbind(facsCl,as.character(clustersF[kpSamp2])),1, function(lab,node_labels){
  rowSums(node_labels[,lab])==2
},node_labels)
colnames(agreeNodes) = as.character(clustersF[kpSamp2])

# Which clusters are in each nodes?
isInNodes = t(apply(node_labels,1, function(node,cl,dend){
  is.element(cl,labels(dend)[node])
},as.character(clustersF[kpSamp2]),dend))
colnames(isInNodes) = as.character(clustersF[kpSamp2])

# For each node, what fraction of cells match?
fracAgree = rowSums(agreeNodes) / rowSums(isInNodes)
pdf(paste0(outputFolder,"FACS_mapping_accuracy.pdf"),height=8,width=20)
par(mar=c(12,5,5,5))
dend %>% set("nodes_cex",0) %>% set("branches_col", "grey") %>% plot
text(get_nodes_xy(dend)[,1],get_nodes_xy(dend)[,2],round(fracAgree*100))
dev.off()

# Based on results above above, build a new vector of clusters for FACS comparison
classNew = c(rep("LAMP5- Inh",2),rep("LAMP5+ Inh",4),rep("LAMP5- Inh",39),labels(dend)[46:48],"L5-6",
  labels(dend)[50],rep("RORB+ L3-5",4),rep("RORB+ L4-5",5),labels(dend)[60],rep("L5-6",9),rep("Non-Neuronal",6))
names(classNew) = labels(dend)
medianClass = findFromGroups(t(medianExpr[gns,labels(dend)]),classNew,median)
medianClass = medianClass[,unique(classNew)[c(2,1,3,5,4,7,8,9,6,11)]]  # 9,10,6
   # Omit Exc L4−5 FEZF2 SCN4B, which seems to be a catch-all bucket, and which we know is VERY rare in human


##################################################################
print("Read in the data and order consistently.")

qprob  = 0.9    # For scaling mFISH to FACS
thresh = 0       # Set counts less than or equal to thresh to 0
fish   = list()
par(mfrow=c(3,3))
nm2 = NULL
for (i in 1:length(nm)){
  fn = paste0(inputFolder,nm[i],".csv")
  fishIn  <- fishTmp <- read.csv(fn,row.names=1)
  
  fishTmp$area_um2 <- fishTmp$area / 100
  fishTmp$total_density <- fishTmp$totalReads / fishTmp$area_um2
  filter = fishTmp$area_um2 > 100 & fishTmp$total_density > 0 & fishTmp$total_density < 1 & fishTmp[,"layerData"]==1
  isExc  = (fishTmp[,"SLC17A7"] >= quantile(fishTmp[filter,"SLC17A7"],0.95))
  isExc  = isExc & fishTmp$area_um2 > 100 & fishTmp$total_density > 0 & fishTmp$total_density < 1 
  
  fishIn$area[fishIn$area==0] = round(mean(fishIn$area))  # To avoid errors later
  fishDat <- as.matrix(t(fishIn[,intersect(colnames(fishIn),rownames(medianClass))]))
  fishDat[fishDat<=thresh] = 0 
	
  fishDat <- t(t(fishDat)/fishIn$area)*mean(fishIn$area) # Account for spot area in gene expression calculation
  fishDat <- log2(fishDat+1)   # Do log transform (It looks worse if we don't)

   plot(density(fishDat["SLC17A7",]))
   plot(density(fishDat["GAD2",]))
   plot(density(fishIn$area))
  
  # Only map to excitatory types found in layers 2 and 3
  isAll   <- colnames(medianClass)  
  is23    <- colnames(medianClass)[3:9]
  
  ###############################################################
  # Scale fish data to roughly match RNA-Seq values
  rn = rownames(fishDat)
  valFACS  = apply(medianClass[rn,is23],1,max)
  qs = apply(fishDat[,isExc],1,quantile,qprob,na.rm=TRUE)
  fishDat2 = fishDat
  for (r in rn) fishDat2[r,] = fishDat2[r,]*valFACS[r]/pmax(qs[r],1)  # Value may need to be adjusted on a case by case basis.
  fishScale=fishDat2
  
  # Scale all values to 0-1 before mapping
  medianScale = medianClass/apply(medianClass,1,max)
  ###############################################################
  
  fishCor <- corTreeMapping(medianScale,fishScale)
  fishCor[is.na(fishCor)] = 0
  fishCl  <- rownames(fishCor)[apply(fishCor,2,which.max)]
  fishType<- ifelse(isExc,"Excitatory","Other")
  fishMax <- apply(fishCor,2,max)
  fishSec <- apply(fishCor,2,function(x,y){
    v = which.max(x)[1]
	return(max(x[y!=y[v]]))
  },isAll)
  fish2 <- data.frame(cluster = fishCl, maxCor = fishMax, confidence = fishMax-fishSec, type = fishType)
  fish2 <- cbind(fish2,fishIn)
  fish2$x <- fishIn$lateralCoordinate
  fish2$y <- fishIn$layerCoordinate

  # Scale to (0,1)
  fish2$x = fish2$x-min(fish2$x)
  fish2$x = fish2$x/max(fish2$x)
  fish2$y = fish2$y-min(fish2$y)
  fish2$y = fish2$y/max(fish2$y)

  # Save results
  fish[[i]] = fish2
}

# Use class colors rather than cluster colors
cols = as.character(clusterInfo$cluster_color)
names(cols) <- rownames(clusterInfo)



#################################################################################33
## Cluster-based distribution plots - EXTENDED DATA FIGURE 9D

pdf(paste0(outputFolder,"FigED9d_mFISH_classDistribution.pdf"),height=8,width=16)
for (i in 1:length(fish)){
 fishh = fish[[i]][fish[[i]]$type=="Excitatory",]
 clusts = intersect(is23,fishh$cluster)
 par(mfrow=c(1,length(clusts)))
 for (j in clusts){
  kp = fishh$cluster==j
  plot(fishh$x[kp],fishh$y[kp],xlab="x",ylab=paste("y -",j),main="Human ISH",
    xlim=c(0,1),ylim=c(1,0),col="white")
  text(fishh$x[kp],fishh$y[kp],fishh$layerData[kp])
 }
 print(paste("Fraction of excitatory cells mapping to excitatory types:",signif(mean(is.element(fishh$cluster,is23)),3)))
}
dev.off()
#[1] "Fraction of excitatory cells mapping to excitatory types: 0.984"
#[1] "Fraction of excitatory cells mapping to excitatory types: 0.988"
#[1] "Fraction of excitatory cells mapping to excitatory types: 0.996"


#################################################################################33
## Heatmap

gnH = c("GAD2","SLC17A7","CUX2","LAMP5","CBLN2","PENK","COL5A2","CARTPT","RXFP1","confidence")
cap = 0.001
colorset = c("darkblue", "dodgerblue", "gray80", "orange", "orangered")
heat_colors <- colorRampPalette(colorset)(1001)

pdf(paste0(outputFolder,"mFISH_heatmap.pdf"),width=12,height=4)
for (i in 1:length(fish)){
 fishh = fish[[i]][fish[[i]]$type=="Excitatory",]
 datTmp = pmin(as.matrix(fishh[,gnH])/fishh[,"area"],cap)
 datTmp[,"confidence"] = sqrt(fishh[,"confidence"])*cap
 datTmp = datTmp[order(factor(fishh$cluster,levels = clusts),-fishh$x),]
 seps   = cumsum(table(fishh$cluster)[clusts])
 heatmap.2(t(datTmp),Rowv=FALSE,Colv=FALSE,dendrogram="none",trace="none",margins = c(3, 10),
  rowsep=9,colsep=seps,key=TRUE,main=paste("Human (by layer), counts/area capped at",cap),col=heat_colors)
 datTmp = pmin(as.matrix(fishh[,gnH])/fishh[,"area"],cap)
 datTmp[,"confidence"] = sqrt(fishh[,"confidence"])*cap
 datTmp = datTmp[order(factor(fishh$cluster,levels = clusts),fishh$confidence),]
 heatmap.2(t(datTmp),Rowv=FALSE,Colv=FALSE,dendrogram="none",trace="none",margins = c(3, 10),
  rowsep=9,colsep=seps,key=TRUE,main=paste("Human (by confidence), counts/area capped at",cap),col=heat_colors)
}
dev.off()



#################################################################################
## Comparison with RNA-Seq data

fishh     = rbind(fish[[1]],fish[[2]],fish[[3]])
fishh     = fishh[(fishh$type=="Excitatory")&(is.element(fishh$cluster,is23)),]
fishProps = table(fishh$cluster,fishh$layerData)[clusts,]
colnames(fishProps) = paste0("L",colnames(fishProps))
fishProps = fishProps[,-1]
kpL23     = is.element(anno$cluster_label,names(classNew)[is.element(classNew,clusts)])
rSeqProps = table(classNew[anno$cluster_label[kpL23]],paste0("L",anno$brain_subregion_id[kpL23]))[clusts,colnames(fishProps)]
cols2     = clusterInfo[clusts,"final.cluster.color"]
d1        = dim(fishProps)[1]
d2        = dim(fishProps)[2]

pdf(paste0(outputFolder,"mFISH_classByLayerBarplot_human.pdf"),width=9,height=9)
plot(fishProps,main="Human (HCR)",ylab="Layer",col=standardColors(),las=2)
plot(rSeqProps,main="Human (RNA-Seq)",ylab="Layer",col=standardColors(),las=2)
plot(t(fishProps),main="Human (HCR)",ylab="Layer",col=cols2,las=2)
plot(t(rSeqProps),main="Human (RNA-Seq)",ylab="Layer",col=cols2,las=2)

barplot(t(t(t(fishProps))/colSums(t(fishProps)))[d2:1,],main="Human (HCR)",ylab="Layer",col=standardColors(),las=2,legend=TRUE,ylim=c(0,2))
barplot(t(t(t(rSeqProps))/colSums(t(rSeqProps)))[d2:1,],main="Human (RNA-Seq)",ylab="Layer",col=standardColors(),las=2,legend=TRUE,ylim=c(0,2))
barplot(t(t(fishProps)/colSums(fishProps))[d1:1,],main="Human (HCR)",ylab="Layer",col=cols2,las=2,legend=TRUE,ylim=c(0,2))
barplot(t(t(rSeqProps)/colSums(rSeqProps))[d1:1,],main="Human (RNA-Seq)",ylab="Layer",col=cols2,las=2,legend=TRUE,ylim=c(0,2))
dev.off()


##########################################
print("Output the results")

fishOut = NULL
for (i in 1:length(fish))
  fishOut = rbind(fishOut,cbind(fish[[i]],nm[i]))
colnames(fishOut) = c(colnames(fish[[1]]),"Experiment")
write.csv(fishOut,paste0(outputFolder,"L23_FISH_output.csv"),row.names=FALSE)



##########################################
print("Re-read in the FISH data and scale spatial coordinates to experiment smFISH_Human_L2-L3_2018_04")

fish <- read_csv(file = paste0(outputFolder,"L23_FISH_output.csv"))

fix.exp <- which(fish$experimentName == "smFISH_Human_L2-L3_2018_04")
fish$layerCoordinate[fix.exp] <- fish$layerCoordinate[fix.exp] - min(fish$layerCoordinate[fix.exp])
for (exp1 in unique(fish$experimentName)) {
  lat1 <- fish$lateralCoordinate[fish$experimentName == exp1]
  fish$lateralCoordinate[fish$experimentName == exp1] <- lat1 - min(lat1)
}
fish$LateralPosition <- abs(fish$lateralCoordinate / 10)
fish$CorticalDepth <- -fish$layerCoordinate / 10
fish$area_um2 <- fish$area / 100
fish$total_density <- fish$totalReads / fish$area_um2


##########################################
print("Plot the FISH results for ALL clusters")


plot_fish <- function(dat, grey.pal.rep = 0) {
  dat.sort <- dat[order(dat$density), ]
  density.breaks <- c(0, 1, 3, 10, 30)
  ggplot(data = dat.sort) +
    aes(x = LateralPosition, y = CorticalDepth, 
        color = density, size = density) +
    geom_point() +
    scale_color_gradientn(colors = c(rep("grey90", grey.pal.rep), viridis(10)), 
                          trans = "sqrt", breaks = density.breaks,
                          name = expression(paste("Spots/100", mu, plain(m)^2))) +
    scale_size_area(breaks = density.breaks, name = "") +
    guides(color = guide_colorbar(order = 1),
           size = guide_legend(order = 0)) +
    ggtitle(first(dat$gene)) +
    xlab(expression(paste("Lateral position (", mu, "m)"))) +
    ylab(expression(paste("Cortical depth (", mu, "m)"))) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
}


fish.genes <- c("SLC17A7", "LAMP5", "PENK", 
                "CUX2", "COL5A2", "CARTPT",
                "CBLN2", "RXFP1", "GAD2")
fish.subset <- fish %>%
  filter(area_um2 > 100 & total_density > 0 & total_density < 1) %>%
  gather(gene, count, CUX2:SLC17A7) %>%
  mutate(gene = factor(gene, levels = fish.genes)) %>%  # Reorder plots
  group_by(experimentName, gene) %>%
  mutate(density = count/area_um2 * 100) # Probe density per 100um^2
  
g.all <- fish.subset %>% 
  do(plots = plot_fish(., grey.pal.rep = 1))


marrangeGrob(g.all$plots, nrow = 3, ncol = 3, top = "")

pdf(file = paste0(outputFolder,"FigED9cd_L23_fish_probes.pdf"), width = 10.5, height = 12)
marrangeGrob(g.all$plots, nrow = 3, ncol = 3, top = "")
dev.off()


##########################################
print("Plot the data only for the FREM3 cluster.")

frem3.genes <- c("LAMP5", "COL5A2")

fish.frem3 <- fish %>%
  filter(area_um2 > 100 & total_density > 0 & total_density < 1) %>%
  filter(SLC17A7 > 9.75) %>% 
  filter(cluster == "Exc L2-3 LINC00507 FREM3") %>% 
  gather(gene, count, CUX2:SLC17A7) %>%
  filter(gene %in% frem3.genes) %>% 
  mutate(gene = factor(gene, levels = frem3.genes)) %>%  # Reorder plots
  group_by(experimentName, gene) %>%
  mutate(density = count/area_um2 * 100) # Probe density per 100um^2

g.frem3 <-  fish.frem3 %>%
  group_by(experimentName, gene) %>% 
  do(plots = plot_fish(., grey.pal.rep = 1) +
       ylim(c(min(fish.subset$CorticalDepth), 
              max(fish.subset$CorticalDepth)))) # 


marrangeGrob(grobs = g.frem3$plots, nrow = 1, ncol = 1, top = "")

pdf(file = paste0(outputFolder,"Fig3c_L23_FREM3_fish_probes.pdf"), width = 3.5, height = 5)
marrangeGrob(g.frem3$plots, nrow = 1, ncol = 1, top = "")
dev.off()


