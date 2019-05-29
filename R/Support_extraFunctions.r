plotTSNE_using_hicat <- function(norm.dat, anno, seed=123, useClusterInfoForDexGenes = FALSE, 
  n.markers=50, fn.size=3, cex=1, ...){

  library(ggplot2)
  library(M3Drop)
  library(scrattch.hicat)

  anno      <- anno[match(colnames(norm.dat),anno$sample_id),]
  ref.cl.df <- as.data.frame(unique(anno[, c("cluster_label", "cluster_id", "cluster_color")]))
  ref.cl.df <- ref.cl.df[order(ref.cl.df$cluster_id),]
  row.names(ref.cl.df) <- ref.cl.df$cluster_id
  ref.cl    <- droplevels(setNames(factor(anno$cluster_id), anno$sample_id))

  #### Select markers ####
  set.seed(seed)
  if((!useClusterInfoForDexGenes)|(length(levels(ref.cl))<2)){
    select.markers <- BrenneckeGetVariableGenes(2^norm.dat-1,suppress.plot=TRUE)
  } else {
    select.markers <- select_markers(norm.dat, ref.cl, n.markers=n.markers)$markers
  }

  #### Plot tSNE ####
  set.seed(seed+1)
  tsne.result <- plot_tsne_cl(norm.dat, select.markers, ref.cl, ref.cl.df, fn.size=fn.size, cex=cex, ...)
  tsne.result$g

}


findFromGroups <- function(datExpr,groupVector,fn="mean"){
  groups   = names(table(groupVector))
  fn       = match.fun(fn)
  datMeans = matrix(0,nrow=dim(datExpr)[2],ncol=length(groups))
  for (i in 1:length(groups)){
    datIn = datExpr[groupVector==groups[i],]
    if (is.null(dim(datIn)[1])) { datMeans[,i] = as.numeric(datIn)
    } else { datMeans[,i] = as.numeric(apply(datIn,2,fn)) }
  };    colnames(datMeans)  = groups;
  rownames(datMeans) = colnames(datExpr)
  return(datMeans)
}


error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
 if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
 stop("vectors must be same length")
 arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}


errorBarPlot <- function(vals,sampleType,col=standardColors(),legend=TRUE,elwd=2,ylim=NA,xlim=NA,length=0.1,...){
 if(is.null(dim(vals))) vals = cbind(vals,vals)
 yy <- t(findFromGroups(vals,sampleType))
 col = col[1:dim(yy)[1]]
 ee <- t(findFromGroups(vals,sampleType,sd))
 if(is.na(ylim[1])) ylim = c(0,max(ee+yy))
 if(is.na(xlim[1])) xlim = c(0,((dim(ee)[2]+1)*dim(ee)[1])*1.4)
 barx <- barplot(yy, beside=TRUE,col=col,legend.text=legend,ylim=ylim,xlim=xlim,...)
 error.bar(barx,yy,ee,lwd=elwd,length=length)
}

mouse2human2 <- function (mouse, m2h){
 # Convert mouse to human symbols
 rownames(m2h) = m2h$Mou
 noHumn = which(!(mouse%in%m2h$Mouse_Symbol))
 humn = which((mouse%in%m2h$Mouse_Symbol))
 mouse[humn] = as.character(m2h[mouse[humn],1])
 mouse[noHumn] = toupper(mouse[noHumn])
 return(mouse)
}

t.test.l <- function(x){
  l  = length(x)
  tt = t.test(x[1:(l/2)],x[(l/2+1):l],paired=FALSE)
  out = c(tt$est,tt$stat,tt$p.val)
  if(is.na(out[2])) out[2] = 0
  if(is.na(out[3])) out[3] = 1
  return(out)
}

getAnovaPvalforApply <- function(x,varLabels,varWeights=NULL){
  anovadat  = as.data.frame(cbind(varLabels,x))
  aov.out   = summary(aov(as.numeric(anovadat[,2])~anovadat[,1],data=anovadat,weights=varWeights))  
  return(aov.out[[1]]$'Pr(>F)'[1])
}



meanEx <- function(x) {if(sum(x)==0) return(0); return(mean(x[x>0]));}


t.test.l.paired <- function(x){
  l  = length(x)
  tt = t.test(x[1:(l/2)],x[(l/2+1):l],paired=TRUE)
  out = c(tt$est,tt$stat,tt$p.val)
  if(is.na(out[2])) out[2] = 0
  if(is.na(out[3])) out[3] = 1
  return(out)
}

  
calc_beta <- function(y, spec.exp = 2) {
  d1 <- as.matrix(dist(y))
  eps1 <- 1e-10
  # Marker score is combination of specificity and sparsity
  score1 <- sum(d1^spec.exp) / (sum(d1) + eps1)
  return(score1)
}

#####################################################################
# FUNCTIONS FOR BUILDING AND PLOTTING THE TREE ARE BELOW

getDend <- function(input,distFun = function(x) return(as.dist(1-cor(x)))){
 distCor  = distFun(input) 
 avgClust = hclust(distCor,method="average")
 dend = as.dendrogram(avgClust)
 dend = labelDend(dend)[[1]]
}

labelDend <- function(dend,n=1)
  {  
    if(is.null(attr(dend,"label"))){
      attr(dend, "label") =paste0("n",n)
      n= n +1
    }
    if(length(dend)>1){
      for(i in 1:length(dend)){
        tmp = labelDend(dend[[i]], n)
        dend[[i]] = tmp[[1]]
        n = tmp[[2]]
      }
    }
    return(list(dend, n))
  }

reorder.dend <- function(dend, l.rank,verbose=FALSE)
  {
    tmp.dend = dend
    sc=sapply(1:length(dend), function(i){
      l = dend[[i]] %>% labels
      mean(l.rank[dend[[i]] %>% labels])
    })
    ord = order(sc)
	if(verbose){
      print(sc)
	  print(ord)
    }
	if(length(dend)>1){
      for(i in 1:length(dend)){
        if(ord[i]!=i){
          dend[[i]]= tmp.dend[[ord[i]]]
        }
        if(length(dend[[i]])>1){
          dend[[i]]=reorder.dend(dend[[i]],l.rank)
        }
      }
    }
    return(dend)
  }

  
  
updateAndOrderClusters <- function(sampleInfo, # Subsetted Samp.dat file for samples from non-outlier clusters
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
  {
  
  ## Define clusterInfo variable to store cluster-level info
  kpColumns = unique(c("cluster_id","cluster_label","cluster_color",classNameColumn,matchNameColumn,regionNameColumn,otherColumns))
  clusterInfo = t(sampleInfo[,kpColumns])
  colnames(clusterInfo) = clusterInfo["cluster_label",]
  clusterInfo = clusterInfo[,unique(colnames(clusterInfo))]
  clusterInfo = t(clusterInfo)
  clusterInfo = as.data.frame(clusterInfo)
  clusterInfo$old_cluster_label = clusterInfo$cluster_label 
  rownames(clusterInfo) = 1:dim(clusterInfo)[1]
  if(!is.null(classLevels)) clusterInfo[,classNameColumn] = factor(clusterInfo[,classNameColumn],levels=classLevels)
  classLabel = clusterInfo[,classNameColumn]
  if(is.factor(classLabel)) classLabel = droplevels(classLabel)

  
  ## Add additional cluster information for renaming
  cl3 = sampleInfo[,"cluster_id"]
  names(cl3) <- sampleInfo$sample_id
  cl3 = factor(cl3,levels=as.numeric(clusterInfo$cluster_id))
  
  if(is.null(newColorNameColumn)) newColorNameColumn = "cluster_color"
  colorVec = as.character(tapply(names(cl3), cl3, function(x) { 
   col = as.factor(sampleInfo[,newColorNameColumn])
   names(col) = sampleInfo$sample_id
   return(names(sort(-table(col[x])))[1])
  }))
  clusterInfo$cluster_color = colorVec
    
  if(!is.null(regionNameColumn)){
   regionVec = as.character(tapply(names(cl3), cl3, function(x) { 
    rg = as.factor(sampleInfo[,regionNameColumn])
	names(rg) = sampleInfo$sample_id
	rg = rg[x]
    rg = table(rg)/table(sampleInfo[,regionNameColumn])
    rg = -sort(-round(100*rg/sum(rg)))[1]
    return(paste(names(rg),rg,sep="~"))
   }))
   clusterInfo$region = regionVec
  }
  clusterInfo$layer = 0
  if(!is.null(layerNameColumn)){
   clLayer <- getLayer(sampleInfo[,layerNameColumn])
   names(clLayer) = names(cl3)
   layerVec <- as.character(tapply(names(cl3), cl3, function(x) {
    lyy    = factor(clLayer)[x] 
    layTab = cbind(as.numeric(names(table(lyy))),table(lyy),table(clLayer))
	lyy    = clLayer[x] 
    #ly     = signif(mean(lyy),4)  # Only use 1 decimal
	ly     = signif(mean(lyy),2)
    ly     = paste(ly,paste(rep(0,3-nchar(ly)),collapse="",sep=""),sep="")
    ly     = gsub("00",".0",ly)
    ly     = paste("L",ly,sep="")
	
	return(ly)
   }))
   clusterInfo$layer = as.numeric(gsub("L","",layerVec))
  }
  if(!is.null(matchNameColumn)){
   matchVec <- as.character(tapply(names(cl3), cl3, function(x) {
    y  = is.element(sampleInfo$sample_id,x)
    nm = -sort(-table(sampleInfo[y,matchNameColumn]))
    return(names(nm)[1])
   }))
   clusterInfo$topMatch = matchVec
  }

  ## Rename the clusters based on the above info
  for (i in 1:dim(clusterInfo)[1]){
    id  = as.character(as.numeric(clusterInfo[i,"cluster_id"]))
	id2 = paste0("cl",id)
	broad = ifelse(is.na(classNameColumn),"",substr(clusterInfo[i,classNameColumn],1,1))
    cn  = paste0(broad,sum(cl3==id))
	lab = paste(id2,cn,sep=sep)
    lab = ifelse(is.null(topGenes),lab,paste(lab,topGenes[as.character(id)],sep=sep))
    lab = ifelse(is.null(matchNameColumn),lab,paste(lab,clusterInfo$topMatch[i],sep=sep))
    lab = ifelse(is.null(regionNameColumn),lab,paste(lab,clusterInfo$region[i],sep=sep))
	lab = ifelse(is.null(layerNameColumn),lab,paste(lab,layerVec[i],sep=sep))
	lab = gsub("-","~",lab)               # To avoid shiny crashing
	clusterInfo[i,"cluster_label"] = lab
  }

  ## Determine a new optimal order based on excitatory layer followed by mapped mouse type
  classLabel2 = factor(classLabel,levels=unique(classLabel))
  tmpLayer = clusterInfo[,"layer"]
  tmpLayer[classLabel2=="inh"] = tmpLayer[classLabel2=="inh"]-100 #0
  tmpLayer[classLabel2=="glia"] = tmpLayer[classLabel2=="glia"] +100 #10
  tmpMouse = tmpLayer
  if(!is.null(matchNameColumn)){
   tmpMouse = as.character(clusterInfo[,"topMatch"])
   tmpMouse = gsub("Vip","NVip",tmpMouse)
   tmpMouse = gsub("Smad","ASmad",tmpMouse)
  }
  ordNew   = order(tmpLayer,tmpMouse,clusterInfo$cluster_id)

  clusterInfo  = clusterInfo[ordNew,]
  rownames(clusterInfo) <- clusterInfo$lrank <- 1:dim(clusterInfo)[1]
  clusterInfo$cluster_id = as.character(as.numeric(clusterInfo$cluster_id))
  
  ## Return clusterInfo 
  return(clusterInfo)
}


# Newer function
getTopMarkersByPropNew <- function(propExpr, medianExpr, propDiff = 0, propMin=0.5, medianFC = 1, 
  excludeGenes = NULL, sortByMedian=TRUE){
  specGenes = rep("none",dim(propExpr)[2])
  names(specGenes) = colnames(propExpr)
   propSort = t(apply(propExpr,1,function(x) return(-sort(-x))))
  propWhich= t(apply(propExpr,1,function(x,y) return(y[order(-x)]),colnames(propExpr)))
  medianDif= apply(cbind(as.numeric(propWhich[,1]),medianExpr),1,function(x,y) {
   wIn = y==as.character(x[1])
   mIn = x[2:length(x)][wIn]
   mOut= max(x[2:length(x)][!wIn])
   return(mIn-mOut)
  }, colnames(propExpr))
  keepProp = (propSort[,1]>=propMin)&((propSort[,1]-propSort[,2])>propDiff)&(medianDif>=medianFC)&(!is.element(rownames(propExpr),excludeGenes))
  propSort = propSort[keepProp,]
  propWhich= propWhich[keepProp,]
  ord      = order(-medianDif[keepProp]*ifelse(sortByMedian,1,0), propSort[,2]-propSort[,1])
  propSort = propSort[ord,]
  propWhich= propWhich[ord,]
  while(sum(keepProp)>1){
   keepProp = !is.element(propWhich[,1],names(specGenes)[specGenes!="none"])
   if(sum(keepProp)<=1) break
   tmp = propWhich[keepProp,]
   specGenes[tmp[1,1]] = rownames(tmp)[1]
  }
  return(specGenes)
}

#newCluster10 = renameClusters(sampleInfo,clusterInfo, propExpr[keepGenes,], medianExpr[keepGenes,], 
#    propDiff = 0, propMin=0.5, medianFC = 1, propLayer=0.1,excludeGenes = excludeGenes,
#	majorGenes=majorGenes, majorLabels=majorLabels, broadGenes=broadGenes)
	

renameClusters <- function(sampleInfo,clusterInfo, propExpr, medianExpr, propDiff = 0, propMin=0.5, medianFC = 1, 
                           excludeGenes = NULL, majorGenes=c("GAD1","SLC17A7","SLC1A3"), majorLabels= majorGenes, 
						   broadGenes = majorGenes, propLayer=0.3, layerNameColumn="layer_label", 
						   getLayer = function(x) return(as.numeric(substr(as.character(x),2,2)))){
  # Layer determination
  clLayer <- getLayer(sampleInfo[,layerNameColumn])
  cl3 = sampleInfo[,"cluster_id"]
  names(cl3) <- sampleInfo$sample_id
  cl3 = factor(cl3,levels=sort(as.numeric(clusterInfo$cluster_id)))
  names(clLayer) = names(cl3)
  layerVec <- (tapply(names(cl3), cl3, function(x) {
    lyy    = factor(clLayer)[x] 
    layTab = cbind(as.numeric(names(table(lyy))),table(lyy),table(clLayer))
	return(((layTab[,2]/layTab[,3])/sum(layTab[,2]/layTab[,3])))  # replace max with sum?
  }))
  rn       = names(layerVec)
  layerVec = matrix(unlist(layerVec), ncol = 6, byrow = TRUE)
  rownames(layerVec) = rn
  colnames(layerVec) = 1:6
  layLab = apply(layerVec,1,function(x,y){
     z = as.numeric(colnames(layerVec)[x>=y])
	 if(length(z)==1) return(z)
	 return(paste(range(z),collapse="-"))
  }, propLayer)
  layLab = paste0("L",layLab)

  # Genes determination
  majorLab   = majorGenes[apply(propExpr[majorGenes,],2,which.max)]
  names(majorLabels) <- majorGenes
  
  broadLab   = broadGenes[apply(propExpr[broadGenes,],2,which.max)]
  broadProp  = apply(propExpr[broadGenes,],2,max)
  broadLab[broadProp<propMin] = majorLab[broadProp<propMin]
  names(broadLab) <- colnames(propExpr)
  
  betaScore  = getBetaScore(propExpr)
  kpGn       = rep(TRUE,dim(propExpr)[1])  #betaScore>=minBeta
  specGenes  = getTopMarkersByPropNew(propExpr=propExpr[kpGn,], medianExpr=medianExpr[kpGn,], propDiff = propDiff, propMin=propMin, 
               medianFC = medianFC, excludeGenes = excludeGenes)
  specGenes0 = getTopMarkersByPropNew(propExpr=propExpr[kpGn,], medianExpr=medianExpr[kpGn,], propDiff = 0, propMin=propMin, 
               medianFC = 0, excludeGenes = excludeGenes)
  for (s in colnames(propExpr)[(specGenes=="none")]){
    kp = propExpr[broadLab[s],]>=propMin
	specGenesTmp = getTopMarkersByPropNew(propExpr=propExpr[kpGn,kp], medianExpr=medianExpr[kpGn,kp], propDiff = propDiff, propMin=propMin, 
               medianFC = medianFC, excludeGenes = excludeGenes)
	if(specGenesTmp[s]!="none") specGenes[s] = specGenesTmp[s]
    if((specGenes0[s]!="none")&(specGenes[s]=="none")){
	  kp = (propExpr[specGenes0[s],] > quantile(propExpr[specGenes0[s],],0.8)) & (broadLab==broadLab[s])
	  specGenesTmp = getTopMarkersByPropNew(propExpr=propExpr[kpGn,kp], medianExpr=medianExpr[kpGn,kp], propDiff = propDiff, propMin=propMin, 
               medianFC = medianFC, excludeGenes = excludeGenes)
	  specGenes[s] = specGenes0[s]
	  if (length(grep(broadLab[s],broadGenes))==0) paste(specGenes[s],specGenesTmp[s])
	}
  }
  
  # Name construction
  nss       = names(specGenes)
  specGenes = paste0(" ",specGenes)
  specGenes[is.element(broadLab,intersect(majorGenes,broadGenes))] = ""
  clNames   = paste(majorLabels[majorLab]," ",layLab," ",broadLab,specGenes,sep="")
  names(clNames) = nss
  clNames = gsub(" none","",clNames)
  
  # Ensure that no excitatory clusters are included in layer 1
  clNames2 = clNames
  clNames = gsub("Exc L1","Exc L2",clNames)
  clNames = gsub("L2-2","L2",clNames)
  
  return(clNames)
}

# This is nearly identical to getTopMarkersByPropNew

getTopMarkersByProp <- function(...) suppressWarnings(getTopMarkersByProp2(...)) # Warnings are not useful from this function

getTopMarkersByProp2 <- function(propExpr,n=1,excludeGenes = NULL,medianExpr=NA,fcThresh=0.5,minProp=0){
 prop   = propExpr[!is.element(rownames(propExpr),excludeGenes),]
 cn     = colnames(prop)
 rn     = rownames(prop)
 wmProp = cn[apply(prop,1,function(x) return(which.max(x)[1]))]
 maxProp= rowMax(prop)
 dfProp = apply(prop,1,function(x) return(max(x)-max(x[x<max(x)])))
 dfProp[dfProp==Inf] = 0
 nmProp = apply(prop,1,function(x) return(sum(x==max(x))))
 dfMed  = dfProp*0+1
 kpI <- kp <- rep(TRUE,length(dfProp))&(nmProp==1)&(maxProp>=minProp)
 if(!is.na(medianExpr[1])){
  med   = medianExpr[!is.element(rownames(medianExpr),excludeGenes),]
  med   = med[rn,]
  wmMed = colnames(med)[apply(med,1,function(x) return(which.max(x)[1]))]
  dfMed = apply(med,1,function(x) return(max(x)-max(x[x<max(x)])))
  kp    = (wmProp==wmMed)&(dfMed>=fcThresh)&kpI&(maxProp>=minProp)
 }
 markCount = table(factor(wmProp[kp],levels=colnames(propExpr)))
 kpVal = kp+kpI
 outGenes = matrix(nrow=length(cn),ncol=n)
 colnames(outGenes) = paste0("Gene_",1:n)
 rownames(outGenes) = cn
 for (cc in cn){
  kk  = wmProp==cc
  ord = order(-kpVal[kk],-dfMed[kk])  #(dfMed[kk]*dfProp[kk]))
  outGenes[cc,] = rn[kk][ord][1:n]
 }
 nm = rownames(outGenes)
 outGenes = as.character(outGenes)
 names(outGenes) = nm
 return(outGenes)
}



updateSampDat <- function(Samp.dat,clusterInfo){
  lab = as.character(Samp.dat$cluster_label)
  id  = as.numeric(Samp.dat$cluster_id)
  for (i in 1:dim(clusterInfo)[1]){
    kpId = as.numeric(clusterInfo[i,"cluster_id"])
	lab[id==kpId] = clusterInfo[i,"cluster_label"]
  }  
  Samp.dat$cluster_label = lab
  return(Samp.dat)
}


getBetaScore <- function(propExpr,returnScore=TRUE,spec.exp = 2){
  # Wrapper for calc_beta
  # propExpr = proportions of cells in a given cluster with CPM/FPKM > 1 (or 0, HCT uses 1)
  betaScore <- apply(propExpr, 1, calc_beta)
  betaScore[is.na(betaScore)] <- 0
  if(returnScore) return(betaScore)
  scoreRank = rank(-betaScore)
  return(scoreRank)
}

rankCompleteDendPlot <- function(input=NULL,l.rank=NULL,dend=NULL,label_color=NULL,node_size=3,
    main="Tree",distFun = function(x) return(dist(1-(1+cor(x))/2)),...){
 if(is.null(dend)) dend = getRankedDend(input,l.rank,distFun)
 plot_dend(dend,node_size=node_size,label_color=label_color,main=main)
}

# Allows plot_dend to work properly in a for loop
plot_dend <- function(...) print(plot_dend2(...))
  
plot_dend2 <- function(dend, dendro_data=NULL,node_size=1,r=NULL,label_color=NULL,main="",rMin=-0.6)  # r=c(-0.1,1)
  {
    require(dendextend)
    require(ggplot2)
    if(is.null(dendro_data)){
      dendro_data = as.ggdend(dend)
      dendro_data$nodes$label =get_nodes_attr(dend, "label")
      dendro_data$nodes = dendro_data$nodes[is.na(dendro_data$nodes$leaf),]
    }
    node_data = dendro_data$nodes
    label_data <- dendro_data$labels
    segment_data <- dendro_data$segments
    if(is.null(node_data$node_color)){
      node_data$node_color="black"
    }
	if(!is.null(label_color)){
      label_data$col=label_color
    }
	rMax = max(node_data$height)*1.05
	main = paste(main,"- Max height =",signif(rMax,3))  # Add/update the title
	if(is.null(r)) r=c(rMin*rMax,rMax)  # Dinamically update the height
    ggplot() + ggtitle(main) + 
      geom_text(data = node_data, aes(x = x, y = y, label = label,color=node_color),size=node_size,vjust = 1) +
        geom_segment(data = segment_data,
                     aes(x=x,xend=xend,y=y,yend=yend), color="gray50") +
                       geom_text(data = label_data, aes(x = x, y = -0.01, label = label, color = col),size=node_size,angle = 90, hjust = 1) +
                           scale_color_identity() +
                             theme_dendro() +
                                scale_y_continuous(limits = r)
    
  }


plot_annotation_comparison <- function (cl, ref.cl, anno, title="plot") 
{
    if(!is.factor(ref.cl)) ref.cl <- as.factor(ref.cl)
    tb <- table(cl = cl, ref.cl = ref.cl)
	tb <- round(100*tb/rowSums(tb))
	
    tb.df <- as.data.frame(tb)
    tb.df <- tb.df[tb.df$Freq > 0, ]
	tb.df$Color <- anno$cluster_color[match(tb.df$cl,anno$cluster_label)]
    g <- ggplot(tb.df, aes(x = cl, y = ref.cl)) + 
         geom_point(aes(size = Freq, color = Color)) + 
         theme(axis.text.x = element_text(vjust = 0.1, 
            hjust = 0.2, angle = 90, size = 7), axis.text.y = element_text(size = 6)) + 
         scale_size(range = c(0, 3)) +
		 ggtitle(title) +
		 scale_color_identity() +
		 theme(legend.position = "none") 
    g
}