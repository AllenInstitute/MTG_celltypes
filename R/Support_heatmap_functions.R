######################################################################################
# CLUSTER HEATMAPS

group_heatmap_plot <- function (genes = c("Hspa8", "Snap25", "Gad2", "Vip"), clusters = 1:10, 
    group_by = "final", calculation = "mean", data_source = "internal", 
    normalize_rows = FALSE, logscale = T, fontsize = 7, labelheight = 25, 
    max_width = 10, showcounts = T, rotatecounts = F, maxval = "auto", 
    colorset = c("darkblue", "dodgerblue", "gray80", "orange", 
        "orangered")) 
{
    library(dplyr)
    library(ggplot2)
    genes <- rev(genes)
    if (data_source == "internal") {
        data <- get_internal_data(genes, group_by, clusters)
    }
    else if (is.list(data_source)) {
        data <- get_list_data(data_source, genes, group_by, clusters)
    }
    else if (grepl("\\\\.db$", data_source)) {
        data <- get_db_data(data_source, genes, group_by, clusters)
    }
    else if (file.exists(paste0(data_source, "/anno.feather"))) {
        data <- get_feather_data(data_source, genes, group_by, 
            clusters)
    }
    else {
        stop("Cannot identify data_source.")
    }
    data <- data %>% select(-xpos) %>% mutate(xpos = plot_id)
    genes <- sub("-", ".", genes)
    genes[grepl("^[0-9]", genes)] <- paste0("X", genes[grepl("^[0-9]", 
        genes)])
    names(data)[grepl("^[0-9]", genes)] <- paste0("X", names(data)[grepl("^[0-9]", 
        genes)])
    genes <- genes[genes %in% names(data)]
    ngenes <- length(genes)
    nclust <- length(unique(data$plot_id))
    nsamples <- nrow(data)
    header_labels <- build_header_labels(data = data, ngenes = ngenes, 
        nsamples = 1, nclust = nclust, labelheight = labelheight, 
        labeltype = "simple")
    heat_data <- data %>% select(plot_id, xpos)
    for (i in 1:length(genes)) {
        gene <- genes[i]
        if (calculation == "mean") {
            gene_func <- paste0("mean(", gene, ")")
        }
        else if (calculation == "trimmed_mean") {
            gene_func <- paste0("mean(", gene, ",trim = 0.25)")
        }
        else if (calculation == "percent") {
            gene_func <- paste0("sum(", gene, " > 0)/length(", 
                gene, ")")
        }
        else if (calculation == "median") {
            gene_func <- paste0("stats::median(", gene, ")")
        }
        gene_data <- data %>% select(one_of(c("plot_id", gene))) %>% 
            group_by(plot_id) %>% summarize_(result = gene_func)
        names(gene_data)[2] <- gene
        heat_data <- heat_data %>% left_join(gene_data, by = "plot_id")
    }
    heat_data <- unique(heat_data)
    max_vals <- heat_data %>% select(one_of(genes)) %>% summarise_each(funs(max)) %>% 
        unlist()
    max_labels <- data.frame(x = (nclust + 0.5) * 1.01, y = 1:ngenes + 
        0.5, label = sci_label(max_vals))
    max_header <- data.frame(x = (nclust + 0.5) * 1.01, y = ngenes + 
        1, label = "Max value")
    max_width <- nclust * (max_width/100)/(1 - max_width/100)
    if (logscale) {
        heat_data[genes] <- log10(heat_data[genes] + 1)
    }
    heat_colors <- colorRampPalette(colorset)(1001)
    if (maxval == "auto") {
        data_max <- max(unlist(heat_data[genes]))
    }
    else {
        data_max <- maxval
    }
    for (gene in genes) {
        if (normalize_rows == T) {
            heat_data[gene] <- heat_colors[unlist(round(heat_data[gene]/max(heat_data[gene]) * 
                1000 + 1, 0))]
        }
        else {
            color_pos <- unlist(round(heat_data[gene]/data_max * 
                1000 + 1, 0))
            color_pos[color_pos > 1001] <- 1001
            heat_data[gene] <- heat_colors[color_pos]
        }
    }
    label_y_size <- max(header_labels$ymax) - min(header_labels$ymin)
    cluster_data <- data %>% group_by(plot_label, plot_color, 
        plot_id) %>% summarise(cn = n()) %>% as.data.frame(stringsAsFactors = F) %>% 
        arrange(plot_id) %>% mutate(labely = ngenes + label_y_size * 
        0.05, cny = max(header_labels$ymax) - 0.1 * label_y_size, 
        xpos = plot_id)
    p <- ggplot() + scale_fill_identity() + scale_y_continuous("", 
        breaks = 1:length(genes) + 0.5, labels = genes, expand = c(0, 
            0)) + scale_x_continuous("", expand = c(0, 0)) + 
        theme_classic(fontsize) + theme(axis.text = element_text(size = rel(1)), 
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        legend.position = "none")
    for (i in 1:length(genes)) {
        p <- p + geom_rect(data = heat_data, aes_string(xmin = "xpos - 0.5", 
            xmax = "xpos + 0.5", ymin = i, ymax = i + 1, fill = genes[i]))
    }
    p <- p + geom_rect(data = header_labels, aes(xmin = xmin, 
        ymin = ymin, xmax = xmax, ymax = ymax, fill = color)) + 
        geom_text(data = header_labels, aes(x = (xmin + xmax)/2, 
            y = ymin + 0.05, label = label), angle = 90, vjust = 0.35, 
            hjust = 0, size = pt2mm(fontsize)) + geom_rect(aes(xmin = nclust + 
        0.5, xmax = nclust + 0.5 + max_width, ymin = 1, ymax = max(header_labels$ymax)), 
        fill = "white") + geom_text(data = max_header, aes(x = x, 
        y = y, label = label), angle = 90, hjust = 0, vjust = 0.35, 
        size = pt2mm(fontsize)) + geom_text(data = max_labels, 
        aes(x = x, y = y, label = label), hjust = 0, vjust = 0.35, 
        size = pt2mm(fontsize), parse = TRUE)
    if (showcounts) {
        if (rotatecounts) {
            p <- p + geom_text(data = cluster_data, aes(y = cny, 
                x = xpos, label = cn), angle = 90, vjust = 0.35, 
                hjust = 1, size = pt2mm(fontsize))
        }
        else {
            p <- p + geom_text(data = cluster_data, aes(y = cny, 
                x = xpos, label = cn), size = pt2mm(fontsize))
        }
    }
    return(p)
}




######################################################################################
# CELL BY CELL HEATMAPS

sample_heatmap_plot <- function(data_source, 
                                genes = c("Hspa8","Snap25","Gad2","Vip"), 
                                group_by = "cluster", 
                                groups = 1:10,
                                top_labels = "layer",
                                sample_n = 0,
                                scale_mode = "scale.log",
                                showall = F,
                                autorange = "auto", 
                                minrange = 0, 
                                maxrange = 10,
                                pfontsize = 14,
                                expand = F,
                                rotatelabels = F,
                                showlines = F,
                                showids = T) {
  
  library(feather)
  library(dplyr)
  library(ggplot2)
  
  data <- read_feather(file.path(data_source, "data.feather"), columns = c("sample_id",genes))
  anno <- read_feather(file.path(data_source, "anno.feather"))
  desc_table <- read_feather(file.path(data_source, "desc.feather"))
  
  primary <- list(base = group_by,
                  id = paste0(group_by,"_id"),
                  label = paste0(group_by,"_label"),
                  color = paste0(group_by,"_color"))
  secondary <- list(base = top_labels,
                    id = paste0(top_labels,"_id"),
                    label = paste0(top_labels,"_label"),
                    color = paste0(top_labels,"_color"))
  
  genes[grepl("^[0-9]",genes)] <- paste0("X",genes[grepl("^[0-9]",genes)])
  
  genes.df <- data

  ###############
  ## Filtering ##
  ###############
  
  # Join the annotation and genes data frames
  sub.df <- anno %>% 
    # Filter for the selected clusters
    filter_(paste0(primary$id, " %in% c(", paste(groups, collapse = ","), ")")) %>%
    left_join(genes.df)
  
  # Subsample data per primary group if Sample N is set to something other than 0
  if(sample_n > 0) {
    sub.df <- sub.df %>%
      group_by_(primary$id) %>%
      sample_n(sample_n, replace = T) %>%
      unique() %>%
      ungroup() %>%
      as.data.frame()
  }
  
  #############
  ## Scaling ##
  #############

  ## If the y-axis is plotted on a log scale, add 1 to the data values to plot data + 1
  if(scale_mode == "scale.log") {
    for(gene in genes) {
      sub.df[,gene] <- log10(sub.df[,gene] + 1)
    }
  }
  if(scale_mode == "scale.rel") {
    for(gene in genes) {
      sub.df[,gene] <- sub.df[,gene]/max(sub.df[,gene])
    }
  }
  if(scale_mode == "scale.log.rel") {
    for(gene in genes) {
      sub.df[,gene] <- log10(sub.df[,gene]+1)/log10((max(sub.df[,gene]+1)))
    }
  }
  
  #############
  ## Sorting ##
  #############
  
  cluster_order <- data.frame(clust = groups,
                              plot_order = 1:length(groups))
  
  names(cluster_order)[1] <- primary$id
  
  sub.df <- sub.df %>% 
    left_join(cluster_order, by = primary$id)
  
  sort.df <- sub.df %>% 
    arrange_(.dots = c("plot_order", secondary$id)) %>% 
    mutate(xpos = 1:nrow(sub.df))
  
  
  
  # Start buildplot
  genes <- sub("-", ".", genes)
  genes[grepl("^[0-9]",genes)] <- paste0("X",genes[grepl("^[0-9]",genes)])
  
  names(sort.df) <- sub("-",".",names(sort.df))
  
  colors <- colorRampPalette(c("darkblue", "white", "red"))(1001)
  
  if(autorange == "auto") {
    min.val <- 0
    max.val <- max(unlist(sort.df[ ,genes]))
  } else if (autorange == "manual") {
    min.val <- as.numeric(minrange)
    max.val <- as.numeric(maxrange)
  }
  
  ## Convert data to geom_rect() compatible table
  plot.df <- data.frame(xmin=numeric(),xmax=numeric(),ymin=numeric(),ymax=numeric(),fill=character())
  
  for(i in 1:length(genes)) {
    
    fill_ids <- unlist(round( (sort.df[,genes[i]] - min.val) / (max.val - min.val) * 1000 ) + 1)
    fill_ids[fill_ids < 1] <- 1
    fill_ids[fill_ids > 1001] <- 1001
    
    gene.plot <- data.frame(xmin = sort.df$xpos - 1,
                            xmax = sort.df$xpos,
                            ymin = length(genes) - i,
                            ymax = length(genes) - i + 1,
                            fill = colors[fill_ids])
    
    plot.df <- rbind(plot.df, gene.plot)
    
  }
  
  primary.plot <- data.frame(xmin = sort.df$xpos - 1,
                             xmax = sort.df$xpos,
                             ymin = -0.5,
                             ymax = 0,
                             fill = unlist(sort.df[ ,primary$color]))
  
  
  
  plot.df <- rbind(plot.df, primary.plot)
  
  ## add additional secondary color bars
  all.desc <- desc_table
  primary.name   <- all.desc$name[all.desc$base == primary$base]
  secondary.name <- all.desc$name[all.desc$base %in% secondary$base]
  
  primary.desc <- all.desc[all.desc$base == primary$base,]
  secondary.desc <- all.desc[all.desc$base %in% secondary$base,]
  other.desc <- all.desc[!all.desc$base %in% c(primary$base,secondary$base),]
  
  sec.color <- paste0(secondary$base, "_color")
  anno.color <- paste0(other.desc$base, "_color")
  anno_y_labels <- data.frame(breaks = numeric(), labels = character())
  
  if(showall) {
    
    # scale the plot so it's evenly divided between annotations and genes
    #s <- length(genes)/nrow(other.desc)*0.75
    s <- 1
    
    for(j in 1:length(secondary$base)) {
      anno.plot <- data.frame(xmin = sort.df$xpos - 1,
                              xmax = sort.df$xpos,
                              ymin = length(genes) + (j - 1) * s,
                              ymax = length(genes) + (j - 1) * s + s,
                              fill = unlist(sort.df[ ,sec.color[j]]))
      plot.df <- rbind(plot.df, anno.plot)
      
      anno_y <- data.frame(breaks = length(genes) + (j - 1)*s + 0.5*s,
                           labels = secondary.desc$name[secondary.desc$base == secondary$base[j]])
      
      anno_y_labels <- rbind(anno_y_labels,anno_y)
    }
    
    sec_top <- length(secondary$base)
    
    for(j in 1:nrow(other.desc)) {
      anno.plot <- data.frame(xmin = sort.df$xpos - 1,
                              xmax = sort.df$xpos,
                              ymin = length(genes) + (j + sec_top - 1) * s,
                              ymax = length(genes) + (j + sec_top - 1) * s + s,
                              fill = unlist(sort.df[ ,anno.color[j]]))
      plot.df <- rbind(plot.df, anno.plot)
      
      anno_y <- data.frame(breaks = length(genes) + (j + sec_top - 1)*s + 0.5*s,
                           labels = other.desc$name[j])
      
      anno_y_labels <- rbind(anno_y_labels,anno_y)
    }
  } else {
    
    if(length(secondary) > 0) {
      
      for(j in 1:length(secondary$base)) {
        anno.plot <- data.frame(xmin = sort.df$xpos - 1,
                                xmax = sort.df$xpos,
                                ymin = length(genes) + (j - 1) * 0.5,
                                ymax = length(genes) + (j - 1) * 0.5 + 0.5,
                                fill = unlist(sort.df[,sec.color[j]]))
        plot.df <- rbind(plot.df,anno.plot)
        
        anno_y <- data.frame(breaks = length(genes) + (j - 1) * 0.5 + 0.25,
                             labels = secondary.desc$name[secondary.desc$base == secondary$base[j]])
        anno_y_labels <- rbind(anno_y_labels,anno_y)
      }
    }
  }
  
  ## build new, more complex y-axis labels
  y_labels <- data.frame(breaks = (1:length(genes) - 0.5),
                         labels = rev(genes))
  primary_y_label <- data.frame(breaks = -0.25,
                                labels = primary.name)
  # secondary_y_label <- data.frame(breaks = mean(c(anno.plot$ymin, secondary.plot$ymax)),
  #                                 labels = secondary.name)
  y_labels <- rbind(y_labels, 
                    primary_y_label, 
                    # secondary_y_label,
                    anno_y_labels)
  
  hlines <- data.frame(yintercept = c(-0.5, max(plot.df$ymax)))
  
  sort.lab <- sort.df %>%
    group_by_(primary$id, primary$label) %>%
    summarise(xmean = mean(c(min(xpos) - 1, xpos)),
              y = length(genes) + 1,
              angle = 90) %>%
    ungroup() %>%
    select_(primary$id, primary$label,"xmean","y","angle")
  names(sort.lab)[1:2] <- c("primary_id","primary_label")
  
  if(showids) {
    sort.lab$primary_label <- paste(sort.lab$primary_id, sort.lab$primary_label)
  }
  
  if(rotatelabels) {
    sort.lab$primary_label <- gsub("[_|;| ]","\n",sort.lab$primary_label)
    sort.lab$angle <- 0
  }
  
  # Segments that divide groups
  segment_lines <- sort.df %>%
    group_by_(primary$id, primary$label) %>%
    summarise(x = max(xpos)) %>%
    mutate(xend = x,
           y = -0.5,
           yend = max(plot.df$ymax)) %>%
    select(x, xend, y , yend) %>%
    as.data.frame()
  
  ##############
  ## Plotting ##
  ##############
  
  p <- ggplot() + 
    # Main heatmap
    geom_rect(data = plot.df,
              aes(xmin = xmin, 
                  xmax = xmax, 
                  ymin = ymin, 
                  ymax = ymax, 
                  fill = fill)) +
    # Axis labels
    scale_x_continuous(breaks = sort.lab$xmean,
                       labels = sort.lab$primary_label,
                       expand = c(0, 0)) +
    scale_y_continuous(breaks = y_labels$breaks,
                       labels = y_labels$labels,
                       expand = c(0, 0)) +
    # fill and theme
    scale_fill_identity(guide=F) +
    theme_classic(base_size=pfontsize) +
    theme(axis.ticks = element_line(size=0.2),
          legend.title = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.line = element_blank())
  
  
  # cluster guide lines
  if(showlines) {
    
    p <- p + 
      geom_segment(data = segment_lines,
                   aes(x = x, 
                       xend = xend, 
                       y = y, 
                       yend = yend),
                   size = 0.2) +
      geom_vline(aes(xintercept = 0), 
                 size = 0.2) + 
      geom_hline(data = hlines,
                 aes(yintercept = yintercept), 
                 size = 0.2)
    
  }
  
  
  # Rotate labels checkbox
  if(rotatelabels) {
    p <- p + 
      theme(axis.text.x = element_text(angle = 0, 
                                       hjust = 0.5, 
                                       vjust=1))
  } else {
    p <- p + 
      theme(axis.text.x = element_text(angle = 90, 
                                       hjust = 1, 
                                       vjust = 0.5))
  }
  
  p
  
}


build_legend_plot <- function(data_source, 
                              genes = c("Hspa8","Snap25","Gad2","Vip"), 
                              autorange = "auto", 
                              minrange = 0, 
                              maxrange = 10,
                              scale_type = "scale.log",
                              pfontsize = 14) {
  
  library(dplyr)
  library(ggplot2)
  library(feather)
  
  data <- read_feather(file.path(data_source,"data.feather"), columns = genes)
  
  genes <- sub("-", ".", genes)
  genes[grepl("^[0-9]",genes)] <- paste0("X",genes[grepl("^[0-9]",genes)])
  
  names(data) <- sub("-",".",names(data))
  colors <- colorRampPalette(c("darkblue","white","red"))(1001)
  
  if(autorange == "auto") {
    min.val <- 0
    max.val <- max(unlist(data[, genes]))
  } else if (autorange == "manual") {
    min.val <- as.numeric(minrange)
    max.val <- as.numeric(maxrange)
  }
  
  ## Build geom_rect() compatible table
  legend_data <- data.frame(xmin = 1:1001,
                            xmax = 1:1001+1,
                            ymin = 0,
                            ymax = 1,
                            fill = colors)
  
  if(scale_type == "scale.abs") {
    scale_name <- "RPKM"
  } else if(scale_type == "scale.log") {
    scale_name <- "log10(RPKM + 1)"
    min.val <- log10(min.val + 1)
    max.val <- log10(max.val + 1)
  } else if(scale_type == "scale.rel") {
    scale_name <- "RPKM/max(RPKM)"
    min.val <- min.val/max.val
    max.val <- 1
  } else if(scale_type == "scale.log.rel") {
    scale_name <- "log10(RPKM + 1)/max(log10(RPKM + 1))"
    min.val <- log10(min.val + 1)
    max.val <- log10(max.val + 1)
    min.val <- min.val/max.val
    max.val <- 1
  }
  
  segment_data <- data.frame()
  
  legend_plot <- ggplot(legend_data) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill)) +
    geom_segment(aes(x = min(xmin), xend = max(xmax), y = 0, yend = 0)) +
    scale_fill_identity() +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(scale_name, breaks=c(0,250,500,750,1000),
                       labels=round(seq(min.val, max.val, by = (max.val-min.val)/4),2)) +
    theme_classic(base_size = pfontsize) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          axis.title.y = element_blank(),
          axis.line.x = element_blank())
  
  return(legend_plot)
}  




##########################################################################################################
### FUNCTIONS FOR MAKING DOT PLOTS


# Function to subset a dendrogram by retrieving a node with 
# a given attribute.
get_node_dend <- function(x, match_attr, match_value) {
  
  output <- NULL
  
  for(i in seq_len(length(x))) {
    
    if(attr(x[[i]], match_attr) == match_value) {
      
      output <- x[[i]]
      
    } else {
      if(is.list(x[[i]])) {
        nest <- get_node_dend(x[[i]], match_attr, match_value)
        if(!is.null(nest)) {
          output <- nest
        }
      }
    }
    
  }
  return(output)
  
}

# Function to subsample cells
subsampleCells <- function(clustersF, subSamp=25, seed=5){
  # Returns a vector of TRUE false for choosing a maximum of subsamp cells in each cluster 
  # clustersF = vector of cluster labels in factor format
  kpSamp = rep(FALSE,length(clustersF))
  for (cli in unique(as.character(clustersF))){
    set.seed(seed)
    seed   = seed+1
    kp     = which(clustersF==cli)
    kpSamp[kp[sample(1:length(kp),min(length(kp),subSamp))]] = TRUE
  }
  return(kpSamp)
}

# Function to build the jittered layer plot + dendrogram
build_layer_plot <- function(anno,
                             dend,
                             cocl,
                             cluster_ids,
                             seed_val = 42,
                             right_pad = 10,
							 fillColor = c("#ECE09C","#F7F2DA","#A7D7DF","#C1E5E7","#D4EDED","#B3E3E3"),
							 textSize = 2,
							 maxPerCluster = Inf,  # How many cells to include per cluster (default is all)
							 seed = 1) {
  
  # cluster annotations
  anno$cl <- anno$dendcluster_id <- anno$cluster_id
  #anno$layer_label <- as.numeric(substr(anno$layer_label,2,2))  # NOT HELPFUL
  
  cluster_anno <- anno %>%
    # dendcluster_id = dendrogram order
    select(cl, dendcluster_id, cluster_id, cluster_label, cluster_color) %>%
    unique()
  
  # layer annotations to retain for the plot
  keep_layers <- c("L1","L2","L3","L4","L5","L6")
  
  # padding between layers
  xpad <- 0.1
  ypad <- 0.05
  
  # ranges to use in the y-dimension for jittering
  layer_ranges <- data.frame(layer_label = rev(keep_layers),
                             ymin = (1:6 - 1) + ypad,
                             ymax = 1:6 - ypad)
  
  # filter for cells with the selected layer annotations
  filtered_anno <- anno %>%
    filter(cluster_id %in% cluster_ids) %>%
    filter(layer_label %in% keep_layers)
  
  # Subsample cells
  subSamp <- subsampleCells(filtered_anno$cluster_id,maxPerCluster,seed)
  filtered_anno <- filtered_anno[subSamp,]
  
  # ranges to use in the x-dimension for jittering
  cluster_ranges <- filtered_anno %>%
    select(cluster_id, cluster_color, cluster_label, dendcluster_id) %>%
    unique() %>%
    arrange(dendcluster_id) %>%
    mutate(xmin = 1:n() - 1 + xpad,
           xmax = 1:n()     - xpad,
           xmid = 1:n() - 0.5)
  
  set.seed(seed_val)
  
  # Assign random positions within the x and y limits
  # for each cell, as defined above.
  plot_data <- filtered_anno %>%
    select(sample_id,
           dendcluster_id, cluster_color, cluster_label,
           layer_id, layer_color, layer_label) %>%
    left_join(layer_ranges) %>%
    left_join(cluster_ranges) %>%
    group_by(dendcluster_id, layer_id) %>%
    mutate(x = runif(n(),xmin + xpad, xmax - xpad),
           y = runif(n(),ymin + ypad, ymax - ypad),
           fill_color = cluster_color)
  
  # Layer background rectangles
  layer_rects <- layer_ranges %>%
    mutate(xmin = min(cluster_ranges$xmin) - xpad, xmax = max(cluster_ranges$xmax) + xpad) %>%
    mutate(fill = fillColor)  # LAST COLOR ADDED
  
  # Cluster color highlights at bottom of the plot
  cluster_rects <- cluster_ranges %>%
    mutate(ymin = -ypad, ymax = ypad)
  
  # Filter the dendrogram for just the clusters that are present
  # in the plot
  prune_dend_labels <- labels(dend)[!labels(dend) %in% filtered_anno$cluster_label]
  filtered_dend <- dend %>%
    prune(prune_dend_labels)
  # Convert with ggdend to segment coordinates, and rescale the plot
  dend_seg <- as.ggdend(filtered_dend)$segments %>%
    mutate(y = (y/max(y))*3 + max(layer_rects$ymax) + ypad,
           yend = (yend/max(yend))*3 + max(layer_rects$ymax) + ypad,
           x = x - 0.5,
           xend = xend - 0.5)
  dend_seg$col = "black"
  dend_seg$lwd = 1
  dend_seg$lty = 1
  
  # padding rectangle to align with the violin plots
  pad_rect <- data.frame(ymin = min(layer_rects$ymin),
                         ymax = max(layer_rects$ymax),
                         xmin = max(layer_rects$xmax),
                         xmax = max(layer_rects$xmax)*(1 + right_pad/100))
  
  p <- ggplot() +
    # right side padding for alignment
    # with the violin plots
    geom_rect(data = pad_rect,
              aes(xmin = xmin, xmax = xmax,
                  ymin = ymin, ymax = ymax,
                  fill = "#FFFFFF",
                  color = "#FFFFFF")) +
    # dendrogram segments
    geom_segment(data = dend_seg,
                 aes(x = x, xend = xend,
                     y = y, yend = yend,
                     size = lwd,
                     color = col),
                 lineend = "square") +
    # layer background rectangles
    geom_rect(data = layer_rects,
              aes(xmin = xmin, xmax = xmax,
                  ymin = ymin, ymax = ymax,
                  fill = fill)) +
    # cluster label rectangles
    geom_rect(data = cluster_rects,
              aes(xmin = xmin, xmax = xmax,
                  ymin = ymin, ymax = ymax,
                  fill = cluster_color)) +
    geom_rect(data = cluster_rects,
              aes(xmin = xmin, xmax = xmax,
                  ymin = -2 - ypad, ymax = -2,
                  fill = cluster_color)) +
    # jittered cell points
    geom_point(data = plot_data,
               aes(x = x,
                   y = y,
                   color = cluster_color),
               size = 0.1) +
    # cluster label tag rectangles
    geom_rect(data = cluster_ranges,
              aes(xmin = xmid - 0.5 + xpad/2,
                  xmax = xmid + 0.5 - xpad/2,
                  ymax = 0 - ypad,
                  ymin = -2),
              fill = "#CAD7D7")+
    # cluster label text
    geom_text(data = cluster_ranges,
              aes(x = xmid,
                  y = -2 + ypad,
                  label = cluster_label),
              angle = 90,
              vjust = 0.3,
              hjust = 0,
              size = textSize) +
    # Plot settings
    scale_color_identity() +
    scale_size(range = c(0.5,1), guide = FALSE) +
    scale_fill_identity() +
    scale_y_continuous(limits = c(-2.1, 9)) +
    scale_x_continuous(expand = c(0,0)) +
    theme_void(base_size = 7) +
    theme(text = element_text(size = 8),
          legend.box.spacing = unit(0,"pt"))
  
  p
}
