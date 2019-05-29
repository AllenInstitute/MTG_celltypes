#' Read data from a directory of feather files
#' 
get_feather_data <- function(feather_dir, genes, group_by, group_ids) {
  
  library(dplyr)
  library(feather)
  
  data_file <- paste0(feather_dir, "/data.feather")
  anno_file <- paste0(feather_dir, "/anno.feather")
  
  data <- feather(data_file)
  
  # Read annotations and convert factors
  anno <- read_feather(anno_file) %>%
    mutate_if(is.factor, as.character)
  
  # If an _id column was a factor, it's now a character. Convert to numeric for sorting.
  id_cols <- names(anno)[grepl("_id$", names(anno)) & names(anno) != "sample_id"]
  anno[id_cols] <- lapply(anno[id_cols], as.numeric)
  
  # Check the provided genes against the column names in data_file
  data_names <- names(data)
  
  if(sum(genes %in% data_names) != length(genes)) {
    # Report if names don't match after ignorning case
    not_found <- genes[!toupper(genes) %in% toupper(data_names)]
    
    warning(paste(paste0(not_found, collapse = ", "), "not found in feather data!"))
    
    # Update genes to use names as formatted in data
    genes <- data_names[toupper(data_names) %in% toupper(genes)]
  }
  
  # Find column indexes for sample_id and the matched genes
  # This seems to be faster than calling data[,c("sample_id",genes)] directly
  data_cols <- which(data_names %in% c("sample_id", genes))
  
  # Read the data from the data feather file into memory
  gene_data <- data[,data_cols]
  
  # Change - to . in column names and genes
  colnames(gene_data) <- gsub("-",".",colnames(gene_data))
  genes <- gsub("-",".",genes)
  
  # rename the _id, _label, and _color for the group_by values for use in plotting
  all_anno <- anno %>%
    rename_("plot_id" = paste0(group_by,"_id"),
            "plot_label" = paste0(group_by,"_label"),
            "plot_color" = paste0(group_by,"_color"))
  
  # use the group_ids to retain the order provided by the group_ids argument
  cluster_order <- data.frame(group_ids = group_ids) %>%
    mutate(cluster_x = 1:n())
  
  # Filter and order the rows
  data <- left_join(all_anno, gene_data, by = "sample_id") %>%
    filter(plot_id %in% group_ids) %>%
    left_join(cluster_order, by = c("plot_id" = "group_ids")) %>%
    arrange(cluster_x) %>%
    mutate(xpos = 1:n()) %>%
    select(-plot_id) %>%
    rename_("plot_id" = "cluster_x")
  
  return(data)
}

#' Convert integers to scientific notation labels
#' 
#' @param in_num a numeric vector
#' @param sig_figs a number indicating how many significant figures should be displayed.
#' @return a character vector with numeric values reformatted in 1.2E3 format
#' 
#' @examples
#' my_numbers <- c(100,15.359,32687,.000468)
#' 
#' sci_label(my_numbers)
#' 
#' sci_label(my_numbers,sig_figs=3)
sci_label <- function(in_num, sig_figs = 2, type = "plot") {
  labels <- character()
  for(i in 1:length(in_num)) {
    x <- in_num[i]
    if(x == 0) {
      first <- paste0("0", ".", paste0(rep("0", sig_figs - 1), collapse="") )
    } else if(log10(x) %% 1 == 0) {
      first <- substr(x, 1, 1)
      if(sig_figs > 1) {
        first <- paste0(first, ".", paste0(rep("0", sig_figs - 1), collapse=""))
      }
    } else {
      first <- round(x / (10 ^ floor(log10(x))), sig_figs - 1)
    }
    if(x == 0) {
      if(type == "plot") {
        label <- paste0(first, "%*%10^0" )
      } else if(type == "datatable") {
        label <- paste0(first, "\u271510<sup>0</sup>" )
      }
    } else {
      if(type == "plot") {
        label <- paste0(first, "%*%10^", floor(log10(x)))
      } else if(type == "datatable") {
        label <- paste0(first, "\u271510<sup>", floor(log10(x)),"</sup>")
      }
    }
    labels <- c(labels, label)
  }
  return(labels)
}

#' Convert font sizes in pt to mm
#' @param pt A numeric font size in pt.
#' @return A numeric font size in mm.
#' 
#' @examples
#' pt2mm(12)
#' 
#' ggplot(mtcars) +
#'   geom_text(aes(x = mpg, y = wt, label = rownames(mtcars)),
#'             size = pt2mm(7))
pt2mm <- function(pt) {
  mm <- pt / 2.834645669
  return(mm)
}

#' Build colorful, rectangular labels for plot headers in plot space
#' 
build_header_labels <- function(data, ngenes, nsamples, nclust, labelheight = 25, labeltype = "simple") {
  
  # Three label types: 
  # simple, which is for use with cluster-based plots
  # angle, for cell-based plots with "angle"-type polygonal labels
  # square, for cell-based plots with "square"-type labels
  
  ## Note on plot dimensions
  # The range of the plot area (not including labels) will be
  # y-axis: 1:(ngenes + 1)
  # x-axis: 0:(nsamples) (for cell-based plots)
  # x-axis: 1:(nclust + 1) (for cluster-based plots)
  
  labheight <- ngenes*(labelheight/100)/(1-labelheight/100)
  
  data <- data %>%
    group_by(plot_id,plot_label,plot_color) %>%
    summarise(minx = min(xpos),
              maxx = max(xpos))
  
  if(labeltype == "simple") {
    xlab.rect <- data.frame(xmin = 1:nclust - 0.5,
                            xmax = 1:nclust + 0.5,
                            ymin = ngenes + 1,
                            ymax = ngenes + 1 + labheight,
                            color = data$plot_color,
                            label = data$plot_label )
  }
  
  if(labeltype == "angle") {
    xlab.rect <- data.frame(xmin = (nsamples) * (1:nclust - 1) / nclust,
                            xmax = (nsamples) * (1:nclust) / nclust,
                            # 10% of the label height is reserved for angled polygons
                            ymin = ngenes + 1 + labheight*0.1,
                            ymax = ngenes + 1 + labheight,
                            color = data$plot_color,
                            label = data$plot_label )
  }
  if(labeltype == "square") {
    xlab.rect <- data %>% 
      group_by(plot_id) %>%
      summarise(xmin = minx - 1,
                xmax = maxx,
                ymin = ngenes + 1 + labheight * 0.1,
                ymax = ngenes + 1 + labheight,
                color = plot_color[1],
                label = plot_label[1])
  }
  
  xlab.rect  
}


#' Violin plots of gene expression for clusters
#' 
#' This function will generate plots similar to Figure 1c of Tasic, et al. (2015).
#' Warning: this is currently only able to work with internally-supplied datasets (v1_data and v1_anno).
#' Extension to user-supplied datasets will come soon.
#' 
#' @param genes A character vector containing gene symbols to be plotted
#' @param clusters A numeric vector containing clusters to plot (for v1_anno, the range is 1:49)
#' @param data_source A character object defining where the data is stored. Currently only works with "internal"
#' @param logscale Logical object, determines if data is log scaled before plotting.
#' @param fontsize numeric object, the font size (in pts) used to make the plot.
#' @param labelheight numeric object, Percent of the plot height that should be used for the labels (0 to 100).
#' 
#' @return a ggplot2 plot object
#' 
#' @examples
#' pottery_plot()
#' 
#' my_genes <- c("Ercc6","Ercc8","Trp53","Pgbd5")
#' my_clusters <- c(1,5,9,10,24,37)
#' pottery_plot(my_genes,my_clusters,logscale=T,fontsize=14)
group_violin_plot <- function(genes = c("Hspa8","Snap25","Gad2","Vip"),
                              group_by = "cluster", clusters = 1:10,
                              data_source,
                              sort = F, logscale = F, showcounts = T, rotatecounts = F,
                              fontsize = 7, labelheight = 25,
                              max_width = 10) {
  library(dplyr)
  library(feather)
  library(ggplot2)
  
  genes <- rev(genes)
  
  # get_feather_data() from data_formatting.R
  data <- get_feather_data(data_source,genes,group_by,clusters)
  
  genes <- genes[genes %in% names(data)]
  
  data <- data %>%
    select(-xpos) %>%
    mutate(xpos = plot_id)
  
  genes <- sub("-",".",genes)
  genes[grepl("^[0-9]",genes)] <- paste0("X",genes[grepl("^[0-9]",genes)])
  names(data)[grepl("^[0-9]",genes)] <- paste0("X",names(data)[grepl("^[0-9]",genes)])
  
  ngenes <- length(genes)
  nclust <- length(unique(data$plot_id))
  nsamples <- nrow(data)
  
  # Compute maximum values before scaling to plot space
  max_vals <- data %>% 
    select(one_of(genes)) %>% 
    summarise_all(max) %>% 
    unlist()
  
  # Variance injection
  # geom_violin() requires some variance, so I add a vanishingly small random number to each data value
  data[genes] <- data[genes] + runif(nrow(data),0,0.00001)
  
  # Scale the data between i and i + 0.9
  for(i in 1:length(genes)) {
    gene <- genes[i]
    gene_max <- max_vals[i]
    if(logscale) {
      data[gene] <- log10(data[gene] + 1) / log10(gene_max + 1) * 0.9 + i
    } else {
      data[gene] <- data[gene] / gene_max * 0.9 + i
    }
  }
  
  header_labels <- build_header_labels(data = data, ngenes = ngenes, nsamples = 1, nclust = nclust, labelheight = labelheight, labeltype = "simple")
  
  # Build the maximum value labels for the right edge
  max_labels <- data.frame(x = (nclust + 0.5) * 1.01,
                           y = 1:ngenes + 0.5,
                           label = sci_label(max_vals) )
  max_header <- data.frame(x = (nclust + 0.5) * 1.01,
                           y = ngenes + 1,
                           label = "Max value")
  max_width <- nclust*(max_width/100)/(1-max_width/100)
  
  label_y_size <- max(header_labels$ymax) - min(header_labels$ymin)
  
  cluster_data <- data %>%
    group_by(plot_label,plot_color,plot_id) %>%
    summarise(cn=n()) %>%
    as.data.frame(stringsAsFactors=F) %>%
    arrange(plot_id) %>%
    mutate(labely = ngenes + label_y_size*0.05,
           cny = max(header_labels$ymax) - 0.1*label_y_size,
           xpos = plot_id)
  
  # Plot setup
  p <- ggplot() +
    scale_fill_identity() +
    scale_y_continuous("", breaks = 1:length(genes) + 0.45, labels = genes, expand = c(0, 0)) +
    scale_x_continuous("", expand = c(0, 0)) +
    theme_classic(fontsize) +
    theme(axis.text = element_text(size = rel(1)),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none") +
    geom_hline(aes(yintercept = 1:(ngenes)), size = 0.2)
  
  # plot the violins for each gene
  for(i in 1:length(genes)) {
    p <- p + 
      geom_violin(data = data,
                  aes_string(x = "xpos", y = genes[i], fill = "plot_color"),
                  scale = "width", adjust = 2) +
      stat_summary(data = data,
                   aes_string(x = "xpos", y = genes[i]),
                   fun.y = "median", fun.ymin = "median", fun.ymax = "median", geom = "point", size = 0.7)
  }
  
  # Cluster labels
  p <- p +
    geom_rect(data = header_labels, 
              aes(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax, fill = color)) +
    geom_text(data = header_labels,
              aes(x = (xmin + xmax) / 2, y = ymin + 0.05, label = label),
              angle = 90, vjust = 0.35, hjust = 0, size = pt2mm(fontsize)) +
    # Maximum value labels on right side of plot
    geom_rect(aes(xmin = nclust + 0.5, xmax = nclust + 0.5 + max_width, ymin = 1, ymax = max(header_labels$ymax)),
              fill = "white") +
    geom_text(data = max_header,
              aes(x = x, y = y, label = label),
              angle = 90, hjust = 0, vjust = 0.35, size = pt2mm(fontsize)) +
    geom_text(data = max_labels,
              aes(x = x, y = y, label = label),
              hjust = 0, vjust = 0.35, size = pt2mm(fontsize), parse = TRUE)
  
  # Cluster counts
  if(showcounts) {
    if(rotatecounts) {
      p <- p + geom_text(data = cluster_data,
                         aes(y = cny, x = xpos, label = cn),
                         angle = 90,
                         vjust = 0.35,
                         hjust = 1,
                         size = pt2mm(fontsize))
    } else {
      p <- p + geom_text(data = cluster_data,
                         aes(y = cny, x = xpos,label = cn),
                         size = pt2mm(fontsize))
    }
  }
  
  p
}
