#' Automatically format an annotation file.  From libary(scrattch.io)
#'
#' This takes an anno file as input at any stage and properly annotates it for compatability with
#'   shiny and other scrattch functions.  In particular, it ensures that columns have a label,
#'   an id, and a color, and that there are no factors.  It won't overwrite columns that have
#'   already been properly process.
#'
#' @param anno an existing annotation data frame
#' @param scale_num should color scaling of numeric values be "predicted" (default and highly recommended;
#'   will return either "linear" or "log10" depending on scaling), "linear","log10","log2", or "zscore".
#' @param na_val_num The value to use to replace NAs for numeric columns. default = 0.
#' @param colorset_num A vector of colors to use for the color gradient.
#'   default = c("darkblue","white","red")
#' @param sort_label_cat a logical value to determine if the data in category columns
#'   should be arranged alphanumerically before ids are assigned. default = T.
#' @param na_val_cat The value to use to replace NAs in category and factor variables.
#'   default = "ZZ_Missing".
#' @param colorset_cat The colorset to use for assigning category and factor colors.
#'   Options are "varibow" (default), "rainbow","viridis","inferno","magma", and "terrain"
#' @param color_order_cat The order in which colors should be assigned for cat and
#'   factor variables. Options are "sort" and "random". "sort" (default) assigns colors
#'   in order; "random" will randomly assign colors.
#'
#' @return an updated data frame that has been automatically annotated properly
#'
#' @export
auto_annotate <- function (anno, scale_num = "predicted", na_val_num = 0, colorset_num = c("darkblue", 
    "white", "red"), sort_label_cat = TRUE, na_val_cat = "ZZ_Missing", 
    colorset_cat = "varibow", color_order_cat = "sort") 
{
    anno_out <- anno
    if (!is.element("sample_name", colnames(anno_out))) {
        colnames(anno_out) <- gsub("sample_id", "sample_name", 
            colnames(anno_out))
    }
    cn <- colnames(anno_out)
    convertColumns <- cn[(!grepl("_label", cn)) & (!grepl("_id", 
        cn)) & (!grepl("_color", cn))]
    convertColumns <- setdiff(convertColumns, "sample_name")
    for (cc in convertColumns) {
        value <- anno_out[, cc]
        if (is.numeric(value)) {
            if (is.element(scale_num, c("linear", "log10", "log2", 
                "zscore"))) {
                anno_out <- annotate_num(df = anno_out, col = cc, 
                  scale = scale_num, na_val = na_val_num, colorset = colorset_num)
            }
            else {
                scalePred <- ifelse(min(value) < 0, "linear", 
                  "log10")
                if ((max(value + 1)/min(value + 1)) < 100) {
                  scalePred <- "linear"
                }
                if (mean((value - min(value))/diff(range(value))) < 
                  0.01) {
                  scalePred <- "log10"
                }
                anno_out <- annotate_num(df = anno_out, col = cc, 
                  scale = scalePred, na_val = na_val_num, colorset = colorset_num)
            }
        }
        else {
            if (is.factor(value)) {
                anno_out <- annotate_factor(df = anno_out, col = cc, 
                  base = cc, na_val = na_val_cat, colorset = colorset_cat, 
                  color_order = color_order_cat)
            }
            else {
                anno_out <- annotate_cat(df = anno_out, col = cc, 
                  base = cc, na_val = na_val_cat, colorset = colorset_cat, 
                  color_order = color_order_cat, sort_label = sort_label_cat)
            }
        }
    }
    anno_out <- group_annotations(anno_out)
    anno_out
}


#' Generate colors and ids for numeric annotations.  From library(scrattch.io)
#'
#' @param df data frame to annotate
#' @param col name of the numeric column to annotate
#' @param base base name for the annotation, which wil be used in the desc table. If not provided, will use col as base.
#' @param scale The scale to use for assigning colors. Options are "linear","log10","log2, and "zscore"
#' @param na_val The value to use to replace NAs. default = 0.
#' @param colorset A vector of colors to use for the color gradient. default = c("darkblue","white","red")
#'
#' @return A modified data frame: the annotated column will be renamed base_label, and base_id and base_color columns will be appended
#'
#' @export
annotate_num <- function (df, col = NULL, base = NULL, scale = "log10", na_val = 0, 
    colorset = c("darkblue", "white", "red")) 
{
    if (class(try(is.character(col), silent = T)) == "try-error") {
        col <- lazyeval::expr_text(col)
    }
    else if (class(col) == "NULL") {
        stop("Specify a column (col) to annotate.")
    }
    if (class(try(is.character(base), silent = T)) == "try-error") {
        base <- lazyeval::expr_text(base)
    }
    else if (class(base) == "NULL") {
        base <- col
    }
    if (!is.numeric(df[[col]])) {
        df[[col]] <- as.numeric(df[[col]])
    }
    df[[col]][is.na(df[[col]])] <- na_val
    x <- df[[col]]
    annotations <- data.frame(label = unique(x)) %>% dplyr::arrange(label) %>% 
        dplyr::mutate(id = 1:dplyr::n())
    if (scale == "log10") {
        colors <- values_to_colors(log10(annotations$label + 
            1), colorset = colorset)
    }
    else if (scale == "log2") {
        colors <- values_to_colors(log2(annotations$label + 1), 
            colorset = colorset)
    }
    else if (scale == "zscore") {
        colors <- values_to_colors(scale(annotations$label), 
            colorset = colorset)
    }
    else if (scale == "linear") {
        colors <- values_to_colors(annotations$label, colorset = colorset)
    }
    annotations <- mutate(annotations, color = colors)
    names(annotations) <- paste0(base, c("_label", "_id", "_color"))
    names(df)[names(df) == col] <- paste0(base, "_label")
    df <- dplyr::left_join(df, annotations, by = paste0(base, 
        "_label"))
    df
}


#' Generate colors and ids for categorical annotations.  From library(scrattch.io)
#'
#' @param df data frame to annotate
#' @param col name of the character column to annotate
#' @param base base name for the annotation, which wil be used in the desc
#'   table. If not provided, will use col as base.
#' @param sort_label a logical value to determine if the data in col should be
#'   arranged alphanumerically before ids are assigned. default = T.
#' @param na_val The value to use to replace NAs. default = "ZZ_Missing".
#' @param colorset The colorset to use for assigning category colors. Options
#'   are "varibow" (default), "rainbow","viridis","inferno","magma", and "terrain"
#' @param color_order The order in which colors should be assigned. Options are
#'   "sort" and "random". "sort" assigns colors in order; "random" will randomly
#'   assign colors.
#'
#' @return A modified data frame: the annotated column will be renamed
#'   base_label, and base_id and base_color columns will be appended
#'
#' @export
annotate_cat <- function (df, col = NULL, base = NULL, sort_label = T, na_val = "ZZ_Missing", 
    colorset = "varibow", color_order = "sort") 
{
    if (class(try(is.character(col), silent = T)) == "try-error") {
        col <- lazyeval::expr_text(col)
    }
    else if (class(col) == "NULL") {
        stop("Specify a column (col) to annotate.")
    }
    if (class(try(is.character(base), silent = T)) == "try-error") {
        base <- lazyeval::expr_text(base)
    }
    else if (class(base) == "NULL") {
        base <- col
    }
    if (!is.character(df[[col]])) {
        df[[col]] <- as.character(df[[col]])
    }
    df[[col]][is.na(df[[col]])] <- na_val
    x <- df[[col]]
    annotations <- data.frame(label = unique(x), stringsAsFactors = F)
    if (sort_label) {
        annotations <- annotations %>% dplyr::arrange(label)
    }
    annotations <- annotations %>% dplyr::mutate(id = 1:n())
    if (colorset == "varibow") {
        colors <- varibow(nrow(annotations))
    }
    else if (colorset == "rainbow") {
        colors <- sub("FF$", "", grDevices::rainbow(nrow(annotations)))
    }
    else if (colorset == "viridis") {
        colors <- sub("FF$", "", viridisLite::viridis(nrow(annotations)))
    }
    else if (colorset == "magma") {
        colors <- sub("FF$", "", viridisLite::magma(nrow(annotations)))
    }
    else if (colorset == "inferno") {
        colors <- sub("FF$", "", viridisLite::inferno(nrow(annotations)))
    }
    else if (colorset == "plasma") {
        colors <- sub("FF$", "", viridisLite::plasma(nrow(annotations)))
    }
    else if (colorset == "terrain") {
        colors <- sub("FF$", "", grDevices::terrain.colors(nrow(annotations)))
    }
    else if (is.character(colorset)) {
        colors <- (grDevices::colorRampPalette(colorset))(nrow(annotations))
    }
    if (color_order == "random") {
        colors <- sample(colors, length(colors))
    }
    annotations <- dplyr::mutate(annotations, color = colors)
    names(annotations) <- paste0(base, c("_label", "_id", "_color"))
    names(df)[names(df) == col] <- paste0(base, "_label")
    df <- dplyr::left_join(df, annotations, by = paste0(base, 
        "_label"))
    df
}


#' Generate colors and ids for categorical annotations that are factors. From libary(scrattch.io)
#'
#' @param df data frame to annotate
#' @param col name of the factor column to annotate
#' @param base base name for the annotation, which wil be used in the desc
#'   table. If not provided, will use col as base.
#' @param na_val The value to use to replace NAs. default = "ZZ_Missing".
#' @param colorset The colorset to use for assigning category colors. Options
#'   are "varibow" (default), "rainbow","viridis","inferno","magma", and "terrain"
#' @param color_order The order in which colors should be assigned. Options are
#'   "sort" and "random". "sort" assigns colors in order; "random" will randomly
#'   assign colors.
#'
#' @return A modified data frame: the annotated column will be renamed
#'   base_label, and base_id and base_color columns will be appended
#'
#' @export
annotate_factor <- function (df, col = NULL, base = NULL, na_val = "ZZ_Missing", 
    colorset = "varibow", color_order = "sort") 
{
    if (class(try(is.character(col), silent = T)) == "try-error") {
        col <- lazyeval::expr_text(col)
    }
    else if (class(col) == "NULL") {
        stop("Specify a column (col) to annotate.")
    }
    if (class(try(is.character(base), silent = T)) == "try-error") {
        base <- lazyeval::expr_text(base)
    }
    else if (class(base) == "NULL") {
        base <- col
    }
    if (!is.factor(df[[col]])) {
        df[[col]] <- as.factor(df[[col]])
    }
    if (sum(is.na(df[[col]])) > 0) {
        lev <- c(levels(df[[col]]), na_val)
        levels(df[[col]]) <- lev
        df[[col]][is.na(df[[col]])] <- na_val
    }
    x <- df[[col]]
    annotations <- data.frame(label = as.character(levels(x)), 
        stringsAsFactors = F)
    annotations <- annotations %>% dplyr::mutate(id = 1:n())
    if (colorset == "varibow") {
        colors <- varibow(nrow(annotations))
    }
    else if (colorset == "rainbow") {
        colors <- sub("FF$", "", grDevices::rainbow(nrow(annotations)))
    }
    else if (colorset == "viridis") {
        colors <- sub("FF$", "", viridisLite::viridis(nrow(annotations)))
    }
    else if (colorset == "magma") {
        colors <- sub("FF$", "", viridisLite::magma(nrow(annotations)))
    }
    else if (colorset == "inferno") {
        colors <- sub("FF$", "", viridisLite::inferno(nrow(annotations)))
    }
    else if (colorset == "plasma") {
        colors <- sub("FF$", "", viridisLite::plasma(nrow(annotations)))
    }
    else if (colorset == "terrain") {
        colors <- sub("FF$", "", grDevices::terrain.colors(nrow(annotations)))
    }
    else if (is.character(colorset)) {
        colors <- (grDevices::colorRampPalette(colorset))(nrow(annotations))
    }
    if (color_order == "random") {
        colors <- sample(colors, length(colors))
    }
    annotations <- dplyr::mutate(annotations, color = colors)
    names(annotations) <- paste0(base, c("_label", "_id", "_color"))
    names(df)[names(df) == col] <- paste0(base, "_label")
    df[[col]] <- as.character(df[[col]])
    df <- dplyr::left_join(df, annotations, by = paste0(base, 
        "_label"))
    df
}


#' Generate a rainbow palette with variation in saturation and value.  From library(scrattch.io)
#'
#' @param n_colors The number of colors to generate
#'
#' @export
varibow <- function (n_colors) 
{
    sats <- rep_len(c(0.55, 0.7, 0.85, 1), length.out = n_colors)
    vals <- rep_len(c(1, 0.8, 0.6), length.out = n_colors)
    sub("FF$", "", grDevices::rainbow(n_colors, s = sats, v = vals))
}


#' Convert values to colors along a color ramp.  From library(scrattch.io)
#'
#' @param x a numeric vector to be converted to colors
#' @param min_val a number that's used to set the low end of the color scale (default = 0)
#' @param max_val a number that's used to set the high end of the color scale. If NULL (default),
#' use the highest value in x
#' @param colorset a set of colors to interpolate between using colorRampPalette
#' (default = c("darkblue","dodgerblue","gray80","orangered","red"))
#' @param missing_color a color to use for missing (NA) values.
#' 
#' @return a character vector of hex color values generated by colorRampPalette. Color values will
#' remain in the same order as x.
#' 
#' @export
values_to_colors <- function (x, min_val = NULL, max_val = NULL, colorset = c("darkblue", 
    "dodgerblue", "gray80", "orange", "orangered"), missing_color = "black") 
{
    heat_colors <- (grDevices::colorRampPalette(colorset))(1001)
    if (is.null(max_val)) {
        max_val <- max(x, na.rm = T)
    }
    else {
        x[x > max_val] <- max_val
    }
    if (is.null(min_val)) {
        min_val <- min(x, na.rm = T)
    }
    else {
        x[x < min_val] <- min_val
    }
    if (sum(x == min_val, na.rm = TRUE) == length(x)) {
        colors <- rep(heat_colors[1], length(x))
    }
    else {
        if (length(x) > 1) {
            if (var(x, na.rm = TRUE) == 0) {
                colors <- rep(heat_colors[500], length(x))
            }
            else {
                heat_positions <- unlist(round((x - min_val)/(max_val - 
                  min_val) * 1000 + 1, 0))
                colors <- heat_colors[heat_positions]
            }
        }
        else {
            colors <- heat_colors[500]
        }
    }
    if (!is.null(missing_color)) {
        colors[is.na(colors)] <- grDevices::rgb(t(grDevices::col2rgb(missing_color)/255))
    }
    colors
}


#' Group annotation columns.  From library(scrattch.io)
#'
#' @param df the annotation dataframe to arrange
#' @param sample_col the column with unique sample ids. Default is "sample_name".
#' @param keep_order a logical value. If FALSE, will sort the annotations alphanumerically by base.
#'
#' @return an annotation data frame with reordered columns
#'
#' @export
group_annotations <- function (df, sample_col = "sample_name", keep_order = TRUE) 
{
    labels <- names(df)[grepl("_label", names(df))]
    if (!keep_order) {
        labels <- labels[order(labels)]
    }
    bases <- sub("_label", "", labels)
    anno_cols <- c(paste0(rep(bases, each = 3), c("_id", "_label", 
        "_color")))
    extras <- setdiff(names(df), anno_cols)
    anno <- select(df, one_of(c(sample_col, anno_cols, extras)))
}


