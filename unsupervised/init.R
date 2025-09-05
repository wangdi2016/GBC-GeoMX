### Check path ###

#parent_dir <- dirname(getwd())
#parent_dir <- "/Users/wangdi/work/DWX/089_sDAS_NCI_Lawrence_Report_DW2"
#parent_dir <- "/Users/wangdi/work/DWX/geomx/2024-08-05/Rose-GeoMX/089_sDAS_NCI_Lawrence_Report_DW2"
parent_dir <- "/Users/wangdi/work/DWX/geomx/2024-09-23/089_sDAS_NCI_Lawrence_Report_DW5-ALL-AOIs"
knitr::opts_knit$set(root.dir = parent_dir)
knitr::opts_chunk$set(warning=FALSE, message=FALSE, echo=FALSE)

# Start fresh
rm(list=base::setdiff(ls(), "parent_dir"))

# List of packages for session
.packages = c("easyalluvial", "plyr", "dplyr", "tibble","gridExtra", "FactoMineR", "EnvStats", "ggrepel", "parcats", "pheatmap", "plotly", "scales", "lme4", "ggplot2", "ggthemes", "GLMMadaptive", "MCMCglmm", "GGally", "ggallin", "stringi", "parallel", "openxlsx", "lubridate", "parcats", "ggforce", "tidyr", "RColorBrewer","viridis", "DT","ComplexHeatmap","GeomxTools", "gridExtra","circlize")

# Install CRAN packages (if not already installed)
#.inst <- .packages %in% installed.packages()
#if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst],repos = "http://cran.us.r-project.org")

# Load packages into session 
null <- lapply(.packages, require, character.only=TRUE)
#library(sDAS)

# Set color pal

pal <- c("#3A6CA1", "#FFD861", "#CF4244", "#47BAB4", "#474747", "#EB739E", "#318026", "#A66293", "#F28E2B", "#8F6954")

# Initiate a date tag & directory naming
date_tag <- format(Sys.Date(), "%Y%m%d")
outdir <- paste0("output_089_NCI_", date_tag)
outdata_dir <- file.path(outdir, "data")
fig_dir <- file.path(outdir, "figures")
object_dir <- file.path(outdir, "R_objects")
de_dir <- file.path(outdir, "de_dir")

# Initialize output directory for static images, tables, and serialized R objects.
dir.create(outdir, recursive = TRUE)
dir.create(fig_dir, recursive = TRUE)
dir.create(de_dir, recursive = TRUE)
dir.create(object_dir, recursive = TRUE)
dir.create(outdata_dir, recursive = TRUE)

#User defined thresholds
min_nc <- 15      # min nucleic count??
min_area <- 1000  # min area??
sd_multiplier <- 1 # 1 sd ?
nf_threshold <- 5 # normalization factor threshold
top_dim <- 25 # number of top features (in terms of CV)
# Set Data Table parameters
dt_params = 
   list(dom = "lfBtip",
        buttons = list(list(extend = "copy"),
                       list(extend = "csv", filename = "ExampleDataSummary.csv"),
                       list(extend = "excel", filename = "ExampleDataSummary.xlsx")),
        autoWidth = FALSE,
        searching = TRUE,
        scrollX = TRUE,
        pagingType = "simple",
        scrollCollapse = TRUE,
        fixedColumns = list(leftColumns = 1))

#Some functions for the analysis
# left pad numbers to the maximum number of characters
pad_zeros <- function(vec){
  x <- as.character(vec)
  max_chars <- max(nchar(x))
  res <- as.character(sapply(x, function(y){
    return(paste0(c(rep("0", max_chars - nchar(y)), y), collapse=""))
  }))
  return(res)
}
# Function for computing SNR
get_snr <- function(exp, backgrounds, x=1, type="SBR", n_processors=4){
  # exp has features as columns
  # type needs to be:
  #      "SBR": signal (S) to background (B) ratio = signal / background: 
  #      "SSR": signal to standard deviation ratio (S-B)/sigma_b 
  # n_processors number of processors to use
  # backgrounds are the control molecules.
  # singals are the signal molecules.
  if(!type %in% c("SBR", "SSR")){
    stop("Needs SBR or SSR for type")
  }
  background_features <- exp[, backgrounds]
  features <- exp
  # We want to keep the backgrounds for comparison
  # to_keep <- setdiff(colnames(features), backgrounds)
  # features <- features[,to_keep]
  
  if(class(backgrounds)=="numeric"){
    background_geoMean <- background_features
    background_sd <- rep(0, length(background_features))
  } else {
    background_geoMean <- apply(background_features, 1, EnvStats::geoMean)
    background_sd <- apply(background_features, 1, sd)
    # Rare cases, if sd = 0, SSR will not work. Replace with the next lowest SD
    background_sd[background_sd==0] <- as.numeric(sort(background_sd[background_sd>0])[1])
  }
  cl <- makeCluster(n_processors)
  clusterExport(cl=cl, varlist=c("x", "type", "features", "type", "background_geoMean", 
                                 "background_sd"), envir=environment())
  
  inner_res <- parLapply(cl, 1:ncol(features), function(j){
    if(type=="SBR"){
      out <- as.numeric((features[,j] / background_geoMean))
    } else if(type=="SSR"){
      out <- as.numeric((signif(features[,j], 5) - signif(background_geoMean, 5))/(x*background_sd)) # sig. digits.
    } else {
      stop("something went wrong. error code 1.")
    }
    #snr[,j] <<- as.numeric((exp[,j] - sample_LOD) / sample_geo_means)
    out <- data.frame(out)
    colnames(out) <- colnames(features)[j]
    rownames(out) <- rownames(features)
    return(out)
  })
  
  stopCluster(cl)
  the_snr <- do.call(cbind, inner_res)
  
  if(all(colnames(the_snr) == colnames(features)) == FALSE){
    warning("Column names were not in the original order and will be adjusted.")
    the_snr <- the_snr %>% dplyr::select(eval(colnames(features)))
  }
  
  if(any(is.na(the_snr))){
    stop("Some genes not computed. Check function.")
  } else {
    return(data.matrix(the_snr))
  }
}
# making colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
split_time <- function(x){
  strsplit(x, split="-")[[1]][1]
}
plot_one <- function(an, ex, the_feature, the_factor, scale_log2=TRUE){
  # an <- dfs[[3]]
  # ex <- dfs[[2]]
  # the_feature <- "CD34"
  # the_factor <- "aoi_type"
  ex <- as.data.frame(ex[,the_feature])
  colnames(ex) <- the_feature
  ex <- ex %>% as_tibble() %>% add_column(Sample_ID=row.names(ex), .before=1)
  
  an <- as.data.frame(an) %>% dplyr::select(eval(the_factor))
  an <- an %>% as_tibble() %>% add_column(Sample_ID=row.names(an), .before=1)
  df <- base::merge(ex, an, by="Sample_ID") 
  
  if(scale_log2){
    df$the_feature_col <- log2(df[,the_feature])
  } else {
    df$the_feature_col <- df[,the_feature]
  }
  df$the_factor_col <- df[,the_factor]
  p <- ggplot(df, aes(x=the_factor_col, y=the_feature_col)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(width=0.2, height=0) + 
    theme_bw() + 
    xlab(paste0(the_factor))
  if(scale_log2){
    p <- p + 
      ylab(paste0("log2 ", the_feature))
  } else {
    p <- p + 
      ylab(paste0(the_feature))
  }
  return(p)
}
headn <- function(x, n=3){
  x[1:n, 1:n]
}
filter_samples <- function(DFS, samples_to_keep){
  if(nrow(DFS[[1]]) != nrow(DFS[[2]]) | nrow(DFS[[1]]) != nrow(DFS[[3]])){
    stop("Check dimensions. Error 1.")
  }
  before <- nrow(DFS[[1]])
  DFS[[1]] <- DFS[[1]][samples_to_keep,]
  DFS[[2]] <- DFS[[2]][samples_to_keep,]
  DFS[[3]] <- DFS[[3]][samples_to_keep,]
  after <- nrow(DFS[[1]])
  message(paste0("Removed ", before - after, " samples/AOIs."))
  if(nrow(DFS[[1]]) != nrow(DFS[[2]]) | nrow(DFS[[1]]) != nrow(DFS[[3]])){
    stop("Check dimensions. Error 2.")
  }
  return(DFS)
}


#' @name getColorPalette
#' @rdname getColorPalette
#' @title Set Color Palette
#' @description Parse and set colors to match annotations and be consistent throughout the report.
#' @param \code{input} a \code{data.frame} consisting of factors of interest with levels to be assigned colors
#' @param \code{start} what color in the palette to start with, default 1
#' @param \code{custom} Custom color palette to use instead of the predefined palette (a \code{vector} of mode \code{character})
#' @param \code{n} Number of colors needed
#' @param \code{method} Type of color palette to return: 'Main', 'Map', or 'Other'
#' @details  This function parses and sets colors to match annotations and be consistent throughout the report.
#' @return  A list of length equal to the number of columns in \code{input}. Names of the elements of the list are equal to the column names in \code{input}. Within each element of the list, a named character vector giving the discrete color value.
#'
#' @export
#' @examples
#' library(GeomxTools)
#' library(RColorBrewer)
#' color_pal <- getColorPalette(pData(sDAS::kidney_targets)[,c('class','segment','region')])
#' custom_pal <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D")
#' color_pal2 <- getColorPalette(pData(sDAS::kidney_targets)[,c('class','segment','region')], custom=custom_pal)
#' color_pal3 <- getColorPalette(pData(sDAS::kidney_targets)[,c('class','segment','region')], custom=brewer.pal(7,"Set1"))

getColorPalette <- function(input,              # Input must be data.frame
                            start = 1 ,         # Which color to start with
                            custom = NULL,      # A custom color palette, vector of type character
                            method = "Map"){    # Color palette type (this could change in the future)

## Check input
   # input should be a data.frame
   if(!inherits(input, "data.frame")){
      stop("input must be a data.frame")
   }

   # start must be a positive numeric integer.
   if(!is.null(start)) {
      if (!inherits(start, "numeric")) {
        stop("start must be numeric.")
      }
      if (start %% 1 != 0) {
        stop(paste0("The start given, ", start, " is not an integer."))
      }
      if (start < 1) {
        stop("start must be >= 1.")
      }
   }

   # custom must be a vector of mode character
   if(!is.null(custom)) {
      if(!inherits(custom, "character")){
        stop("custom color palette must be a vector of mode character")
      }
      if(length(custom) < 1){
        stop("custom must be at least length 1")
      }
      if(is.vector(custom) & is.list(custom)){
        stop("custom must be a vector of mode character, not a list")
      }
   }

   # method mustbe a character either: "Map", "Main", or "Other"; "Map" is default
   if(!method %in% c("Map","Main","Other")){
    stop("method must be one of: 'Main', 'Map', or 'Other'")
   }


# Make a palette data frame for reference at start of report
           l <- start
           input <- as.data.frame(input)
           color_pal <- list()
           input[] <- lapply(input, factor)
           for (i in 1:ncol(input)) {
            input[, i] <-
             factor(input[, i], levels = levels(factor(input[, i])), order = TRUE)
            pal <- report_pals(method = method,
                               n = length(levels(factor(input[, i]))),
                               start = l,
                               custom = custom)
            names(pal) <- levels(input[, i])
            l <- l + length(levels(input[, i]))
            color_pal[[length(color_pal) + 1]] <- pal
           }
           names(color_pal) <- names(input)

  return(color_pal)
}


#Helper function:
# Grab colors depending on which method: Main, Map, other
report_pals <- function(method = NULL,          # which palette type (Main, Map, or Other)
                        n = NULL,                   # how many colors to include in the palette
                        start = NULL,                # what color to start on
                        custom = NULL) {          # list of custom colors to be used as a palette
 if(!is.null(custom)) {
  pal <- custom
 } else {
  if(method == "Main") {
   # define palette for use in boxplots / scatter plots with annotations, n = 10 colors
   pal <- c("#3A6CA1", "#FFD861", "#CF4244", "#47BAB4", "#474747", "#EB739E", "#318026", "#A66293", "#F28E2B", "#8F6954")
   #defaults: blues  ,  yellows ,   reds   ,   teals  ,   grays  ,   pinks  ,   greens ,  purples ,  oranges ,  browns
  } else if(method == "Map") {
   # define palette to be passed to annotation maps (heatmap), n = 20 colors
   pal <- c("#3A6CA1", "#FFD861", "#D86769", "#AEE8E2", "#999999", "#FABFD2", "#318026", "#A66293", "#F28E2B", "#E3C0AC",
            "#A0CBE8", "#9E7E20", "#FFBFBD", "#2CABA3", "#474747", "#EB739E", "#A0E391", "#E8C1DE", "#FCB36A", "#B0846B")
  } else if(method == "Other") {
   # define palette to be passed to other types of factors that could be useful in the future: cell types, genes, etc, n = 20 colors
   pal <- c("#3A6CA1", "#FFD861", "#D86769", "#AEE8E2", "#999999", "#FABFD2", "#318026", "#A66293", "#F28E2B", "#E3C0AC",
            "#A0CBE8", "#9E7E20", "#FFBFBD", "#2CABA3", "#474747", "#EB739E", "#A0E391", "#E8C1DE", "#FCB36A", "#B0846B")
  }
 }

 ntot <- n + start - 1
 if(ntot > length(pal)) {
  pal <- rep(pal, ceiling(ntot / length(pal)))
 }
 return(pal[seq(start, ntot)])
}

