###
source("./init.dw.R")

### load image from step1
load(file=file.path(object_dir,"step2.qc.norm.RData"))

#####################################
# 3 Exploratory Data Analysis (EDA) #
#####################################

# The following section is used to understand the general structure of the data and whether the hypotheses that were defined initially are likely to identify significant findings. This qualitative exploration can also help identify any study design considerations for downstream analyses.

compute_cv <- function(x){
  return(sd(x) / mean(x))
}
make_pca <- function(DFS, use_top_cv=TRUE, top_dim = NULL, use_zscores=TRUE, use_log2=FALSE){
  if(use_log2){
    dfs2 <- log2(DFS[[2]])
  } else {
    dfs2 <- DFS[[2]]
  }
  if(is.null(top_dim)){
    top_dim <- ncol(DFS[[2]])
  }
  # dfs2 <- t(dfs2)

  if(use_top_cv){
    cv <- apply(dfs2, 2, compute_cv)
    top_cv_features <- names(sort(cv, decreasing = TRUE)[1:top_dim])
    dfs2 <- dfs2[,eval(top_cv_features)]
  }
  if(use_zscores){
    for(i in 1:nrow(dfs2)){
      dfs2[i,] <- as.numeric(scale(as.numeric(dfs2[i,])), center = TRUE, scale=TRUE)
    }
  }
  pca.res <- PCA(dfs2, ncp = 3, graph = FALSE, scale.unit = TRUE)
  coords <- as.data.frame(pca.res$ind$coord[,1:3])
  pvar <- as.data.frame(pca.res$eig)[1:3,2]
  coords$Sample_ID <- row.names(coords)

  dfs3 <- as.data.frame(DFS[[3]])
  dfs3$Sample_ID <- row.names(dfs3)
  row.names(dfs3) <- NULL
  dfs3_pca <- base::merge(dfs3, coords, by="Sample_ID")
  return(list(dfs3_pca, pvar))
}

# 3.1 PCA
# In this section, we use principal component analysis, PCA, to identify broad patterns in the expression data.

# Analyst notes: The data separate primarily by PanCK status.

###################
# 3.1.0. make PCA #
###################
pca1 <- make_pca(DFS=dfs, 
                 use_top_cv = FALSE, 
                 use_zscores = FALSE, use_log2 = TRUE)
pvar <- pca1[[2]]
pca1 <- pca1[[1]]

###################################################################
# 3.1.1. Here, we color the points by PanCK (negative vs. positive). #
###################################################################
p <- ggplot(pca1, aes(x=Dim.1, y=Dim.2, label=row.names(dfs[[2]]))) +
     geom_point(aes(colour=PanCK), size=3, alpha=0.6) +
     xlab(paste0("PCA1: (Variation Explained ", round(pvar[1], 1), "%)")) +
     ylab(paste0("PCA2: (Variation Explained ", round(pvar[2], 1), "%)")) +
     theme_bw() + 
     scale_color_manual(values=pal_main[["PanCK"]])
#ggsave(p, filename=file.path(fig_dir, "PCA_all_features_CK.svg"), width=8, height=6)
ggsave(p, filename=file.path(fig_dir, "PCA_all_features_PanCK.png"), width=8, height=6)
p

#################################################################################
# 3.1.2 Here, we color the points by TumorEdge (Non tumor edge vs. Tumor edge). #
#################################################################################
#pca1 <- make_pca(DFS=dfs,
#                 use_top_cv = FALSE,
#                 use_zscores = FALSE, use_log2 = TRUE)
#pvar <- pca1[[2]]
#pca1 <- pca1[[1]]
p <- ggplot(pca1, aes(x=Dim.1, y=Dim.2, label=row.names(dfs[[2]]))) +
     geom_point(aes(colour=CancerType), size=3, alpha=0.6) +
     xlab(paste0("PCA1: (Variation Explained ", round(pvar[1], 1), "%)")) +
     ylab(paste0("PCA2: (Variation Explained ", round(pvar[2], 1), "%)")) +
     theme_bw() + 
     scale_color_manual(values=pal_main[["CancerType"]])
#ggsave(p, filename=file.path(fig_dir, "PCA_all_features_TumorEdge.svg"), width=8, height=6)
ggsave(p, filename=file.path(fig_dir, "PCA_all_features_CancerType.png"), width=8, height=6)
p

############################################################################################
# 3.1.3 Here, we color the points by Morphology (Combined stroma vs. Predominantly tumor). #
############################################################################################
#pca1 <- make_pca(DFS=dfs,
#                 use_top_cv = FALSE,
#                 use_zscores = FALSE, use_log2 = TRUE)
#pvar <- pca1[[2]]
#pca1 <- pca1[[1]]
p <- ggplot(pca1, aes(x=Dim.1, y=Dim.2, label=row.names(dfs[[2]]))) +
     geom_point(aes(colour=EpithelialStatus), size=3, alpha=0.6) +
     xlab(paste0("PCA1: (Variation Explained ", round(pvar[1], 1), "%)")) +
     ylab(paste0("PCA2: (Variation Explained ", round(pvar[2], 1), "%)")) +
     theme_bw() + 
     scale_color_manual(values=pal_main[["EpithelialStatus"]])
#ggsave(p, filename=file.path(fig_dir, "PCA_all_features_Morphology.svg"), width=8, height=6)
ggsave(p, filename=file.path(fig_dir, "PCA_all_features_EpithelialStatus.png"), width=8, height=6)
p

###########################
# 3.2 Heatmap of all AOIs #
###########################
hm_ann <- dfs[[3]][,factors_of_interest]
heatmap_dat <- t((dfs[[2]]))

## print dimension
dim(heatmap_dat)

label_hm <- rownames(heatmap_dat)
# rescale each protein (Z-scores)
for(i in 1:nrow(heatmap_dat)){
  heatmap_dat[i,] <- as.numeric(scale(log2(as.numeric(heatmap_dat[i,])), center = TRUE, scale=TRUE))
}
heatmap <- 
 ComplexHeatmap::Heatmap(heatmap_dat,
         col = circlize::colorRamp2(c(seq(-3, 3, 0.05)), 
                          c(colorRampPalette(c("#0092b5", "white", "#a6ce39"))(121))),
         name = 'z-score', use_raster = TRUE,
         clustering_distance_rows = "pearson", clustering_method_rows = "average",
         clustering_distance_columns = "pearson", clustering_method_columns = "average",
         #cluster_columns = FALSE,
         row_split = 2, column_split = 2,
         border_gp = grid::gpar(col = "darkgray"),
         show_row_names = FALSE, show_column_names = FALSE,
         right_annotation =
          ComplexHeatmap::rowAnnotation(foo=anno_mark(at=match(label_hm, rownames(heatmap_dat)),
                                      labels=label_hm)),
         top_annotation = 
          ComplexHeatmap::HeatmapAnnotation(df= hm_ann,
                            col = pal_main))  #,
                            #gp = grid::gpar(col = "gray")))
svg(filename = file.path(fig_dir, "Heatmap_norm_zscores_allFeatures.svg"), width=14, height=16)
h <- ComplexHeatmap::draw(heatmap, merge_legend = TRUE, heatmap_legend_side = "right", 
      annotation_legend_side = "right", adjust_annotation_extension = TRUE)
invisible(capture.output(dev.off())) 

# Save the heatmap as a PNG file
png(filename = file.path(fig_dir, "Heatmap_norm_zscores_allFeatures.png"), width = 14, height = 16, units = "in", res = 300)
draw(h)  # Draw the saved heatmap object h
invisible(dev.off())  # Close the PNG device

## save image
save.image(file=file.path(object_dir,"step3.EDA.RData"))

message("== step3 complete ==")

##
sessionInfo()
### End of Exploratory Data Analysis ###
