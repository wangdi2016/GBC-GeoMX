###
source("./init.R")

### load image from step1
load(file=file.path(object_dir,"step3.EDA.RData"))

# ```
# # Differential Expression Analyses
# 
# Differential expression (DE) is a common method used to compare protein expression between sample groups and to quantify the significance of these differences. We perform DE on a per-protein basis where the normalized expression is modeled using a linear mixed-effect model (LMM) to account for the sampling of multiple ROI/AOI segments per patient/tissue.
# 
# For a given comparison, we will subset the data down to the relevant segments and perform DE. We provide tables, violin plots, heatmaps, and volcano plots of the top differentially expressed proteins.
# 
# 
# This analysis report focuses on these key experimental comparisons (all comparisons are in the CK- compartment unless otherwise noted): 
# 
# **Question 1:** PanCK+ vs PanCK- AOIs
# 
# **Question 3:** Tumor edge vs non-tumor edge
# 
# **Question 5:** Proximal vs distal PanCK 
# 
# **Question 6:** Proximal vs distal CD45
# 
# **Question 7:** Luminal B vs Luminal A
# 
# **Question 12:** Younger vs older women
# 
# ```{r compute_lmer_func}
#Here are a few functions for computing DE and for plotting the results.
compute_lmer <- function(d, 
                         the_formula='protein_i ~ 1 + fixed + (1|random) ', 
                         the_factor="segment", the_random, the_baseline){
  # d is a list with d[[2]] the norm data and d[[3]] the annotations.
  to_return <- do.call(rbind, lapply(1:ncol(d[[2]]), function(i){
    protein_i <- colnames(d[[2]])[i]
    dat <- d[[3]]
    dat$protein_i <- log2(d[[2]][,i])
    dat$fixed <- factor(d[[3]][,the_factor], order=FALSE)
    dat$fixed <- relevel(dat$fixed, eval(the_baseline))
    if(length(levels(dat$fixed))!=2) stop("two levels are required.")
    
    dat$random <- factor(d[[3]][,the_random], order=FALSE)
    mod <- lmerTest::lmer(as.formula(the_formula), data=dat)
    coeffs <- as.data.frame(summary(mod)$coefficients)
    coeffs$beta <- rownames(coeffs)
    rownames(coeffs) <- NULL  
    out <- cbind(data.frame("Protein"=rep(protein_i, nrow(coeffs))), coeffs)
    colnames(out) <- c("Protein", "Est", "SE", "df", "t", "pval", "beta")
    # add singularity flag
    out$isSingular <- ifelse(isSingular(mod), 1, 0)
    return(out)
    }))
  return(to_return)
}
plot_volcano <- function(de_results, left, right, main, the_file=NULL){
  # de_results are from compute_lmer (only the beta you need)
  # left is the 0 factor name
  # right is the 1 factor name
  # main is the title
  # the_file is not NULL if you want to save the plots under the_file name
  
  
  de_results$padj <- p.adjust(de_results$pval, method="BH")
  de_results$FDR <- ifelse(de_results$padj < 0.05, "FDR<0.05", "FDR>=0.05")
  
  max_val <- max(abs(de_results$Est))*1.1
  p <- ggplot() + 
    geom_vline(xintercept = c(-1, 1), colour="grey", lty="dashed") +
    geom_point(data=de_results, 
               aes(x=Est, y=-log10(pval), colour=FDR)) +
    scale_color_manual(values=c("FDR<0.05"="black", "FDR>=0.05"="grey")) +
    geom_text_repel(data=de_results, 
                    aes(x=Est, y=-log10(pval), label=Protein), 
                    size=3, colour="grey20", max.overlaps = 32) + 
    scale_x_continuous(limits=c(-max_val, max_val)) + 
    theme_bw() + 
    ylab("Significance: -log10(P)") + 
    xlab(paste0(left, " << log2(FC) >> ", right)) + 
    ggtitle(main)
  if(!is.null(the_file)){
    ggsave(p, filename=the_file, width=5, height=5)
  }
  return(list(p, de_results))
}
v_plot2 <- function(df, to_label, label_color="grey", 
                    to_lowlight, lowlight_color, fc_threshold=0.58,
                    WIDTH=9, HEIGHT=5, SCALE=1.4, LWD=1,
                    out_dir=de_dir,
                    the_name, 
                    left="left", 
                    right="right", 
                    sizes){
  fc_cutoff = fc_threshold # 2^0.5849621 = 0.5 => 50% change.
  pval_cutoff = 0.05
  fileType = "svg"
  require(plyr)
  require(dplyr)
  
  dir.create(out_dir, recursive = TRUE)
  
  dfx <- df
  dfx$padj <- p.adjust(dfx$pval, method="BH")
  # Now plot the volcano
  # #####################
  # Basic plot of type ty
  # #####################
  dfx$sig_level <- ifelse(dfx$padj < 0.05, "FDR < 0.05", ifelse(dfx$pval < 0.05, "P < 0.05", "NS"))
  dfx$sig_level <- factor(dfx$sig_level, 
                          levels=c("NS", "P < 0.05", "FDR < 0.05"), order=TRUE)
  
  #the_x_lab <- paste0(left, " ", sprintf('\u2190'), " log", sprintf('\u2082'), "(FC) ", sprintf('\u2192'), " ", right)
  the_x_lab <- substitute(a %<-% s %->% b, list(a = paste0(left, "  "), b = paste0("  ", right), s = bquote("log"[2]*"FC")))
  

  
  if(!is.null(sizes)){
    colnames(sizes)[1] <- "Protein"
    dfx <- base::merge(dfx, sizes, by="Protein")
    plt <- ggplot(data=dfx, aes(x=Est, y=-log10(pval), label=Protein)) + 
    # geom_point(aes(size=-log10(pvalue)), colour="grey", alpha=0.7)
      geom_point(aes(colour=sig_level, size=log2(Median)), alpha=0.5) + 
      scale_colour_manual(values=c("NS"="grey", "P < 0.05"="dodgerblue", "FDR < 0.05"="red")) + 
      guides(colour=guide_legend(title="Significance\nLevel"), size=guide_legend(title="log2 UQ SNR"))
  } else {
    plt <- ggplot(data=dfx, aes(x=Est, y=-log10(pval), label=Protein)) + 
    # geom_point(aes(size=-log10(pvalue)), colour="grey", alpha=0.7)
      geom_point(aes(colour=sig_level), alpha=0.5) + 
      scale_colour_manual(values=c("NS"="grey", "P < 0.05"="dodgerblue", "FDR < 0.10"="red")) + 
      guides(colour=guide_legend(title="Significance\nLevel"))
  }
  # Dress plot
  plt <- plt + 
    geom_vline(xintercept=c(-fc_cutoff,fc_cutoff), lty = "dashed", lwd = 0.5, color = "gray") + 
    geom_hline(yintercept = -log10(pval_cutoff), lty = "dashed", lwd = 0.5, color = "grey") + 
    theme_bw(base_size = 14) + 
    ylab("-log10(P-value)") + 
    xlab(the_x_lab) + 
    scale_x_continuous(limits=c(min(dfx$Est)*SCALE, max(dfx$Est)*SCALE)) + 
    ggtitle(label = paste0(the_name))
    
    
  ignore_list <- c()
  if(!is.null(to_label)){
    
    if(is.numeric(to_label[1])){
     dfx_favs <- filter(dfx, 
                           Protein %in% unique(c(
                             (dfx %>% arrange(pval) %>% head(to_label[1]))$Protein)))
      top_favs <- filter(dfx_favs, abs(Est)>=fc_threshold & pval < 0.05)$Protein
      dfx_favs <- filter(dfx_favs, Protein %in% top_favs)
      
    } else {
      dfx_favs <- filter(dfx, Protein %in% to_label)
    }
    
    fav_need_plotting <- TRUE
      pltfav <- ggrepel::geom_text_repel(
        data=dfx_favs, aes(x=Est, y=-log10(pval), label=Protein),
        color = label_color, box.padding = 0.6, point.padding = 0.2,
        min.segment.length = 0.25, fontface = "bold", max.overlaps = 40)
    
    ignore_list <- c(ignore_list, as.character(dfx_favs$Protein))
  }
  
  to_add <- c()
  # Plot the to_lowlight genes
  if(!is.null(to_lowlight)){
    if(any(to_lowlight %in% dfx$Protein)){
      dfx_to_lowlight <- filter(dfx, Protein %in% to_lowlight)
      plt <- plt + 
        geom_point(data=dfx_to_lowlight, 
                   aes(x=Est, y=-log10(pval)), 
                   pch=1, size=2.5, stroke=LWD, colour=lowlight_color)
      # to_add <- paste0(" (", nrow(dfx_to_lowlight), "/", length(to_lowlight), ")")
    }
  }
  # We want the favorite proteins drawn on top of any other annotations.
  if(fav_need_plotting==TRUE){
    plt <- plt + pltfav + ggtitle(paste0(the_name, to_add))
  }
  
  dfx$relativeFC <- 2^abs(dfx$Est)*ifelse(dfx$Est<0, -1, 1) # -1.15 => 15% change; +4 = 400 percent change
  # Send to disk.
  prefix <- gsub("/", "_", gsub(" ", "_", the_name))
  out_file_name <- paste0(out_dir, "/", prefix, '_volcano', '.', fileType)
  ggsave(plt, filename = out_file_name, width=WIDTH, height=HEIGHT)
  
  out_file_name2 <- paste0(out_dir, "/", prefix, '.tsv')
  write.table(dfx, file = out_file_name2, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
  
  return(list(plt, dfx))
  
} 

message("==========================================")
message("== Question  1:** PanCK+ vs PanCK- AOIs ==")
message("==========================================")

# ```
## Contrast 1: PanCK+ vs PanCK- AOIs

# We will now compare PanCK+ vs PanCK- AOIs. The following formula was used to model differences for a given protein:

# $$ protein \sim PanCK + (1|Scan\_ID) $$

# We adjust for the multiple sampling of ROI segments per patient with the $Scan\_ID$ variable.
# ```{r de_contrast1, message=FALSE, warning=FALSE, echo=TRUE}
write.table(dfs, file="express.protein.tsv",sep="\t")
dfs_contrast1 <- dfs
to_keep_contrast1 <- row.names(dfs[[3]] %>% filter(PanCK %in% c("negative", "positive")))
dfs_contrast1 <- filter_samples(DFS=dfs_contrast1, samples_to_keep = to_keep_contrast1)
dt1 <- ddply(dfs_contrast1[[3]], .(PanCK), summarize, 
       nPatients=length(unique(Scan_ID)), 
       nAOIs=length(PanCK))
knitr::kable(dt1, format="html")

de_contrast1 <- compute_lmer(d=dfs_contrast1, the_formula='protein_i ~ 1 + fixed + (1 |random)  ', 
                         the_factor="PanCK", 
                         the_random = "Scan_ID",
                         the_baseline="negative")

de_contrast1 <- de_contrast1 %>% filter(isSingular==0, beta!="(Intercept)")
#de_contrast1
de_contrast1$FDR <- p.adjust(de_contrast1$pval, method = "BH")

#Proteins to label in plots
label_de_contrast1 <- de_contrast1[de_contrast1$pval<0.05,]$Protein


#```

#> Analyst notes: 43 proteins are significantly enriched at an FDR threshold of 0.05. 43 proteins are significantly enriched at a nominal P-value < 0.05. All proteins where P < 0.05 are labeled in the following figures.

#::: {#group_de_contrast1 .tabgroup}
#::::: {.tab}
#<button class="tablinks active" onclick="unrolltab(event, 'de_volcano_contrast1', 'group_de_contrast1')">Volcano</button>
#<button class="tablinks" onclick="unrolltab(event, 'de_heatmap_contrast1', 'group_de_contrast1')">Heatmap</button>
#<button class="tablinks" onclick="unrolltab(event, 'de_violin_contrast1_FOI', 'group_de_contrast1')">Violins</button>
#::::::: {#de_volcano_contrast1 .tabcontent style="display: block"}
#```{r de_volcano_contrast1,  fig.width=8, fig.height=6}
de_results_contrast1 <- v_plot2(df= de_contrast1,
        to_label=label_de_contrast1, 
        label_color="black", 
        fc_threshold = 0.50,
        to_lowlight = c(),
        lowlight_color = "grey", #Change this color to highlight any proteins of interest 
        WIDTH=9, HEIGHT=8, SCALE=1.4, LWD=1, out_dir=de_dir,
        the_name="PanCK (PanCK+ vs. PanCK-)",
        left="Negative", right="Positive",
        sizes=medians %>% dplyr::select(protein, Median)
        )
#de_results_contrast1[[1]]

#```
#
#:::::::
#::::::: {#de_heatmap_contrast1 .tabcontent}

#```{r de_heatmap_contrast1, echo=TRUE, message=FALSE, warning=FALSE, collapse=TRUE, fig.width=10, fig.height=12, eval = TRUE}
heatmap_contrast1_features <- de_contrast1
heatmap_contrast1_features <- heatmap_contrast1_features$Protein
heatmap_norm <- as.data.frame(heatmap_dat)
heatmap_dat_contrast1 <- t(scale(t(heatmap_norm[label_de_contrast1, to_keep_contrast1])))

#dfs_contrast1[[3]][factors_of_interest]

contrast1_de_heatmap <- 
 Heatmap(heatmap_dat_contrast1,
         col = colorRamp2(c(seq(-3, 3, 0.05)), 
                          c(colorRampPalette(c("#0092b5", "white", "#a6ce39"))(121))),
         name = 'z-score', use_raster = TRUE,
         clustering_distance_rows = "pearson", clustering_method_rows = "average",
         clustering_distance_columns = "pearson", clustering_method_columns = "average",
         row_split = 2, column_split = 2,
         border_gp = gpar(col = "darkgray"),
         show_row_names = FALSE, show_column_names = FALSE,
         right_annotation = 
          rowAnnotation(foo = anno_mark(at = match(label_de_contrast1, rownames(heatmap_dat_contrast1)),
                                        labels = label_de_contrast1)),
         top_annotation = 
          HeatmapAnnotation(df=dfs_contrast1[[3]][factors_of_interest],
                            col = pal_main )) #,
                            #gp = gpar(col = "gray")))
svglite::svglite(filename = file.path(de_dir, "contrast1_de_heatmap.svg"),
                width=10, height=12)
draw(contrast1_de_heatmap, merge_legend = TRUE, heatmap_legend_side = "right",
     annotation_legend_side = "right", adjust_annotation_extension = TRUE)
dev.off()

#png(filename = file.path(de_dir, "contrast1_de_heatmap.png"), width=10, height=12)
#draw(contrast1_de_heatmap, merge_legend = TRUE, heatmap_legend_side = "right",
#     annotation_legend_side = "right", adjust_annotation_extension = TRUE)
#dev.off()

#draw(contrast1_de_heatmap, merge_legend = TRUE, heatmap_legend_side = "right", 
#      annotation_legend_side = "right", adjust_annotation_extension = TRUE)
#```
#
#:::::::
#::::::: {#de_violin_contrast1_FOI .tabcontent}
#<!-- Make a series of violin plots of the label_de1 features. -->
#```{r de_violin_contrast1_FOI, eval = TRUE, echo = TRUE, fig.width=18, fig.height=12}
contrast_factor <- "PanCK"
comparison <- "positive vs negative"
p_adjust_class <- "FDR"
p_adjust_method <- "FDR"
violin_df <- cbind(dfs_contrast1[[3]]%>% dplyr::select(eval(contrast_factor)),
                   dfs_contrast1[[2]])
violin_df <- violin_df %>% tidyr::pivot_longer(cols=-1, names_to = "Protein", values_to = "Expression")
colnames(violin_df)[1] <- "contrast_factor"
violin_df <- violin_df %>% filter(Protein %in% label_de_contrast1)


# Format summary stats
violin_p_df <- de_contrast1 %>% filter(isSingular==0, beta!="(Intercept)")
violin_p_df$Comparison <- comparison
violin_p_df <- violin_p_df %>% tidyr::separate(col = Comparison, into=c("group1", "group2"), sep=" vs ")
violin_p_df["FDR"] <- signif(violin_p_df["FDR"], 3)
violin_exp_max <- ddply(violin_df, .(Protein), summarize,
                        y.position=((max(Expression)+1)*1.1)) # +1 for safe log2
violin_p_df <- base::merge(violin_p_df, violin_exp_max, by="Protein")
p <- ggplot(violin_df, 
            aes(x=contrast_factor, y=Expression, fill=contrast_factor)) + 
  geom_violin() +
  geom_jitter(width=0.25, height=0, size = 0.8, alpha=0.3) + 
  scale_fill_manual(values = pal_main[[contrast_factor]]) +
  facet_wrap(~Protein, scales = "free_y") +
  labs(x = eval(contrast_factor), y = "Expression (normalized counts)") +
  scale_y_continuous(trans = "log2", expand = expansion(mult = 0.2)) +
  theme_bw(base_size = 14) + 
  guides(fill=guide_legend(title = eval(contrast_factor)))
p <- p + ggprism::add_pvalue(
  violin_p_df,
  label="{p_adjust_class} = {violin_p_df[[p_adjust_method]]}", label.size = 3.6,
  y.position = violin_p_df$y.position
) + theme(legend.position="bottom")
ggsave(p, filename = file.path(de_dir, "contrast1_violins.svg"), width=18, height=12)
p

#```
#
#
#:::::::
#:::::
#:::

#The searchable table below lists log~2~ fold change estimates and p-values for each protein.

#```{r dt results contrast1}
dt_params$buttons <-   list(list(extend = "copy"),
                       list(extend = "csv", filename = "Contrast1"),
                       list(extend = "excel", filename = "Contrast1"))
dt_contrast1 <- de_contrast1
DT::datatable(
  dt_contrast1[, c(1,2,6, 9)] %>% arrange(pval),
  extensions = c("Buttons", "Scroller", "FixedColumns"),
  options = dt_params,
  rownames = FALSE) %>% DT::formatRound(columns=3, digits=3) %>% DT::formatSignif(columns=2:4)

#```

message("==============================================")
message("== Question  3:** Tumor edge vs non-tumor edge")
message("==============================================")

#```
## Contrast 3: Tumor edge vs non-tumor edge
#We will now compare tumor edge vs non-tumor edge PanCK- AOIs. The following formula was used to model differences for a given protein:
#$$ protein \sim TumorEdge + (1|Scan\_ID) $$
#We adjust for the multiple sampling of ROI segments per patient with the $Scan\_ID$ variable.
#```{r de_contrast3, message=FALSE, warning=FALSE, echo=TRUE}

dfs_contrast3 <- dfs
to_keep_contrast3 <- row.names(dfs[[3]] %>% filter(TumorEdge %in% c("Non tumor edge","Tumor edge"),
                                                   PanCK=="negative"))
dfs_contrast3 <- filter_samples(DFS=dfs_contrast3, samples_to_keep = to_keep_contrast3)
dt3 <- ddply(dfs_contrast3[[3]], .(TumorEdge), summarize, 
       nPatients=length(unique(Scan_ID)), 
       nAOIs=length(TumorEdge))
knitr::kable(dt3, format="html")

de_contrast3 <- compute_lmer(d=dfs_contrast3, the_formula='protein_i ~ 1 + fixed + (1 |random)  ', 
                         the_factor="TumorEdge", 
                         the_random = "Scan_ID",
                         the_baseline="Non tumor edge")

de_contrast3 <- de_contrast3 %>% filter(isSingular==0, beta!="(Intercept)")
de_contrast3$FDR <- p.adjust(de_contrast3$pval, method = "BH")

#Which proteins to label in plots
label_de_contrast3 <- de_contrast3[de_contrast3$pval<0.05,]$Protein

#```
#> Analyst notes: 23 proteins are significantly enriched at an FDR threshold of 0.05. 27 proteins are significantly enriched at a nominal P-value < 0.05. All proteins where P < 0.05 are labeled in the following figures.
#::: {#group_de_contrast3 .tabgroup}
#::::: {.tab}
#<button class="tablinks active" onclick="unrolltab(event, 'de_volcano_contrast3', 'group_de_contrast3')">Volcano</button>
#<button class="tablinks" onclick="unrolltab(event, 'de_heatmap_contrast3', 'group_de_contrast3')">Heatmap</button>
#<button class="tablinks" onclick="unrolltab(event, 'de_violin_contrast3_FOI', 'group_de_contrast3')">Violins</button>
#::::::: {#de_volcano_contrast3 .tabcontent style="display: block"}
#```{r de_volcano_contrast3,  fig.width=8, fig.height=6}

de_results_contrast3 <- v_plot2(df= de_contrast3,
        to_label=label_de_contrast3, 
        label_color="black", 
        fc_threshold = 0.50,
        to_lowlight = c(),
        lowlight_color = "grey", #Change this color to highlight any proteins of interest 
        WIDTH=9, HEIGHT=8, SCALE=1.4, LWD=1, out_dir=de_dir,
        the_name="TumerEdge (Tumor edge vs. Non tumor edge)",
        left="Non tumor edge", right="Tumor edge",
        sizes=medians %>% dplyr::select(protein, Median)
        )
de_results_contrast3[[1]]

#```
#
#:::::::
#::::::: {#de_heatmap_contrast3 .tabcontent}
#
#```{r de_heatmap_contrast3,  echo=TRUE, message=FALSE, warning=FALSE, collapse=TRUE, fig.width=10, fig.height=10, eval = TRUE}
heatmap_contrast3_features <- de_contrast3
heatmap_contrast3_features <- heatmap_contrast3_features$Protein
heatmap_norm <- as.data.frame(heatmap_dat)
heatmap_dat_contrast3 <- t(scale(t(heatmap_norm[label_de_contrast3, to_keep_contrast3])))
contrast3_de_heatmap <- 
 Heatmap(heatmap_dat_contrast3,
         col = colorRamp2(c(seq(-3, 3, 0.05)), 
                          c(colorRampPalette(c("#0092b5", "white", "#a6ce39"))(121))),
         name = 'z-score', use_raster = TRUE,
         clustering_distance_rows = "pearson", clustering_method_rows = "average",
         clustering_distance_columns = "pearson", clustering_method_columns = "average",
         row_split = 2, column_split = 2,
         border_gp = gpar(col = "darkgray"),
         show_row_names = FALSE, show_column_names = FALSE,
         right_annotation = 
          rowAnnotation(foo = anno_mark(at = match(label_de_contrast3, rownames(heatmap_dat_contrast3)),
                                        labels = label_de_contrast3)),
         top_annotation = 
          HeatmapAnnotation(df=dfs_contrast3[[3]][factors_of_interest],
                            col = pal_main))
svglite::svglite(filename = file.path(de_dir, "contrast3_de_heatmap.svg"),
                width=10, height=10)
draw(contrast3_de_heatmap, merge_legend = TRUE, heatmap_legend_side = "right",
     annotation_legend_side = "right", adjust_annotation_extension = TRUE)
dev.off()

#png(filename = file.path(de_dir, "contrast3_de_heatmap.png"), width=10, height=10)
#draw(contrast3_de_heatmap, merge_legend = TRUE, heatmap_legend_side = "right", 
#      annotation_legend_side = "right", adjust_annotation_extension = TRUE)
#dev.off()

#```
#:::::::
#::::::: {#de_violin_contrast3_FOI .tabcontent}
#<!-- Make a series of violin plots of the label_de3 features. -->
#```{r de_violin_contrast3_FOI, eval = TRUE, echo = TRUE, fig.width=18, fig.height=12}
contrast_factor <- "TumorEdge"
comparison <- "Tumor edge vs Non tumor edge"
p_adjust_class <- "FDR"
p_adjust_method <- "FDR"
violin_df <- cbind(dfs_contrast3[[3]]%>% dplyr::select(eval(contrast_factor)),
                   dfs_contrast3[[2]])
violin_df <- violin_df %>% tidyr::pivot_longer(cols=-1, names_to = "Protein", values_to = "Expression")
colnames(violin_df)[1] <- "contrast_factor"
violin_df <- violin_df %>% filter(Protein %in% label_de_contrast3)


# Format summary stats
violin_p_df <- de_contrast3 %>% filter(isSingular==0, beta!="(Intercept)")
violin_p_df$Comparison <- comparison
violin_p_df <- violin_p_df %>% tidyr::separate(col = Comparison, into=c("group1", "group2"), sep=" vs ")
violin_p_df["FDR"] <- signif(violin_p_df["FDR"], 3)
violin_exp_max <- ddply(violin_df, .(Protein), summarize,
                        y.position=((max(Expression)+1)*1.1)) # +1 for safe log2
violin_p_df <- base::merge(violin_p_df, violin_exp_max, by="Protein")
p <- ggplot(violin_df, 
            aes(x=contrast_factor, y=Expression, fill=contrast_factor)) + 
  geom_violin() +
  geom_jitter(width=0.25, height=0, size = 0.8, alpha=0.5) + 
  scale_fill_manual(values = pal_main[[contrast_factor]]) +
  facet_wrap(~Protein, scales = "free_y") +
  labs(x = eval(contrast_factor), y = "Expression (normalized counts)") +
  scale_y_continuous(trans = "log2", expand = expansion(mult = 0.2)) +
  theme_bw(base_size = 14) + 
  guides(fill=guide_legend(title = eval(contrast_factor)))
p <- p + ggprism::add_pvalue(
  violin_p_df,
  label="{p_adjust_class} = {violin_p_df[[p_adjust_method]]}", label.size = 3.6,
  y.position = violin_p_df$y.position
) + theme(legend.position="bottom")
ggsave(p, filename = file.path(de_dir, "contrast3_violins.svg"), width=18, height=12)
p

#```
#:::::::
#:::::
#:::
#The searchable table below lists log~2~ fold change estimates and p-values for each protein.
#```{r dt results contrast3}
dt_params$buttons <-   list(list(extend = "copy"),
                       list(extend = "csv", filename = "Contrast3"),
                       list(extend = "excel", filename = "Contrast3"))
dt_contrast3 <- de_contrast3
DT::datatable(
  dt_contrast3[, c(1,2,6, 9)] %>% arrange(pval),
  extensions = c("Buttons", "Scroller", "FixedColumns"),
  options = dt_params,
  rownames = FALSE) %>% DT::formatRound(columns=3, digits=3) %>% DT::formatSignif(columns=2:4)

message("=================================================================")
message("== Question  5:** ROI_neighborhood_PanCK PanCK-high vs. PanCK-low")
message("=================================================================")

# #```
## Contrast 5: Proximal vs distal PanCK
#We will now compare proximal PanCK vs distal PanCK AOIs from the CK- compartment. The following formula was used to model differences for a given protein:
#$$ protein \sim PanCK\_proximity + (1|Scan\_ID) $$
#We adjust for the multiple sampling of ROI segments per patient with the $Scan\_ID$ variable.
#```{r de_contrast5, message=FALSE, warning=FALSE, echo=TRUE}
dfs_contrast5 <- dfs
to_keep_contrast5 <- row.names(dfs[[3]] %>% filter(ROI_neighborhood_PanCK %in% c("PanCK-high", "PanCK-low"),
                                                   PanCK=="negative"))
dfs_contrast5 <- filter_samples(DFS=dfs_contrast5, samples_to_keep = to_keep_contrast5)
message("== dim dfs_comtrast5 ==")
head(dfs_contrast5)
dim(dfs_contrast5)
dt5 <- ddply(dfs_contrast5[[3]], .(ROI_neighborhood_PanCK), summarize, 
       nPatients=length(unique(Scan_ID)), 
       nAOIs=length(ROI_neighborhood_PanCK))
knitr::kable(dt5, format="html")

de_contrast5 <- compute_lmer(d=dfs_contrast5, the_formula='protein_i ~ 1 + fixed + (1 |random)  ', 
                         the_factor="ROI_neighborhood_PanCK", 
                         the_random = "Scan_ID",
                         the_baseline="PanCK-low")

de_contrast5 <- de_contrast5 %>% filter(isSingular==0, beta!="(Intercept)")
de_contrast5$FDR <- p.adjust(de_contrast5$pval, method = "BH")

#Proteins to label in plots
label_de_contrast5 <- de_contrast5[de_contrast5$pval<0.05,]$Protein

#```
#> Analyst notes: 14 proteins are significantly enriched at an FDR threshold of 0.05. 18 proteins are significantly enriched at a nominal P-value < 0.05. All proteins where P < 0.05 are labeled in the following figures.
#::: {#group_de_contrast5 .tabgroup}
#::::: {.tab}
#<button class="tablinks active" onclick="unrolltab(event, 'de_volcano_contrast5', 'group_de_contrast5')">Volcano</button>
#<button class="tablinks" onclick="unrolltab(event, 'de_heatmap_contrast5', 'group_de_contrast5')">Heatmap</button>
#<button class="tablinks" onclick="unrolltab(event, 'de_violin_contrast5_FOI', 'group_de_contrast5')">Violins</button>
#::::::: {#de_volcano_contrast5 .tabcontent style="display: block"}
#```{r de_volcano_contrast5,  fig.width=10, fig.height=6}
de_results_contrast5 <- v_plot2(df= de_contrast5 ,
        to_label=label_de_contrast5, 
        label_color="black", 
        fc_threshold = 0.50,
        to_lowlight = c(),
        lowlight_color = "grey", #Change this color to highlight any proteins of interest 
        WIDTH=9, HEIGHT=8, SCALE=1.4, LWD=1, out_dir=de_dir,
        the_name="ROI_neighborhood_PanCK (PanCK-high vs. PanCK-low)",
        left="PanCK-low", right="PanCK-high",
        sizes=medians %>% dplyr::select(protein, Median)
        )
#de_results_contrast5[[1]]

#```
#:::::::
#::::::: {#de_heatmap_contrast5 .tabcontent}
#```{r de_heatmap_contrast5, echo=TRUE, message=FALSE, warning=FALSE, collapse=TRUE, fig.width=10, fig.height=10, eval = TRUE}
heatmap_contrast5_features <- de_contrast5
heatmap_contrast5_features <- heatmap_contrast5_features$Protein
heatmap_norm <- as.data.frame(heatmap_dat)
heatmap_dat_contrast5 <- t(scale(t(heatmap_norm[label_de_contrast5, to_keep_contrast5])))
contrast5_de_heatmap <- 
 Heatmap(heatmap_dat_contrast5,
         col = colorRamp2(c(seq(-3, 3, 0.05)), 
                          c(colorRampPalette(c("#0092b5", "white", "#a6ce39"))(121))),
         name = 'z-score', use_raster = TRUE,
         clustering_distance_rows = "pearson", clustering_method_rows = "average",
         clustering_distance_columns = "pearson", clustering_method_columns = "average",
         row_split = 2, column_split = 2,
         border_gp = gpar(col = "darkgray"),
         show_row_names = FALSE, show_column_names = FALSE,
         right_annotation = 
          rowAnnotation(foo = anno_mark(at = match(label_de_contrast5, rownames(heatmap_dat_contrast5)),
                                        labels = label_de_contrast5)),
         top_annotation = 
          HeatmapAnnotation(df=dfs_contrast5[[3]][factors_of_interest],
                            col = pal_main))
svglite::svglite(filename = file.path(de_dir, "contrast5_de_heatmap.svg"),
                width=10, height=10)
draw(contrast5_de_heatmap, merge_legend = TRUE, heatmap_legend_side = "right",
     annotation_legend_side = "right", adjust_annotation_extension = TRUE)
dev.off()

#png(filename = file.path(de_dir, "contrast5_de_heatmap.png"), width=10, height=10)
#draw(contrast5_de_heatmap, merge_legend = TRUE, heatmap_legend_side = "right", 
#      annotation_legend_side = "right", adjust_annotation_extension = TRUE)
#dev.off()
#```
#
#:::::::
#::::::: {#de_violin_contrast5_FOI .tabcontent}
#<!-- Make a series of violin plots of the label_de5 features. -->
#```{r de_violin_contrast5_FOI, eval = TRUE, echo = TRUE, fig.width=18, fig.height=12}
contrast_factor <- "ROI_neighborhood_PanCK"
comparison <- "PanCK-high vs PanCK-low"
p_adjust_class <- "FDR"
p_adjust_method <- "FDR"
violin_df <- cbind(dfs_contrast5[[3]]%>% dplyr::select(eval(contrast_factor)),
                   dfs_contrast5[[2]])
violin_df <- violin_df %>% tidyr::pivot_longer(cols=-1, names_to = "Protein", values_to = "Expression")
colnames(violin_df)[1] <- "contrast_factor"
violin_df <- violin_df %>% filter(Protein %in% label_de_contrast5)


# Format summary stats
violin_p_df <- de_contrast5 %>% filter(isSingular==0, beta!="(Intercept)")
violin_p_df$Comparison <- comparison
violin_p_df <- violin_p_df %>% tidyr::separate(col = Comparison, into=c("group1", "group2"), sep=" vs ")
violin_p_df["FDR"] <- signif(violin_p_df["FDR"], 3)
violin_exp_max <- ddply(violin_df, .(Protein), summarize,
                        y.position=((max(Expression)+1)*1.1)) # +1 for safe log2
violin_p_df <- base::merge(violin_p_df, violin_exp_max, by="Protein")
p <- ggplot(violin_df, 
            aes(x=contrast_factor, y=Expression, fill=contrast_factor)) + 
  geom_violin() +
  geom_jitter(width=0.25, height=0, size = 0.8, alpha=0.5) + 
  scale_fill_manual(values = pal_main[[contrast_factor]]) +
  facet_wrap(~Protein, scales = "free_y") +
  labs(x = eval(contrast_factor), y = "Expression (normalized counts)") +
  scale_y_continuous(trans = "log2", expand = expansion(mult = 0.2)) +
  theme_bw(base_size = 14) + 
  guides(fill=guide_legend(title = eval(contrast_factor)))
p <- p + ggprism::add_pvalue(
  violin_p_df,
  label="{p_adjust_class} = {violin_p_df[[p_adjust_method]]}", label.size = 3.6,
  y.position = violin_p_df$y.position
) + theme(legend.position="bottom")
ggsave(p, filename = file.path(de_dir, "contrast5_violins.svg"), width=18, height=12)
p

#```
#:::::::
#:::::
#:::
#The searchable table below lists log~2~ fold change estimates and p-values for each protein.
#
#```{r dt results contrast5}
dt_params$buttons <-   list(list(extend = "copy"),
                       list(extend = "csv", filename = "contrast5"),
                       list(extend = "excel", filename = "contrast5"))
dt_contrast5 <- de_contrast5
DT::datatable(
  dt_contrast5[, c(1,2,6, 9)] %>% arrange(pval),
  extensions = c("Buttons", "Scroller", "FixedColumns"),
  options = dt_params,
  rownames = FALSE) %>% DT::formatRound(columns=3, digits=3) %>% DT::formatSignif(columns=2:4)

message("==============================================================")
message("== Question  6:** ROI_neighborhood_CD45 CD45-high vs. CD45-low")
message("==============================================================")

#```
## Contrast 6: Proximal vs distal CD45
#We will now compare proximal CD45 vs distal CD45 AOIs in the CK- compartment. The following formula was used to model differences for a given protein:
#$$ protein \sim CD45\_proximity + (1|Scan\_ID) $$
#We adjust for the multiple sampling of ROI segments per patient with the $Scan\_ID$ variable.
#```{r de_contrast6, message=FALSE, warning=FALSE, echo=TRUE}
dfs_contrast6 <- dfs
to_keep_contrast6 <- row.names(dfs[[3]] %>% filter(ROI_neighborhood_CD45 %in% c("CD45-low","CD45-high"),
                                                   PanCK=="negative"))
dfs_contrast6 <- filter_samples(DFS=dfs_contrast6, samples_to_keep = to_keep_contrast6)
dt6 <- ddply(dfs_contrast6[[3]], .(ROI_neighborhood_CD45), summarize, 
       nPatients=length(unique(Scan_ID)), 
       nAOIs=length(ROI_neighborhood_CD45))
knitr::kable(dt6, format="html")

de_contrast6 <- compute_lmer(d=dfs_contrast6, the_formula='protein_i ~ 1 + fixed + (1 |random)  ', 
                         the_factor="ROI_neighborhood_CD45", 
                         the_random = "Scan_ID",
                         the_baseline="CD45-low")
de_contrast6 <- de_contrast6 %>% filter(isSingular==0, beta!="(Intercept)")
de_contrast6$FDR <- p.adjust(de_contrast6$pval, method = "BH")

#Proteins to label in plots
label_de_contrast6 <- de_contrast6[de_contrast6$pval<0.05,]$Protein
#```
#
#> Analyst notes: 18 proteins are significantly enriched at an FDR threshold of 0.05. 21 proteins are significantly enriched at a nominal P-value < 0.05. All proteins where P < 0.05 are labeled in the following figures.

#::: {#group_de_contrast6 .tabgroup}
#::::: {.tab}
#<button class="tablinks active" onclick="unrolltab(event, 'de_volcano_contrast6', 'group_de_contrast6')">Volcano</button>
#<button class="tablinks" onclick="unrolltab(event, 'de_heatmap_contrast6', 'group_de_contrast6')">Heatmap</button>
#<button class="tablinks" onclick="unrolltab(event, 'de_violin_contrast6_FOI', 'group_de_contrast6')">Violins</button>
#::::::: {#de_volcano_contrast6 .tabcontent style="display: block"}
#```{r de_volcano_contrast6,  fig.width=9, fig.height=6}
de_results_contrast6 <- v_plot2(df= de_contrast6,
        to_label=label_de_contrast6, 
        label_color="black", 
        fc_threshold = 0.50,
        to_lowlight = c(),
        lowlight_color = "grey", #Change this color to highlight any proteins of interest 
        WIDTH=9, HEIGHT=8, SCALE=1.4, LWD=1, out_dir=de_dir,
        #the_name="Contrast 6  Proximal vs distal CD45",
        the_name="ROI_neighborhood_CD45 (CD45-high vs. CD45-low)",
        left="CD45-low", right="CD45-high",
        sizes=medians %>% dplyr::select(protein, Median)
        )

#de_results_contrast6[[1]]

#```
#
#:::::::
#::::::: {#de_heatmap_contrast6 .tabcontent}

#```{r de_heatmap_contrast6, echo=TRUE, message=FALSE, warning=FALSE, collapse=TRUE, fig.width=10, fig.height=10, eval = TRUE}
heatmap_contrast6_features <- de_contrast6
heatmap_contrast6_features <- heatmap_contrast6_features$Protein
heatmap_norm <- as.data.frame(heatmap_dat)
heatmap_dat_contrast6 <- t(scale(t(heatmap_norm[label_de_contrast6, to_keep_contrast6])))
contrast6_de_heatmap <- 
 Heatmap(heatmap_dat_contrast6,
         col = colorRamp2(c(seq(-3, 3, 0.05)), 
                          c(colorRampPalette(c("#0092b5", "white", "#a6ce39"))(121))),
         name = 'z-score', use_raster = TRUE,
         clustering_distance_rows = "pearson", clustering_method_rows = "average",
         clustering_distance_columns = "pearson", clustering_method_columns = "average",
         row_split = 2, column_split = 2,
         border_gp = gpar(col = "darkgray"),
         show_row_names = FALSE, show_column_names = FALSE,
         right_annotation = 
          rowAnnotation(foo = anno_mark(at = match(label_de_contrast6, rownames(heatmap_dat_contrast6)),
                                        labels = label_de_contrast6)),
         top_annotation = 
          HeatmapAnnotation(df=dfs_contrast6[[3]][factors_of_interest],
                            col = pal_main))
svglite::svglite(filename = file.path(de_dir, "contrast6_de_heatmap.svg"),
                width=10, height=10)
draw(contrast6_de_heatmap, merge_legend = TRUE, heatmap_legend_side = "right",
     annotation_legend_side = "right", adjust_annotation_extension = TRUE)
dev.off()

#png(filename = file.path(de_dir, "contrast6_de_heatmap.png"), width=10, height=10)
#draw(contrast6_de_heatmap, merge_legend = TRUE, heatmap_legend_side = "right", 
#      annotation_legend_side = "right", adjust_annotation_extension = TRUE)
#dev.off()

#```
#:::::::
#::::::: {#de_violin_contrast6_FOI .tabcontent}
#<!-- Make a series of violin plots of the label_de6 features. -->
#```{r de_violin_contrast6_FOI, eval = TRUE, echo = TRUE, fig.width=18, fig.height=12}
contrast_factor <- "ROI_neighborhood_CD45"
comparison <- "CD45-low vs CD45-high"
p_adjust_class <- "FDR"
p_adjust_method <- "FDR"
violin_df <- cbind(dfs_contrast6[[3]]%>% dplyr::select(eval(contrast_factor)),
                   dfs_contrast6[[2]])
violin_df <- violin_df %>% tidyr::pivot_longer(cols=-1, names_to = "Protein", values_to = "Expression")
colnames(violin_df)[1] <- "contrast_factor"
violin_df <- violin_df %>% filter(Protein %in% label_de_contrast6)


# Format summary stats
violin_p_df <- de_contrast6 %>% filter(isSingular==0, beta!="(Intercept)")
violin_p_df$Comparison <- comparison
violin_p_df <- violin_p_df %>% tidyr::separate(col = Comparison, into=c("group1", "group2"), sep=" vs ")
violin_p_df["FDR"] <- signif(violin_p_df["FDR"], 3)
violin_exp_max <- ddply(violin_df, .(Protein), summarize,
                        y.position=((max(Expression)+1)*1.1)) # +1 for safe log2
violin_p_df <- base::merge(violin_p_df, violin_exp_max, by="Protein")
p <- ggplot(violin_df, 
            aes(x=contrast_factor, y=Expression, fill=contrast_factor)) + 
  geom_violin() +
  geom_jitter(width=0.25, height=0, size = 0.8, alpha=0.5) + 
  scale_fill_manual(values = pal_main[[contrast_factor]]) +
  facet_wrap(~Protein, scales = "free_y") +
  labs(x = eval(contrast_factor), y = "Expression (normalized counts)") +
  scale_y_continuous(trans = "log2", expand = expansion(mult = 0.2)) +
  theme_bw(base_size = 14) + 
  guides(fill=guide_legend(title = eval(contrast_factor)))
p <- p + ggprism::add_pvalue(
  violin_p_df,
  label="{p_adjust_class} = {violin_p_df[[p_adjust_method]]}", label.size = 3.6,
  y.position = violin_p_df$y.position
) + theme(legend.position="bottom")
ggsave(p, filename = file.path(de_dir, "contrast6_violins.svg"), width=18, height=12)
#ggsave(p, filename = file.path(de_dir, "contrast6_violins.png"), width=18, height=12)
p

#```
#:::::::
#:::::
#:::
#The searchable table below lists log~2~ fold change estimates and p-values for each protein.
#```{r dt results contrast6}
dt_params$buttons <-   list(list(extend = "copy"),
                       list(extend = "csv", filename = "contrast6"),
                       list(extend = "excel", filename = "contrast6"))
dt_contrast6 <- de_contrast6
DT::datatable(
  dt_contrast6[, c(1,2,6, 9)] %>% arrange(pval),
  extensions = c("Buttons", "Scroller", "FixedColumns"),
  options = dt_params,
  rownames = FALSE) %>% DT::formatRound(columns=3, digits=3) %>% DT::formatSignif(columns=2:4)

message("========================================")
message("== Question  7:** Luminal B vs luminal A")
message("========================================")

#```
## Contrast 7: Luminal vs non-luminal
#We will now compare luminal vs non-luminal CK- AOIs. The following formula was used to model differences for a given protein:
#$$ protein \sim Luminal + (1|Scan\_ID) $$
#We adjust for the multiple sampling of ROI segments per patient with the $Scan\_ID$ variable.
#```{r de_contrast7, message=FALSE, warning=FALSE, echo=TRUE}
dfs_contrast7 <- dfs
to_keep_contrast7 <- row.names(dfs[[3]] %>% filter(PAM50 %in% c("LumB","LumA"),
                                                   PanCK=="negative"))
message("== Luminal B vs luminal A dim ==")
dim(to_keep_contrast7)
dfs_contrast7 <- filter_samples(DFS=dfs_contrast7, samples_to_keep = to_keep_contrast7)
dt7 <- ddply(dfs_contrast7[[3]], .(PAM50), summarize, 
       nPatients=length(unique(Scan_ID)), 
       nAOIs=length(PAM50))
knitr::kable(dt7, format="html")

de_contrast7 <- compute_lmer(d=dfs_contrast7, the_formula='protein_i ~ 1 + fixed + (1 |random)  ', 
                         the_factor="PAM50", 
                         the_random = "Scan_ID",
                         the_baseline="LumA")
de_contrast7 <- de_contrast7 %>% filter(isSingular==0, beta!="(Intercept)")
de_contrast7$FDR <- p.adjust(de_contrast7$pval, method = "BH")

#Proteins to label in plots
label_de_contrast7 <- de_contrast7[de_contrast7$pval<0.05,]$Protein
#```
#
#> Analyst notes: 19 proteins are significantly enriched at an FDR threshold of 0.05. 23 proteins are significantly enriched at a nominal P-value < 0.05. All proteins where P < 0.05 are labeled in the following figures.

#::: {#group_de_contrast7 .tabgroup}
#::::: {.tab}
#<button class="tablinks active" onclick="unrolltab(event, 'de_volcano_contrast7', 'group_de_contrast7')">Volcano</button>
#<button class="tablinks" onclick="unrolltab(event, 'de_heatmap_contrast7', 'group_de_contrast7')">Heatmap</button>
#<button class="tablinks" onclick="unrolltab(event, 'de_violin_contrast7_FOI', 'group_de_contrast7')">Violins</button>
#::::::: {#de_volcano_contrast7 .tabcontent style="display: block"}
#```{r de_volcano_contrast7,  fig.width=8, fig.height=6}
de_results_contrast7 <- v_plot2(df= de_contrast7,
        to_label=label_de_contrast7, 
        label_color="black", 
        fc_threshold = 0.50,
        to_lowlight = c(),
        lowlight_color = "grey", #Change this color to highlight any proteins of interest 
        WIDTH=9, HEIGHT=8, SCALE=1.4, LWD=1, out_dir=de_dir,
        the_name="PAM50 (LumB vs. LumA)",
        left="LumA", right="LumB",
        sizes=medians %>% dplyr::select(protein, Median)
        )
#de_results_contrast7[[1]]
#```
#
#:::::::
#::::::: {#de_heatmap_contrast7 .tabcontent}
#```{r de_heatmap_contrast7,  echo=TRUE, message=FALSE, warning=FALSE, collapse=TRUE, fig.width=10, fig.height=10, eval = TRUE}
heatmap_contrast7_features <- de_contrast7
heatmap_contrast7_features <- heatmap_contrast7_features$Protein
heatmap_norm <- as.data.frame(heatmap_dat)
heatmap_dat_contrast7 <- t(scale(t(heatmap_norm[label_de_contrast7, to_keep_contrast7])))
contrast7_de_heatmap <- 
 Heatmap(heatmap_dat_contrast7,
         col = colorRamp2(c(seq(-3, 3, 0.05)), 
                          c(colorRampPalette(c("#0092b5", "white", "#a6ce39"))(121))),
         name = 'z-score', use_raster = TRUE,
         clustering_distance_rows = "pearson", clustering_method_rows = "average",
         clustering_distance_columns = "pearson", clustering_method_columns = "average",
         row_split = 2, column_split = 2,
         border_gp = gpar(col = "darkgray"),
         show_row_names = FALSE, show_column_names = FALSE,
         right_annotation = 
          rowAnnotation(foo = anno_mark(at = match(label_de_contrast7, rownames(heatmap_dat_contrast7)),
                                        labels = label_de_contrast7)),
         top_annotation = 
          HeatmapAnnotation(df=dfs_contrast7[[3]][factors_of_interest],
                            col = pal_main))
svglite::svglite(filename = file.path(de_dir, "contrast7_de_heatmap.svg"),
                width=10, height=10)
draw(contrast7_de_heatmap, merge_legend = TRUE, heatmap_legend_side = "right",
     annotation_legend_side = "right", adjust_annotation_extension = TRUE)
dev.off()

#png(filename = file.path(de_dir, "contrast7_de_heatmap.png"), width=10, height=10)
#draw(contrast7_de_heatmap, merge_legend = TRUE, heatmap_legend_side = "right", 
#      annotation_legend_side = "right", adjust_annotation_extension = TRUE)
#dev.off()

#```
#:::::::
#::::::: {#de_violin_contrast7_FOI .tabcontent}
#<!-- Make a series of violin plots of the label_de7 features. -->
#```{r de_violin_contrast7_FOI, eval = TRUE, echo = TRUE, fig.width=14, fig.height=10}
contrast_factor <- "PAM50"
comparison <- "LumB vs LumA"
p_adjust_class <- "FDR"
p_adjust_method <- "FDR"
violin_df <- cbind(dfs_contrast7[[3]]%>% dplyr::select(eval(contrast_factor)),
                   dfs_contrast7[[2]])
violin_df <- violin_df %>% tidyr::pivot_longer(cols=-1, names_to = "Protein", values_to = "Expression")
colnames(violin_df)[1] <- "contrast_factor"
violin_df <- violin_df %>% filter(Protein %in% label_de_contrast7)

# Format summary stats
violin_p_df <- de_contrast7 %>% filter(isSingular==0, beta!="(Intercept)")
violin_p_df$Comparison <- comparison
violin_p_df <- violin_p_df %>% tidyr::separate(col = Comparison, into=c("group1", "group2"), sep=" vs ")
violin_p_df["FDR"] <- signif(violin_p_df["FDR"], 3)
violin_exp_max <- ddply(violin_df, .(Protein), summarize,
                        y.position=((max(Expression)+1)*1.1)) # +1 for safe log2
violin_p_df <- base::merge(violin_p_df, violin_exp_max, by="Protein")
p <- ggplot(violin_df, 
            aes(x=contrast_factor, y=Expression, fill=contrast_factor)) + 
  geom_violin() +
  geom_jitter(width=0.25, height=0, size = 0.8, alpha=0.5) + 
  scale_fill_manual(values = pal_main[[contrast_factor]]) +
  facet_wrap(~Protein, scales = "free_y") +
  labs(x = eval(contrast_factor), y = "Expression (normalized counts)") +
  scale_y_continuous(trans = "log2", expand = expansion(mult = 0.2)) +
  theme_bw(base_size = 14) + 
  guides(fill=guide_legend(title = eval(contrast_factor)))
p <- p + ggprism::add_pvalue(
  violin_p_df,
  label="{p_adjust_class} = {violin_p_df[[p_adjust_method]]}", label.size = 3.6,
  y.position = violin_p_df$y.position
) + theme(legend.position="bottom")
ggsave(p, filename = file.path(de_dir, "contrast7_violins.svg"), width=14, height=10)
#ggsave(p, filename = file.path(de_dir, "contrast7_violins.png"), width=14, height=10)
p

#```
#:::::::
#:::::
#:::
#The searchable table below lists log~2~ fold change estimates and p-values for each protein.
#```{r dt results contrast7}
dt_params$buttons <-   list(list(extend = "copy"),
                       list(extend = "csv", filename = "contrast7"),
                       list(extend = "excel", filename = "contrast7"))
dt_contrast7 <- de_contrast7
DT::datatable(
  dt_contrast7[, c(1,2,6, 9)] %>% arrange(pval),
  extensions = c("Buttons", "Scroller", "FixedColumns"),
  options = dt_params,
  rownames = FALSE) %>% DT::formatRound(columns=3, digits=3) %>% DT::formatSignif(columns=2:4)
#```

message("========================================")
message("== Question  8:** Her2-enriched vs luminal A")
message("========================================")

#```
## Contrast 8: Her2-enriched vs lumA
#We will now compare Her2-enriched vs LumA PanCK- AOIs. The following formula was used to model differences for a given protein:
#$$ protein \sim PAM50 + (1|Scan\_ID) $$
#We adjust for the multiple sampling of ROI segments per patient with the $Scan\_ID$ variable.
#```{r de_contrast8, message=FALSE, warning=FALSE, echo=TRUE}
dfs_contrast8 <- dfs
to_keep_contrast8 <- row.names(dfs[[3]] %>% filter(PAM50 %in% c("Her2-enriched","LumA"),
                                                   PanCK=="negative"))
message("== Her2-enriched vs luminal A dim ==")
dim(to_keep_contrast8)
dfs_contrast8 <- filter_samples(DFS=dfs_contrast8, samples_to_keep = to_keep_contrast8)
dt8 <- ddply(dfs_contrast8[[3]], .(PAM50), summarize, 
       nPatients=length(unique(Scan_ID)), 
       nAOIs=length(PAM50))
knitr::kable(dt8, format="html")

de_contrast8 <- compute_lmer(d=dfs_contrast8, the_formula='protein_i ~ 1 + fixed + (1 |random)  ', 
                         the_factor="PAM50", 
                         the_random = "Scan_ID",
                         the_baseline="LumA")
de_contrast8 <- de_contrast8 %>% filter(isSingular==0, beta!="(Intercept)")
de_contrast8$FDR <- p.adjust(de_contrast8$pval, method = "BH")

#Proteins to label in plots
label_de_contrast8 <- de_contrast8[de_contrast8$pval<0.05,]$Protein
#```
#
#> Analyst notes: 19 proteins are significantly enriched at an FDR threshold of 0.05. 23 proteins are significantly enriched at a nominal P-value < 0.05. All proteins where P < 0.05 are labeled in the following figures.

#::: {#group_de_contrast8 .tabgroup}
#::::: {.tab}
#<button class="tablinks active" onclick="unrolltab(event, 'de_volcano_contrast8', 'group_de_contrast8')">Volcano</button>
#<button class="tablinks" onclick="unrolltab(event, 'de_heatmap_contrast8', 'group_de_contrast8')">Heatmap</button>
#<button class="tablinks" onclick="unrolltab(event, 'de_violin_contrast8_FOI', 'group_de_contrast8')">Violins</button>
#::::::: {#de_volcano_contrast8 .tabcontent style="display: block"}
#```{r de_volcano_contrast8,  fig.width=8, fig.height=6}
de_results_contrast8 <- v_plot2(df= de_contrast8,
        to_label=label_de_contrast8, 
        label_color="black", 
        fc_threshold = 0.50,
        to_lowlight = c(),
        lowlight_color = "grey", #Change this color to highlight any proteins of interest 
        WIDTH=9, HEIGHT=8, SCALE=1.4, LWD=1, out_dir=de_dir,
        the_name="PAM50 (Her2-enriched vs. LumA)",
        left="LumA", right="Her2-enriched",
        sizes=medians %>% dplyr::select(protein, Median)
        )
#de_results_contrast8[[1]]
#```
#
#:::::::
#::::::: {#de_heatmap_contrast8 .tabcontent}
#```{r de_heatmap_contrast8,  echo=TRUE, message=FALSE, warning=FALSE, collapse=TRUE, fig.width=10, fig.height=10, eval = TRUE}
heatmap_contrast8_features <- de_contrast8
heatmap_contrast8_features <- heatmap_contrast8_features$Protein
heatmap_norm <- as.data.frame(heatmap_dat)
heatmap_dat_contrast8 <- t(scale(t(heatmap_norm[label_de_contrast8, to_keep_contrast8])))
contrast8_de_heatmap <- 
 Heatmap(heatmap_dat_contrast8,
         col = colorRamp2(c(seq(-3, 3, 0.05)), 
                          c(colorRampPalette(c("#0092b5", "white", "#a6ce39"))(121))),
         name = 'z-score', use_raster = TRUE,
         clustering_distance_rows = "pearson", clustering_method_rows = "average",
         clustering_distance_columns = "pearson", clustering_method_columns = "average",
         row_split = 2, column_split = 2,
         border_gp = gpar(col = "darkgray"),
         show_row_names = FALSE, show_column_names = FALSE,
         right_annotation = 
          rowAnnotation(foo = anno_mark(at = match(label_de_contrast8, rownames(heatmap_dat_contrast8)),
                                        labels = label_de_contrast8)),
         top_annotation = 
          HeatmapAnnotation(df=dfs_contrast8[[3]][factors_of_interest],
                            col = pal_main))
svglite::svglite(filename = file.path(de_dir, "contrast8_de_heatmap.svg"),
                width=10, height=10)
draw(contrast8_de_heatmap, merge_legend = TRUE, heatmap_legend_side = "right",
     annotation_legend_side = "right", adjust_annotation_extension = TRUE)
dev.off()

#png(filename = file.path(de_dir, "contrast8_de_heatmap.png"), width=10, height=10)
#draw(contrast8_de_heatmap, merge_legend = TRUE, heatmap_legend_side = "right", 
#      annotation_legend_side = "right", adjust_annotation_extension = TRUE)
#dev.off()

#```
#:::::::
#::::::: {#de_violin_contrast8_FOI .tabcontent}
#<!-- Make a series of violin plots of the label_de7 features. -->
#```{r de_violin_contrast8_FOI, eval = TRUE, echo = TRUE, fig.width=14, fig.height=10}
contrast_factor <- "PAM50"
comparison <- "Her2-enriched vs LumA"
p_adjust_class <- "FDR"
p_adjust_method <- "FDR"
violin_df <- cbind(dfs_contrast8[[3]]%>% dplyr::select(eval(contrast_factor)),
                   dfs_contrast8[[2]])
violin_df <- violin_df %>% tidyr::pivot_longer(cols=-1, names_to = "Protein", values_to = "Expression")
colnames(violin_df)[1] <- "contrast_factor"
violin_df <- violin_df %>% filter(Protein %in% label_de_contrast8)

# Format summary stats
violin_p_df <- de_contrast8 %>% filter(isSingular==0, beta!="(Intercept)")
violin_p_df$Comparison <- comparison
violin_p_df <- violin_p_df %>% tidyr::separate(col = Comparison, into=c("group1", "group2"), sep=" vs ")
violin_p_df["FDR"] <- signif(violin_p_df["FDR"], 3)
violin_exp_max <- ddply(violin_df, .(Protein), summarize,
                        y.position=((max(Expression)+1)*1.1)) # +1 for safe log2
violin_p_df <- base::merge(violin_p_df, violin_exp_max, by="Protein")
p <- ggplot(violin_df, 
            aes(x=contrast_factor, y=Expression, fill=contrast_factor)) + 
  geom_violin() +
  geom_jitter(width=0.25, height=0, size = 0.8, alpha=0.5) + 
  scale_fill_manual(values = pal_main[[contrast_factor]]) +
  facet_wrap(~Protein, scales = "free_y") +
  labs(x = eval(contrast_factor), y = "Expression (normalized counts)") +
  scale_y_continuous(trans = "log2", expand = expansion(mult = 0.2)) +
  theme_bw(base_size = 14) + 
  guides(fill=guide_legend(title = eval(contrast_factor)))
p <- p + ggprism::add_pvalue(
  violin_p_df,
  label="{p_adjust_class} = {violin_p_df[[p_adjust_method]]}", label.size = 3.6,
  y.position = violin_p_df$y.position
) + theme(legend.position="bottom")
ggsave(p, filename = file.path(de_dir, "contrast8_violins.svg"), width=14, height=10)
#ggsave(p, filename = file.path(de_dir, "contrast8_violins.png"), width=14, height=10)
p

#```
#:::::::
#:::::
#:::
#The searchable table below lists log~2~ fold change estimates and p-values for each protein.
#```{r dt results contrast8}
dt_params$buttons <-   list(list(extend = "copy"),
                       list(extend = "csv", filename = "contrast8"),
                       list(extend = "excel", filename = "contrast8"))
dt_contrast8 <- de_contrast8
DT::datatable(
  dt_contrast8[, c(1,2,6, 9)] %>% arrange(pval),
  extensions = c("Buttons", "Scroller", "FixedColumns"),
  options = dt_params,
  rownames = FALSE) %>% DT::formatRound(columns=3, digits=3) %>% DT::formatSignif(columns=2:4)
#```

message("================================================")
message("== Question  9:** Basal-like vs luminal A       ")
message("================================================")

#```
## Contrast 9: TN vs lumA
#We will now compare TN vs LumA PanCK- AOIs. The following formula was used to model differences for a given protein:
#$$ protein \sim PAM50 + (1|Scan\_ID) $$
#We adjust for the multiple sampling of ROI segments per patient with the $Scan\_ID$ variable.
#```{r de_contrast9, message=FALSE, warning=FALSE, echo=TRUE}
dfs_contrast9 <- dfs
to_keep_contrast9 <- row.names(dfs[[3]] %>% filter(PAM50 %in% c("Basal-like","LumA"),
                                                   PanCK=="negative"))
message("== Basal-like vs luminal A dim ==")
dim(to_keep_contrast9)
dfs_contrast9 <- filter_samples(DFS=dfs_contrast9, samples_to_keep = to_keep_contrast9)
dt9 <- ddply(dfs_contrast9[[3]], .(PAM50), summarize, 
       nPatients=length(unique(Scan_ID)), 
       nAOIs=length(PAM50))
knitr::kable(dt9, format="html")

de_contrast9 <- compute_lmer(d=dfs_contrast9, the_formula='protein_i ~ 1 + fixed + (1 |random)  ', 
                         the_factor="PAM50", 
                         the_random = "Scan_ID",
                         the_baseline="LumA")
de_contrast9 <- de_contrast9 %>% filter(isSingular==0, beta!="(Intercept)")
de_contrast9$FDR <- p.adjust(de_contrast9$pval, method = "BH")

#Proteins to label in plots
label_de_contrast9 <- de_contrast9[de_contrast9$pval<0.05,]$Protein
#```
#
#> Analyst notes: 19 proteins are significantly enriched at an FDR threshold of 0.05. 23 proteins are significantly enriched at a nominal P-value < 0.05. All proteins where P < 0.05 are labeled in the following figures.

#::: {#group_de_contrast9 .tabgroup}
#::::: {.tab}
#<button class="tablinks active" onclick="unrolltab(event, 'de_volcano_contrast9', 'group_de_contrast9')">Volcano</button>
#<button class="tablinks" onclick="unrolltab(event, 'de_heatmap_contrast9', 'group_de_contrast9')">Heatmap</button>
#<button class="tablinks" onclick="unrolltab(event, 'de_violin_contrast9_FOI', 'group_de_contrast9')">Violins</button>
#::::::: {#de_volcano_contrast9 .tabcontent style="display: block"}
#```{r de_volcano_contrast9,  fig.width=8, fig.height=6}
de_results_contrast9 <- v_plot2(df= de_contrast9,
        to_label=label_de_contrast9, 
        label_color="black", 
        fc_threshold = 0.50,
        to_lowlight = c(),
        lowlight_color = "grey", #Change this color to highlight any proteins of interest 
        WIDTH=9, HEIGHT=8, SCALE=1.4, LWD=1, out_dir=de_dir,
        the_name="PAM50 (Basal-like vs. LumA)",
        left="LumA", right="Basal-like",
        sizes=medians %>% dplyr::select(protein, Median)
        )
#de_results_contrast9[[1]]
#```
#
#:::::::
#::::::: {#de_heatmap_contrast9 .tabcontent}
#```{r de_heatmap_contrast9,  echo=TRUE, message=FALSE, warning=FALSE, collapse=TRUE, fig.width=10, fig.height=10, eval = TRUE}
heatmap_contrast9_features <- de_contrast9
heatmap_contrast9_features <- heatmap_contrast9_features$Protein
heatmap_norm <- as.data.frame(heatmap_dat)
heatmap_dat_contrast9 <- t(scale(t(heatmap_norm[label_de_contrast9, to_keep_contrast9])))
contrast9_de_heatmap <- 
 Heatmap(heatmap_dat_contrast9,
         col = colorRamp2(c(seq(-3, 3, 0.05)), 
                          c(colorRampPalette(c("#0092b5", "white", "#a6ce39"))(121))),
         name = 'z-score', use_raster = TRUE,
         clustering_distance_rows = "pearson", clustering_method_rows = "average",
         clustering_distance_columns = "pearson", clustering_method_columns = "average",
         row_split = 2, column_split = 2,
         border_gp = gpar(col = "darkgray"),
         show_row_names = FALSE, show_column_names = FALSE,
         right_annotation = 
          rowAnnotation(foo = anno_mark(at = match(label_de_contrast9, rownames(heatmap_dat_contrast9)),
                                        labels = label_de_contrast9)),
         top_annotation = 
          HeatmapAnnotation(df=dfs_contrast9[[3]][factors_of_interest],
                            col = pal_main))
svglite::svglite(filename = file.path(de_dir, "contrast9_de_heatmap.svg"),
                width=10, height=10)
draw(contrast9_de_heatmap, merge_legend = TRUE, heatmap_legend_side = "right",
     annotation_legend_side = "right", adjust_annotation_extension = TRUE)
dev.off()

#png(filename = file.path(de_dir, "contrast9_de_heatmap.png"), width=10, height=10)
#draw(contrast9_de_heatmap, merge_legend = TRUE, heatmap_legend_side = "right", 
#      annotation_legend_side = "right", adjust_annotation_extension = TRUE)
#dev.off()

#```
#:::::::
#::::::: {#de_violin_contrast9_FOI .tabcontent}
#<!-- Make a series of violin plots of the label_de7 features. -->
#```{r de_violin_contrast9_FOI, eval = TRUE, echo = TRUE, fig.width=14, fig.height=10}
contrast_factor <- "PAM50"
comparison <- "Basal-like vs LumA"
p_adjust_class <- "FDR"
p_adjust_method <- "FDR"
violin_df <- cbind(dfs_contrast9[[3]]%>% dplyr::select(eval(contrast_factor)),
                   dfs_contrast9[[2]])
violin_df <- violin_df %>% tidyr::pivot_longer(cols=-1, names_to = "Protein", values_to = "Expression")
colnames(violin_df)[1] <- "contrast_factor"
violin_df <- violin_df %>% filter(Protein %in% label_de_contrast9)

# Format summary stats
violin_p_df <- de_contrast9 %>% filter(isSingular==0, beta!="(Intercept)")
violin_p_df$Comparison <- comparison
violin_p_df <- violin_p_df %>% tidyr::separate(col = Comparison, into=c("group1", "group2"), sep=" vs ")
violin_p_df["FDR"] <- signif(violin_p_df["FDR"], 3)
violin_exp_max <- ddply(violin_df, .(Protein), summarize,
                        y.position=((max(Expression)+1)*1.1)) # +1 for safe log2
violin_p_df <- base::merge(violin_p_df, violin_exp_max, by="Protein")
p <- ggplot(violin_df, 
            aes(x=contrast_factor, y=Expression, fill=contrast_factor)) + 
  geom_violin() +
  geom_jitter(width=0.25, height=0, size = 0.8, alpha=0.5) + 
  scale_fill_manual(values = pal_main[[contrast_factor]]) +
  facet_wrap(~Protein, scales = "free_y") +
  labs(x = eval(contrast_factor), y = "Expression (normalized counts)") +
  scale_y_continuous(trans = "log2", expand = expansion(mult = 0.2)) +
  theme_bw(base_size = 14) + 
  guides(fill=guide_legend(title = eval(contrast_factor)))
p <- p + ggprism::add_pvalue(
  violin_p_df,
  label="{p_adjust_class} = {violin_p_df[[p_adjust_method]]}", label.size = 3.6,
  y.position = violin_p_df$y.position
) + theme(legend.position="bottom")
ggsave(p, filename = file.path(de_dir, "contrast9_violins.svg"), width=14, height=10)
#ggsave(p, filename = file.path(de_dir, "contrast9_violins.png"), width=14, height=10)
p

#```
#:::::::
#:::::
#:::
#The searchable table below lists log~2~ fold change estimates and p-values for each protein.
#```{r dt results contrast9}
dt_params$buttons <-   list(list(extend = "copy"),
                       list(extend = "csv", filename = "contrast9"),
                       list(extend = "excel", filename = "contrast9"))
dt_contrast9 <- de_contrast9
DT::datatable(
  dt_contrast9[, c(1,2,6, 9)] %>% arrange(pval),
  extensions = c("Buttons", "Scroller", "FixedColumns"),
  options = dt_params,
  rownames = FALSE) %>% DT::formatRound(columns=3, digits=3) %>% DT::formatSignif(columns=2:4)
#```

message("========================================")
message("== Question 12:** Younger vs older women")
message("========================================")

## Contrast 12: Younger vs older women
#We will now compare CK- AOIs from younger vs older women. The following formula was used to model differences for a given protein:
#$$ protein \sim Age + (1|Scan\_ID) $$
#We adjust for the multiple sampling of ROI segments per patient with the $Scan\_ID$ variable.
#```{r de_contrast12, message=FALSE, warning=FALSE, echo=TRUE}
dfs_contrast12 <- dfs
to_keep_contrast12 <- row.names(dfs[[3]] %>% filter( Age %in% c("young","old"),
                                                     PanCK=="negative"))
dfs_contrast12 <- filter_samples(DFS=dfs_contrast12, samples_to_keep = to_keep_contrast12)
dt12 <- ddply(dfs_contrast12[[3]], .(Age), summarize, 
       nPatients=length(unique(Scan_ID)), 
       nAOIs=length(Age))
knitr::kable(dt12, format="html")

de_contrast12 <- compute_lmer(d=dfs_contrast12, the_formula='protein_i ~ 1 + fixed + (1 |random)  ', 
                         the_factor="Age", 
                         the_random = "Scan_ID",
                         the_baseline="old")
de_contrast12 <- de_contrast12 %>% filter(isSingular==0, beta!="(Intercept)")
de_contrast12$FDR <- p.adjust(de_contrast12$pval, method = "BH")

#Proteins to label in plots
label_de_contrast12 <- de_contrast12[de_contrast12$pval<0.05,]$Protein

#```
#
#> Analyst notes: No proteins are significantly enriched at an FDR threshold of 0.05. 11 proteins are significantly enriched at a nominal P-value < 0.05. Only proteins where P < 0.05 are labeled in the following figures.

#::: {#group_de_contrast12 .tabgroup}
#::::: {.tab}
#<button class="tablinks active" onclick="unrolltab(event, 'de_volcano_contrast12', 'group_de_contrast12')">Volcano</button>
#<button class="tablinks" onclick="unrolltab(event, 'de_heatmap_contrast12', 'group_de_contrast12')">Heatmap</button>
#<button class="tablinks" onclick="unrolltab(event, 'de_violin_contrast12_FOI', 'group_de_contrast12')">Violins</button>
#::::::: {#de_volcano_contrast12 .tabcontent style="display: block"}
#```{r de_volcano_contrast12,  fig.width=8, fig.height=6}
de_results_contrast12 <- v_plot2(df= de_contrast12,
        to_label=label_de_contrast12, 
        label_color="black", 
        fc_threshold = 0.50,
        to_lowlight = c(),
        lowlight_color = "grey", #Change this color to highlight any proteins of interest 
        WIDTH=9, HEIGHT=8, SCALE=1.4, LWD=1, out_dir=de_dir,
        the_name="Age (Young vs. Old)",
        left="Old", right="Young",
        sizes=medians %>% dplyr::select(protein, Median)
        )
#de_results_contrast12[[1]]

#```
#
#:::::::
#::::::: {#de_heatmap_contrast12 .tabcontent}
#
#```{r de_heatmap_contrast12,  echo=TRUE, message=FALSE, warning=FALSE, collapse=TRUE, fig.width=10, fig.height=10, eval = TRUE}
heatmap_contrast12_features <- de_contrast12
heatmap_contrast12_features <- heatmap_contrast12_features$Protein
heatmap_norm <- as.data.frame(heatmap_dat)
heatmap_dat_contrast12 <- t(scale(t(heatmap_norm[label_de_contrast12, to_keep_contrast12])))
contrast12_de_heatmap <- 
 Heatmap(heatmap_dat_contrast12,
         col = colorRamp2(c(seq(-3, 3, 0.05)), 
                          c(colorRampPalette(c("#0092b5", "white", "#a6ce39"))(121))),
         name = 'z-score', use_raster = TRUE,
         clustering_distance_rows = "pearson", clustering_method_rows = "average",
         clustering_distance_columns = "pearson", clustering_method_columns = "average",
         row_split = 2, column_split = 2,
         border_gp = gpar(col = "darkgray"),
         show_row_names = FALSE, show_column_names = FALSE,
         right_annotation = 
          rowAnnotation(foo = anno_mark(at = match(label_de_contrast12, rownames(heatmap_dat_contrast12)),
                                        labels = label_de_contrast12)),
         top_annotation = 
          HeatmapAnnotation(df=dfs_contrast12[[3]][factors_of_interest],
                            col = pal_main))
svglite::svglite(filename = file.path(de_dir, "contrast12_de_heatmap.svg"),
                width=10, height=10)
draw(contrast12_de_heatmap, merge_legend = TRUE, heatmap_legend_side = "right",
     annotation_legend_side = "right", adjust_annotation_extension = TRUE)
dev.off()

#draw(contrast12_de_heatmap, merge_legend = TRUE, heatmap_legend_side = "right", 
#      annotation_legend_side = "right", adjust_annotation_extension = TRUE)

#```
#
#:::::::
#::::::: {#de_violin_contrast12_FOI .tabcontent}
#<!-- Make a series of violin plots of the label_de12 features. -->
#```{r de_violin_contrast12_FOI, eval = TRUE, echo = TRUE, fig.width=12, fig.height=8}

contrast_factor <- "Age"
comparison <- "young vs old"
p_adjust_class <- "FDR"
p_adjust_method <- "FDR"
violin_df <- cbind(dfs_contrast12[[3]]%>% dplyr::select(eval(contrast_factor)),
                   dfs_contrast12[[2]])
violin_df <- violin_df %>% tidyr::pivot_longer(cols=-1, names_to = "Protein", values_to = "Expression")
colnames(violin_df)[1] <- "contrast_factor"
violin_df <- violin_df %>% filter(Protein %in% label_de_contrast12)


# Format summary stats
violin_p_df <- de_contrast12 %>% filter(isSingular==0, beta!="(Intercept)")
violin_p_df$Comparison <- comparison
violin_p_df <- violin_p_df %>% tidyr::separate(col = Comparison, into=c("group1", "group2"), sep=" vs ")
violin_p_df["FDR"] <- signif(violin_p_df["FDR"], 3)
violin_exp_max <- ddply(violin_df, .(Protein), summarize,
                        y.position=((max(Expression)+1)*1.1)) # +1 for safe log2
violin_p_df <- base::merge(violin_p_df, violin_exp_max, by="Protein")
p <- ggplot(violin_df, 
            aes(x=contrast_factor, y=Expression, fill=contrast_factor)) + 
  geom_violin() +
  geom_jitter(width=0.25, height=0, size = 0.8, alpha=0.5) + 
  scale_fill_manual(values = pal_main[[contrast_factor]]) +
  facet_wrap(~Protein, scales = "free_y") +
  labs(x = eval(contrast_factor), y = "Expression (normalized counts)") +
  scale_y_continuous(trans = "log2", expand = expansion(mult = 0.2)) +
  theme_bw(base_size = 14) + 
  guides(fill=guide_legend(title = eval(contrast_factor)))
p <- p + ggprism::add_pvalue(
  violin_p_df,
  label="{p_adjust_class} = {violin_p_df[[p_adjust_method]]}", label.size = 3.6,
  y.position = violin_p_df$y.position
) + theme(legend.position="bottom")
ggsave(p, filename = file.path(de_dir, "contrast12_violins.svg"), width=12, height=8)
p

#```
#:::::::
#:::::
#:::
#The searchable table below lists log~2~ fold change estimates and p-values for each protein.
#```{r dt results contrast12}

dt_params$buttons <-   list(list(extend = "copy"),
                       list(extend = "csv", filename = "contrast12"),
                       list(extend = "excel", filename = "contrast12"))
dt_contrast12 <- de_contrast12
DT::datatable(
  dt_contrast12[, c(1,2,6, 9)] %>% arrange(pval),
  extensions = c("Buttons", "Scroller", "FixedColumns"),
  options = dt_params,
  rownames = FALSE) %>% DT::formatRound(columns=3, digits=3) %>% DT::formatSignif(columns=2:4)

#```
#
# Data and R Information
#Save analysis image and print R session information for reference:
#```{r report_close}
save.image(paste0(file=date_tag, ".RData"))
sessionInfo()
#```
