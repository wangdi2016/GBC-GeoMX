###
#source("./init.R")
###
source("./step4.dea.x3.init.v2.R")

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
# **Question 1:** CK+ vs CK- AOIs
#
# **Question 2:** IM vs CC
#
# **Question 3:** LGD vs IM
#
# **Question 4:** HGD vs LGD
#
# **Question 5:** IC vs HGD
# 
message("========================================")
message("== Question  5:** IC vs HGD             ")
message("========================================")

#```
## Contrast 5: IC vs HGD
#We will now compare IC vs HGD for CK+ and CK- AOIs. The following formula was used to model differences for a given protein:
#$$ protein \sim HGD + (1|Scan\_ID) $$
#We adjust for the multiple sampling of ROI segments per patient with the $Scan\_ID$ variable.
#```{r de_contrast5, message=FALSE, warning=FALSE, echo=TRUE}
dfs_contrast5 <- dfs
#to_keep_contrast5 <- row.names(dfs[[3]] %>% filter(CancerType %in% c("IM","CC"), CK=="negative"))
to_keep_contrast5 <- row.names(dfs[[3]] %>% filter(CancerType %in% c("IC","HGD")))
dfs_contrast5 <- filter_samples(DFS=dfs_contrast5, samples_to_keep = to_keep_contrast5)
dt5 <- ddply(dfs_contrast5[[3]], .(CancerType), summarize, 
       nPatients=length(unique(Scan_ID)), 
       nAOIs=length(CancerType))
knitr::kable(dt5, format="html")

de_contrast5 <- compute_lmer(d=dfs_contrast5, the_formula='protein_i ~ 1 + fixed + (1 |random)  ', 
                         the_factor="CancerType", 
                         the_random = "Scan_ID",
                         the_baseline="HGD")
de_contrast5 <- de_contrast5 %>% filter(isSingular==0, beta!="(Intercept)")
de_contrast5$FDR <- p.adjust(de_contrast5$pval, method = "BH")

#Proteins to label in plots
label_de_contrast5 <- de_contrast5[de_contrast5$pval<0.05,]$Protein
#```
#
#> Analyst notes: 19 proteins are significantly enriched at an FDR threshold of 0.05. 23 proteins are significantly enriched at a nominal P-value < 0.05. All proteins where P < 0.05 are labeled in the following figures.

#::: {#group_de_contrast5 .tabgroup}
#::::: {.tab}
#<button class="tablinks active" onclick="unrolltab(event, 'de_volcano_contrast5', 'group_de_contrast5')">Volcano</button>
#<button class="tablinks" onclick="unrolltab(event, 'de_heatmap_contrast5', 'group_de_contrast5')">Heatmap</button>
#<button class="tablinks" onclick="unrolltab(event, 'de_violin_contrast5_FOI', 'group_de_contrast5')">Violins</button>
#::::::: {#de_volcano_contrast5 .tabcontent style="display: block"}
#```{r de_volcano_contrast5,  fig.width=8, fig.height=6}
de_results_contrast5 <- v_plot2(df= de_contrast5,
        to_label=label_de_contrast5, 
        label_color="black", 
        fc_threshold = 0.50,
        to_lowlight = c(),
        lowlight_color = "grey", #Change this color to highlight any proteins of interest 
        WIDTH=9, HEIGHT=8, SCALE=1.4, LWD=1, out_dir=de_dir,
        the_name="Contrast 5 IC vs HGD",
        left="HGD", right="IC",
        sizes=medians %>% dplyr::select(protein, Median)
        )
de_results_contrast5[[1]]
#```
#
#:::::::
#::::::: {#de_heatmap_contrast5 .tabcontent}
#```{r de_heatmap_contrast5,  echo=TRUE, message=FALSE, warning=FALSE, collapse=TRUE, fig.width=10, fig.height=10, eval = TRUE}
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

## 
pdf(file.path(de_dir, "contrast5_de_heatmap.pdf"), 
    width = 10, height = 10)

draw(contrast5_de_heatmap, 
     merge_legend = TRUE, 
     heatmap_legend_side = "right",
     annotation_legend_side = "right", 
     adjust_annotation_extension = TRUE)

dev.off()

#```
#:::::::
#::::::: {#de_violin_contrast5_FOI .tabcontent}
#<!-- Make a series of violin plots of the label_de7 features. -->
#```{r de_violin_contrast5_FOI, eval = TRUE, echo = TRUE, fig.width=14, fig.height=10}
contrast_factor <- "CancerType"
comparison <- "IC vs HGD"
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
ggsave(p, filename = file.path(de_dir, "contrast5_violins.svg"), width=14, height=10)
ggsave(p, filename = file.path(de_dir, "contrast5_violins.pdf"), width=14, height=10)
p

#```
#:::::::
#:::::
#:::
#The searchable table below lists log~2~ fold change estimates and p-values for each protein.
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
#```

# Data and R Information
#Save analysis image and print R session information for reference:
#```{r report_close}
save.image(paste0(file=date_tag, ".RData"))
sessionInfo()
#```
