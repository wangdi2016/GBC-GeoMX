###
source("./step4.dea.x3.init.v2.R")

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

message("== PanCK+ vs PanCK- ==")
dim(to_keep_contrast1)

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
        the_name="Contrast 1 PanCK+ vs PanCK-",
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

##
pdf(file.path(de_dir, "contrast1_de_heatmap.pdf"), 
    width = 10, height = 12)

draw(contrast1_de_heatmap, 
     merge_legend = TRUE, 
     heatmap_legend_side = "right",
     annotation_legend_side = "right", 
     adjust_annotation_extension = TRUE)

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
ggsave(p, filename = file.path(de_dir, "contrast1_violins.pdf"), width=18, height=12)
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

##
sessionInfo()
