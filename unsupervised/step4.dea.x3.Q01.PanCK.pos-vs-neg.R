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

## load packages
library(dplyr)
library(broom)

##
# force CancerType in this order
desired_order <- c("CC", "IM", "LGD", "HGD", "IC")
contrast_factor <- "CancerType"

violin_df <- cbind(dfs_contrast1[[3]]%>% dplyr::select(eval(contrast_factor)),
                   dfs_contrast1[[2]])
violin_df <- violin_df %>% tidyr::pivot_longer(cols=-1, names_to = "Protein", values_to = "Expression")
#colnames(violin_df)[1] <- "contrast_factor"
violin_df <- violin_df %>% filter(Protein %in% label_de_contrast1)

head(violin_df)

#
# Force factor levels to the desired order
violin_df[[contrast_factor]] <- factor(
  violin_df[[contrast_factor]],
  levels = desired_order
)

# Make sure CancerType is a factor with the baseline you want
violin_df$CancerType <- factor(violin_df$CancerType,
                               levels = c("CC", "IM", "LGD", "HGD", "IC"))

# Run linear regression per protein
lm_coeffs <- violin_df %>%
  group_by(Protein) %>%
  do(tidy(lm(Expression ~ CancerType, data = .))) %>%
  ungroup()

# Inspect
head(lm_coeffs)
write.csv(lm_coeffs, "lm_coeffs.csv", row.names = FALSE)


# Encode CancerType numerically for trend
violin_df$CancerType_num <- as.numeric(factor(violin_df$CancerType,
                                              levels = c("CC", "IM", "LGD", "HGD", "IC")))
corr_coeffs <- violin_df %>%
  group_by(Protein) %>%
  summarise(
    y_pos = max(Expression, na.rm = TRUE) + 50,  # top of y-axis for this protein
    cor_coef = cor(Expression, CancerType_num, method = "pearson"),
    .groups = "drop"
  )

write.csv(corr_coeffs, "correlation_coeffs.csv", row.names = FALSE)

plot_df <- violin_df %>%
  left_join(corr_coeffs, by = "Protein")

### make plot
p2 <- ggplot(plot_df,
            aes(x=.data[[contrast_factor]], y=Expression, fill=.data[[contrast_factor]])) +
  geom_violin() +
  geom_jitter(width=0.25, height=0, size = 0.8, alpha=0.3) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "red", linetype = "dashed") +
  scale_fill_manual(values = pal_main[[contrast_factor]]) +
  #geom_text(aes(x = 3, y = max(Expression, na.rm = TRUE),
  geom_text(data = corr_coeffs, aes(x = 3, y = y_pos,
                label = paste0("r = ", round(cor_coef, 2))),
                inherit.aes = FALSE, hjust = 0.5) +
  facet_wrap(~Protein, scales = "free_y", ncol = 4) +
  labs(x = eval(contrast_factor), y = "Expression (normalized counts)") +
  scale_y_continuous(trans = "log2", expand = expansion(mult = 0.2)) +
  theme_bw(base_size = 14) +
  guides(fill=guide_legend(title = eval(contrast_factor))) +
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 45, hjust = 1) 
  )

ggsave(p2, filename = file.path(de_dir, "contrast1_violins.byCancerType.svg"), width=18, height=32)
ggsave(p2, filename = file.path(de_dir, "contrast1_violins.byCancerType.pdf"), width=18, height=32)

p2

## make new plot

## make plot without labels
##
p3 <- ggplot(violin_df,
          aes(x=.data[[contrast_factor]], y=Expression, fill=.data[[contrast_factor]])) +
  geom_violin() +
  geom_jitter(width=0.25, height=0, size = 0.8, alpha=0.3) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "red", linetype = "dashed") +
  scale_fill_manual(values = pal_main[[contrast_factor]]) +
#  facet_wrap(~Protein, scales = "free_y", ncol = 4)
  facet_wrap(~Protein, scales = "free_y", ncol = 8)

# Build the ggplot object
gb <- ggplot_build(p3)

# Extract panel y-limits
ymax_df <- data.frame(
  Protein = levels(factor(violin_df$Protein)),
  y_max = sapply(gb$layout$panel_scales_y, function(x) x$range$range[2])
)

ymin_df <- data.frame(
  Protein = levels(factor(violin_df$Protein)),
  y_min = sapply(gb$layout$panel_scales_y, function(x) x$range$range[1])
)

## Combine with correlation coefficients
cor_df <- violin_df %>%
  group_by(Protein) %>%
  summarise(cor_coef = cor(Expression, as.numeric(factor(CancerType,
                                                         levels = c("CC","IM","LGD","HGD","IC")))),
            .groups = "drop")

label_df1 <- left_join(cor_df, ymax_df, by = "Protein")
label_df2 <- left_join(cor_df, ymin_df, by = "Protein")

## Add text using computed y_max
p3 <- p3 + geom_text(data = label_df2,
              aes(x = 1, y = y_min*0.5, label = paste0("r = ", sprintf("%.2f", cor_coef))),
              inherit.aes = FALSE, hjust = 0.5) +
  labs(x = eval(contrast_factor), y = "Expression (normalized counts)") +
  scale_y_continuous(trans = "log2", expand = expansion(mult = 0.2)) +
  theme_bw(base_size = 14) +
  guides(fill=guide_legend(title = eval(contrast_factor))) +
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))

#ggsave(p3, filename = file.path(de_dir, "contrast1_violins.byCancerType.v3.svg"), width=18, height=32)
ggsave(p3, filename = file.path(de_dir, "contrast1_violins.byCancerType.v3.svg"), width=36, height=16)
#ggsave(p3, filename = file.path(de_dir, "contrast1_violins.byCancerType.v3.pdf"), width=18, height=32)
ggsave(p3, filename = file.path(de_dir, "contrast1_violins.byCancerType.v3.pdf"), width=36, height=16)

###
# Numeric encoding for CancerType
violin_df$CancerType_num <- as.numeric(factor(violin_df$CancerType,
                                              levels = c("CC","IM","LGD","HGD","IC")))

# Correlation per protein
cor_df <- violin_df %>%
  group_by(Protein) %>%
  summarise(cor_coef = cor(Expression, CancerType_num), .groups = "drop")

# Step 3: Determine number of pages
library(ggforce)
n_proteins <- length(unique(violin_df$Protein))
per_page <- 8   # 4 proteins per page
n_pages <- ceiling(n_proteins / per_page)

pdf(file.path(de_dir, "contrast1_violins.byCancerType.v4.pdf"), width = 12, height = 6)

# Loop over pages and extract y-max per page
for (i in 1:n_pages) {

p4 <- ggplot(violin_df, aes(x = CancerType, y = Expression, fill = CancerType)) +
  geom_violin() +
  geom_jitter(width=0.25, height=0, size = 0.8, alpha=0.3) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "red", linetype = "dashed") +
  facet_wrap_paginate(~Protein, scales = "free_y", ncol = 4, nrow = 2, page = i)

# Extract y-min per facet in this page
  gb <- ggplot_build(p4)
  panel_proteins <- gb$layout$panel_params[[1]]$strip_text  # list of proteins on this page
  y_min <- sapply(gb$layout$panel_scales_y, function(x) x$range$range[1]) ## x$range$range[2] for ymax
  
  # Prepare label dataframe for this page
  label_df <- data.frame(
    Protein = unique(violin_df$Protein)[gb$layout$layout$ROW],  # match panel order
    y_min = y_min
  ) %>% 
    left_join(cor_df, by = "Protein")
  
# Add correlation text
p4 <- p4 + geom_text(data = label_df,
                     aes(x = 1, y = y_min*0.5,
                         label = paste0("r = ", sprintf("%.2f", cor_coef))),
                     inherit.aes = FALSE, hjust = 0.5) +
  labs(x = eval(contrast_factor), y = "Expression (normalized counts)") +
  scale_y_continuous(trans = "log2", expand = expansion(mult = 0.2)) +
  theme_bw(base_size = 14) +
  guides(fill=guide_legend(title = eval(contrast_factor))) +
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))

print(p4)
}

dev.off()

###
### make px ggplot obj then split it
###

px <- ggplot(violin_df,
          aes(x=.data[[contrast_factor]], y=Expression, fill=.data[[contrast_factor]])) +
  geom_violin() +
  geom_jitter(width=0.25, height=0, size = 0.8, alpha=0.3) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "red", linetype = "dashed") +
  scale_fill_manual(values = pal_main[[contrast_factor]])

# Build the ggplot object
gb <- ggplot_build(px)

# Extract panel y-limits
ymax_df <- data.frame(
  Protein = levels(factor(violin_df$Protein)),
  y_max = sapply(gb$layout$panel_scales_y, function(x) x$range$range[2])
)

ymin_df <- data.frame(
  Protein = levels(factor(violin_df$Protein)),
  y_min = sapply(gb$layout$panel_scales_y, function(x) x$range$range[1])
)

## Combine with correlation coefficients
cor_df <- violin_df %>%
  group_by(Protein) %>%
  summarise(cor_coef = cor(Expression, as.numeric(factor(CancerType,
                                                         levels = c("CC","IM","LGD","HGD","IC")))),
            .groups = "drop")

label_df1 <- left_join(cor_df, ymax_df, by = "Protein")
label_df2 <- left_join(cor_df, ymin_df, by = "Protein")

## Add text using computed y_max
px <- px + geom_text(data = label_df2,
              aes(x = 1.5, y = y_min*0.5, label = paste0("r = ", sprintf("%.2f", cor_coef))),
              inherit.aes = FALSE, hjust = 0.5) +
  labs(x = eval(contrast_factor), y = "Expression (normalized counts)") +
  scale_y_continuous(trans = "log2", expand = expansion(mult = 0.2), labels = function(x) sprintf("%.1f", x)) +
  theme_bw(base_size = 14) +
  guides(fill=guide_legend(title = eval(contrast_factor))) +
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(px, filename = file.path(de_dir, "contrast1_violins.byCancerType.v5.pdf"), width=12, height=6)

##
n_proteins <- length(unique(violin_df$Protein))
per_page <- 8  # number of proteins per page
n_pages <- ceiling(n_proteins / per_page)

pdf(file.path(de_dir, "contrast1_violins.byCancerType.vx.pdf"), width = 12, height = 6)  

for (i in 1:n_pages) {
  # Add pagination facet
  p_page <- px + facet_wrap_paginate(~Protein, scales = "free_y",
                                    ncol = 4, nrow = 2, page = i)
  # Optionally: add per-facet correlation labels or other annotations here
  print(p_page)
}

dev.off()


# ####################
# library(dplyr)
# library(broom)
# 
# # Make sure CancerType is a factor with the baseline you want
# violin_df$CancerType <- factor(violin_df$CancerType,
#                                levels = c("CC", "IM", "LGD", "HGD", "IC"))
# 
# # Run linear regression per protein
# lm_coeffs <- violin_df %>%
#   group_by(Protein) %>%
#   do(tidy(lm(Expression ~ CancerType, data = .))) %>%
#   ungroup()
# 
# # Inspect
# head(lm_coeffs)
# write.csv(lm_coeffs, "lm_coeffs.csv", row.names = FALSE)
# 
# 
# # Encode CancerType numerically for trend
# violin_df$CancerType_num <- as.numeric(factor(violin_df$CancerType,
#                                               levels = c("CC", "IM", "LGD", "HGD", "IC")))
# 
# corr_coeffs <- violin_df %>%
#   group_by(Protein) %>%
#   summarise(
#     cor_coef = cor(Expression, CancerType_num, method = "pearson"),
#     .groups = "drop"
#   )
# 
# write.csv(corr_coeffs, "correlation_coeffs.csv", row.names = FALSE)

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
