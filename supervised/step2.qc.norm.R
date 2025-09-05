###
source("./init.dw.R")

### load image from step1
load(file=file.path(object_dir,"step1.checkdata.RData"))

#######################################
# 2 Quality Control and Normalization #
#######################################
#
# We explore quality control data found in the Initial Dataset. We use plots to guide sample filtering.

# Every ROI/AOI segment was tested for:

# ERCC positive controls: Hybridization factors should be between 0.1 and 10. Samples will be normalized using these positive hybridization factors.
# Nuclei: >100 nuclei per segment is generally recommended; however, this cutoff is highly study/tissue dependent and may need to be reduced; what is most important is consistency in the nuclei distribution for segments within the study.
# Area: generally correlates with nuclei; a strict cutoff is not generally applied based on area.

##################
# 2.1 QC Summary #
##################
# One sample was filtered based on high ERCC positive control normalization factors
# 8 AOI segments were removed due to low nuclei counts and low area
# 3 AOI segments were removed due to high normalization factors
# 719 segments were profiled and 707 passed QC and were included in downstream analysis
# 44/52 proteins are retained for downstream analysis
# The 8 proteins removed are ERCC related, background IgG controls, and housekeeper proteins
# Data underwent housekeeper normalization (S6 and Histone H3)

######################
# 2.2 ERCC normalize #
######################
# First, we calculate the hybridization factors and remove any sample with a value >= 10. Next, we perform ERCC normalization: normalize the data using the hybridization factors (only those less than 10). And finally, we remove ERCC-related columns.
message("== dim of raw ==")
dim(raw)
#raw

hyb_factors <- mean(raw[,which(colnames(raw)=="HYB-POS")])/raw[,which(colnames(raw)=="HYB-POS")]
pruned <- raw[hyb_factors<10,]
message("== dim of kept ==")
dim(pruned)
#pruned
colnames(pruned)

removed <- raw[hyb_factors>=10,]
message("== dim of removed ==")
dim(removed)
head(removed)
cat(colnames(removed), sep = "\n")

ercc_norm <- sweep(pruned, 1, hyb_factors[hyb_factors<10], "*")
#ercc_norm_removed <- sweep(removed, 1, hyb_factors[hyb_factors>=10], "*")

message("== Kept after ERCC Normalization ==")
ercc_norm <- ercc_norm[,!grepl("^HYB", colnames(ercc_norm))]
dim(ercc_norm)
cat(colnames(ercc_norm), sep = "\n")

#message("== removed after ERCC Normalization ==")
#dim(ercc_norm_removed)
#ercc_norm_removed

message("== removed HYB-POS and HYB-NEG proteins from 60 protein list")
message("== End of ERCC Normalization ==")

##################
# Analyst notes: #
##################
# There were 719 AOI segments before hybridization pruning.
# There are 718 AOI segments after hybridization pruning.
# There are now 50 molecules including control and HK molecules.
# One sample was removed due to high positive hybridization factors.

######################################################################
# 2.3 Quality control of technical factors - AOI Area & Nuclei Count #
######################################################################

#Remove the individual samples that did not pass ERCC normalization - in this case, no samples were removed.
annots <- annots[row.names(ercc_norm),] 

#Renaming and sorting objects.
dfs3 <- annots #annotations
dfs1 <- ercc_norm[rownames(annots), ] #protein counts, filtered for ERCC norm threshold

###########################
# 2.3.1 QC Visualizations #
###########################
# Analyst notes: 8 AOI segments are below the nuclei threshold of 25 and area threshold of 1000. 8 segments were removed during this step.

########################
# 2.3.1.1 Nuclei count #
########################

# Plotting nuclei count and coloring by factors of interest. Minimum nuclei = 25 (vertical dashed line).

g1 <- ggplot(dfs3, aes(x=Nuclei)) + 
  geom_histogram(aes(fill=PanCK), bins=100) + 
  theme_bw() + 
  geom_vline(xintercept = min_nc, color="grey", lty="dashed") + 
  xlab("Nuclei Counts") + ylab("Number of AOIs") +
  scale_fill_manual(values = pal_main[["PanCK"]])
g2 <- ggplot(dfs3, aes(x=Nuclei)) + 
  geom_histogram(aes(fill=CancerType), bins=100) + 
  theme_bw() + 
  geom_vline(xintercept = min_nc, color="grey", lty="dashed") + 
  xlab("Nuclei Counts") + ylab("Number of AOIs") +
  scale_fill_manual(values = pal_main[["CancerType"]])
g3 <- ggplot(dfs3, aes(x=Nuclei)) + 
  geom_histogram(aes(fill=EpithelialStatus), bins=100) + 
  theme_bw() + 
  geom_vline(xintercept = min_nc, color="grey", lty="dashed") + 
  xlab("Nuclei Counts") + ylab("Number of AOIs") +
  scale_fill_manual(values = pal_main[["EpithelialStatus"]])
#svg(filename= file.path(fig_dir, "QC_NucleiCounts.svg"), width=8, height=8)

g123 <- grid.arrange(g1,g2,g3, nrow=3)
#g123 <- grid.arrange(g1,g2, nrow=2)
#invisible(capture.output(dev.off())) 
#grid.arrange(g1,g2,g3, nrow=3)
ggsave(g123, filename= file.path(fig_dir, "QC_NucleiCounts.png"), width=8, height=8)

####################
# 2.3.1.2 AOI Area #
####################

# Plotting AOI area and coloring by factors of interest. Minimum area = 1000 (vertical dashed line).

g4 <- ggplot(dfs3, aes(x=Area)) + 
  geom_histogram(aes(fill=PanCK), bins=100) + 
  theme_bw() + 
  geom_vline(xintercept = min_area, color="grey", lty="dashed") + 
  xlab("Area") + ylab("Number of AOIs") +
  scale_fill_manual(values = pal_main[["PanCK"]])
g5 <- ggplot(dfs3, aes(x=Area)) + 
  geom_histogram(aes(fill=CancerType), bins=100) + 
  theme_bw() + 
  geom_vline(xintercept = min_area, color="grey", lty="dashed") + 
  xlab("Area") + ylab("Number of AOIs") +
  scale_fill_manual(values = pal_main[["CancerType"]])
g6 <- ggplot(dfs3, aes(x=Area)) + 
  geom_histogram(aes(fill=EpithelialStatus), bins=100) + 
  theme_bw() + 
  geom_vline(xintercept = min_area, color="grey", lty="dashed") + 
  xlab("Area") + ylab("Number of AOIs") +
  scale_fill_manual(values = pal_main[["EpithelialStatus"]])
#svg(filename= file.path(fig_dir, "QC_AOIArea.svg"), width=8, height=8)
#png(filename= file.path(fig_dir, "QC_AOIArea.png"), width=8, height=8)

g456 <- grid.arrange(g4,g5,g6, nrow=3)
ggsave(g456, filename= file.path(fig_dir, "QC_AOIArea.png"), width=8, height=8)

##
##
below_nc <- nrow(filter(dfs3, Nuclei<=min_nc)) #8
#message(paste0("Number of AOIs below the nuclei threshold of ", min_nc, ": ", below_nc, " (", round(100*below_nc/nrow(dfs3), 2), "%)."))
below_area <-nrow(filter(dfs3, Area<=1000)) #6
#message(paste0("Number of AOIs below the area of ", min_area, ": ", below_area, " (", round(100*below_nc/nrow(dfs3), 2), "%)."))
to_filter <-nrow(filter(dfs3, Nuclei <= min_nc | Area<=1000)) #8
message(paste0("Number of AOIs below the area of ", min_area, " and nuclei threshold of ",min_nc, ": ", to_filter, " (", round(100*to_filter/nrow(dfs3), 2), "%)."))

### 8 samples were removed during this analysis ###
### All samples withe low area also have low nuclei counts. Will remove samples with low nuclei counts.
dfs3$low_nc <- ifelse(dfs3$Nuclei<=min_nc, 1, 0)
#dfs3[(dfs3$low_nc == 1),]  

rownames(dfs3[(dfs3$low_nc == 1),])

##                                                        Segment       PanCK      Morphology TumorEdge_LE_vs_rest
## Slide-0989632_AKU040 T12-17 | LE-LQ1 | PanCK positive PanCK positive positive Combined stroma                other
## Slide-0989632_AKU040 T12-17 | MM-LQ1 | PanCK positive PanCK positive positive Combined stroma                other
## Slide-0989648_KAIC004 T1-16 | LE-LQ1 | PanCK positive PanCK positive positive Combined stroma                other
## Slide-0991905_KAIC017 T1-12 | LE-UQ2 | PanCK positive PanCK positive positive Combined stroma                other
## Slide-0991913_KAIC022 T2-12 | LE-LQ1 | PanCK positive PanCK positive positive Combined stroma                other
## Slide-0991913_KAIC022 T2-12 | LE-UQ2 | PanCK positive PanCK positive positive Combined stroma                other
## Slide-0991927_KAIC034 T1-10 | LE-LQ2 | PanCK positive PanCK positive positive Combined stroma                other
## Slide-0991931_KAIC007 T1-5 | MM-LQ1 | PanCK positive  PanCK positive positive Combined stroma                other
##                                                    Distance.to.the.centroid.of.ROIs.that.Xing.generated.and.dichotomized.using.median
## Slide-0989632_AKU040 T12-17 | LE-LQ1 | PanCK positive                                                                         0.67462547
## Slide-0989632_AKU040 T12-17 | MM-LQ1 | PanCK positive                                                                         0.99470400
## Slide-0989648_KAIC004 T1-16 | LE-LQ1 | PanCK positive                                                                         0.21460710
## Slide-0991905_KAIC017 T1-12 | LE-UQ2 | PanCK positive                                                                        -0.97240524
## Slide-0991913_KAIC022 T2-12 | LE-LQ1 | PanCK positive                                                                         0.01847085
## Slide-0991913_KAIC022 T2-12 | LE-UQ2 | PanCK positive                                                                         1.17638270
## Slide-0991927_KAIC034 T1-10 | LE-LQ2 | PanCK positive                                                                        -0.56434401
## Slide-0991931_KAIC007 T1-5 | MM-LQ1 | PanCK positive                                                                          1.08803036
##                                                    DistanceToCentroid PanCK_proximity CD45_proximity Luminal Grade
## Slide-0989632_AKU040 T12-17 | LE-LQ1 | PanCK positive               high            <NA>           <NA> non.Lum  high
## Slide-0989632_AKU040 T12-17 | MM-LQ1 | PanCK positive               high   nonCKenriched   CD45enriched     Lum  high
## Slide-0989648_KAIC004 T1-16 | LE-LQ1 | PanCK positive               high   nonCKenriched   CD45enriched     Lum   low
## Slide-0991905_KAIC017 T1-12 | LE-UQ2 | PanCK positive                low            <NA>           <NA>     Lum   low
## Slide-0991913_KAIC022 T2-12 | LE-LQ1 | PanCK positive               high   nonCKenriched   CD45enriched     Lum  high
## Slide-0991913_KAIC022 T2-12 | LE-UQ2 | PanCK positive               high   nonCKenriched   CD45enriched     Lum  high
## Slide-0991927_KAIC034 T1-10 | LE-LQ2 | PanCK positive                low   nonCKenriched   CD45enriched     Lum  high
## Slide-0991931_KAIC007 T1-5 | MM-LQ1 | PanCK positive                high   nonCKenriched   CD45enriched     Lum  high
##                                                    Parity.status(above.3.cut.off) Parity weight height BMI_value  BMI Waist_to_hip_ratio
## Slide-0989632_AKU040 T12-17 | LE-LQ1 | PanCK positive                              4   high    103  180.0  31.79012 high               high
## Slide-0989632_AKU040 T12-17 | MM-LQ1 | PanCK positive                              4   high    103  180.0  31.79012 high               high
## Slide-0989648_KAIC004 T1-16 | LE-LQ1 | PanCK positive                              4   high     79  152.5  33.96936 high               high
## Slide-0991905_KAIC017 T1-12 | LE-UQ2 | PanCK positive                              0    low     NA  149.3  27.36593 high                low
## Slide-0991913_KAIC022 T2-12 | LE-LQ1 | PanCK positive                              3    low     62  159.6  24.34030  low                low
## Slide-0991913_KAIC022 T2-12 | LE-UQ2 | PanCK positive                              3    low     62  159.6  24.34030  low                low
## Slide-0991927_KAIC034 T1-10 | LE-LQ2 | PanCK positive                              5   high     55  163.1  20.67540  low                low
## Slide-0991931_KAIC007 T1-5 | MM-LQ1 | PanCK positive                              10   high     50  149.0  22.52151  low                low
##                                                    Age_value   Age       Area Nuclei BindingDensity                     Scan_ID ROI_ID
## Slide-0989632_AKU040 T12-17 | LE-LQ1 | PanCK positive        50   old  318.12038     12           0.14 Slide-0989632_AKU040 T12-17 LE-LQ1
## Slide-0989632_AKU040 T12-17 | MM-LQ1 | PanCK positive        50   old   30.14447      0           0.17 Slide-0989632_AKU040 T12-17 MM-LQ1
## Slide-0989648_KAIC004 T1-16 | LE-LQ1 | PanCK positive        59   old  335.43742      2           0.29 Slide-0989648_KAIC004 T1-16 LE-LQ1
## Slide-0991905_KAIC017 T1-12 | LE-UQ2 | PanCK positive        45 young  737.25681      5           0.34 Slide-0991905_KAIC017 T1-12 LE-UQ2
## Slide-0991913_KAIC022 T2-12 | LE-LQ1 | PanCK positive        68   old  724.91041     11           0.26 Slide-0991913_KAIC022 T2-12 LE-LQ1
## Slide-0991913_KAIC022 T2-12 | LE-UQ2 | PanCK positive        68   old  192.09084      0           0.26 Slide-0991913_KAIC022 T2-12 LE-UQ2
## Slide-0991927_KAIC034 T1-10 | LE-LQ2 | PanCK positive        49 young 1208.34436     11           0.28 Slide-0991927_KAIC034 T1-10 LE-LQ2
## Slide-0991931_KAIC007 T1-5 | MM-LQ1 | PanCK positive         50   old 2294.34706      6           0.24  Slide-0991931_KAIC007 T1-5 MM-LQ1
##                                                    LOT_Human_Immune_Cell_Profiling_Protein_Core LOT_Human_IO_Drug_Target_Protein_Module
## Slide-0989632_AKU040 T12-17 | LE-LQ1 | PanCK positive                                    ICPH10002                               DRGH10001
## Slide-0989632_AKU040 T12-17 | MM-LQ1 | PanCK positive                                    ICPH10002                               DRGH10001
## Slide-0989648_KAIC004 T1-16 | LE-LQ1 | PanCK positive                                    ICPH10002                               DRGH10001
## Slide-0991905_KAIC017 T1-12 | LE-UQ2 | PanCK positive                                    ICPH10002                               DRGH10001
## Slide-0991913_KAIC022 T2-12 | LE-LQ1 | PanCK positive                                    ICPH10002                               DRGH10001
## Slide-0991913_KAIC022 T2-12 | LE-UQ2 | PanCK positive                                    ICPH10002                               DRGH10001
## Slide-0991927_KAIC034 T1-10 | LE-LQ2 | PanCK positive                                    ICPH10002                               DRGH10001
## Slide-0991931_KAIC007 T1-5 | MM-LQ1 | PanCK positive                                     ICPH10002                               DRGH10001
##                                                    LOT_Human_Immune_Cell_Typing_Protein_Module LOT_Human_Pan_Tumor_Protein_Module      TumorEdge
## Slide-0989632_AKU040 T12-17 | LE-LQ1 | PanCK positive                                   ICTH10002                             474038 Non tumor edge
## Slide-0989632_AKU040 T12-17 | MM-LQ1 | PanCK positive                                   ICTH10002                             474038 Non tumor edge
## Slide-0989648_KAIC004 T1-16 | LE-LQ1 | PanCK positive                                   ICTH10002                             474038 Non tumor edge
## Slide-0991905_KAIC017 T1-12 | LE-UQ2 | PanCK positive                                   ICTH10002                             474038 Non tumor edge
## Slide-0991913_KAIC022 T2-12 | LE-LQ1 | PanCK positive                                   ICTH10002                             474038 Non tumor edge
## Slide-0991913_KAIC022 T2-12 | LE-UQ2 | PanCK positive                                   ICTH10002                             474038 Non tumor edge
## Slide-0991927_KAIC034 T1-10 | LE-LQ2 | PanCK positive                                   ICTH10002                             474038 Non tumor edge
## Slide-0991931_KAIC007 T1-5 | MM-LQ1 | PanCK positive                                    ICTH10002                             474038 Non tumor edge
##                                                    low_nc
## Slide-0989632_AKU040 T12-17 | LE-LQ1 | PanCK positive      1
## Slide-0989632_AKU040 T12-17 | MM-LQ1 | PanCK positive      1
## Slide-0989648_KAIC004 T1-16 | LE-LQ1 | PanCK positive      1
## Slide-0991905_KAIC017 T1-12 | LE-UQ2 | PanCK positive      1
## Slide-0991913_KAIC022 T2-12 | LE-LQ1 | PanCK positive      1
## Slide-0991913_KAIC022 T2-12 | LE-UQ2 | PanCK positive      1
## Slide-0991927_KAIC034 T1-10 | LE-LQ2 | PanCK positive      1
## Slide-0991931_KAIC007 T1-5 | MM-LQ1 | PanCK positive       1

before <- nrow(dfs1)
dfs1 <- dfs1[dfs3$Nuclei>min_nc,]
dfs3 <- dfs3[dfs3$Nuclei>min_nc,]
after <- nrow(dfs1)
message(paste0("Removed ", before - after , " sample(s) with low nuclei counts or low area."))

#########################################################################
# 2.4 Quality control of housekeeper proteins and background (isotypes) #
#########################################################################

# Before we do normalization, we need to check if the control molecules are themselves correlated with the predictors of interests.
# There are two classes of control molecules in the dataset: isotypes and housekeepers.

# Analyst notes: Housekeeper proteins include S6, Histone H3, and GAPDH. Background isotypes include Ms IgG2a, Ms IgG1, and Rb IgG.
# We see high correlation between the IgGs as well as the housekeeper proteins. Ultimately, Histone H3 and S6 display the highest
# correlation and least amount of noise, so we will use these two housekeeper proteins to perform normalization on the data.

# biological negatives (isotype controls)
iggs = colnames(dfs1)[grepl("IgG", colnames(dfs1))]
message(paste0("All IgGs: ",paste0(iggs, collapse=", ")))
# housekeepers:
hks = c("S6", "Histone H3", "GAPDH")
message(paste0("All HKGs: ",paste0(hks, collapse=", ")))

##########################
# 2.4.1 Isotype Controls #
##########################

# Plot the joint expression of each pairwise igg control molecule to see if:

# There’s evidence that one igg is noisy compared to others
# There’s clustering by predictors

##
# All Samples, CK, CancerType, and EpithelialStatus
##
iggs_df <- log2(as.data.frame(dfs1[,iggs]))
n <- ncol(iggs_df)
#iggs_df <- cbind(iggs_df, dfs3 %>% dplyr::select(PanCK, TumorEdge, Morphology))
iggs_df <- cbind(iggs_df, dfs3 %>% dplyr::select(PanCK, CancerType, EpithelialStatus))
p0 <- GGally::ggpairs(data=iggs_df, columns=1:n, progress=FALSE,
            ggplot2::aes(alpha=0.1)) + theme_bw()
#ggsave(p0, file=file.path(fig_dir, "QC_pairs_background.svg"), width=n*2.5, height=n*2.5)
ggsave(p0, file=file.path(fig_dir, "QC_pairs_background.png"), width=n*2.5, height=n*2.5)
pairs_custom <- function(dat, column_of_interest){
  require(rlang)
  n <- ncol(dat)
  p <- dat %>% ggpairs(., 
    mapping = ggplot2::aes(colour=as.character({{column_of_interest}}), alpha=0.1),
    columns=1:n, progress=FALSE,
    lower = list(continuous = wrap("smooth", alpha=0.3, size=0.3),
                 combo=wrap("facethist", bins=30))
  ) + theme_bw() 
  return(p)
}
p1 <- pairs_custom(iggs_df %>% dplyr::select(eval(iggs), PanCK), PanCK)
#ggsave(p1, file=file.path(fig_dir, "QC_pairs_background_CK.svg"), width=n*2.8, height=n*2.8)
ggsave(p1, file=file.path(fig_dir, "QC_pairs_background_PanCK.png"), width=n*2.8, height=n*2.8)

#p2 <- pairs_custom(iggs_df %>% dplyr::select(eval(iggs), TumorEdge), TumorEdge)
p2 <- pairs_custom(iggs_df %>% dplyr::select(eval(iggs), CancerType), CancerType)
#ggsave(p2, file=file.path(fig_dir, "QC_pairs_background_TumorEdge.svg"), width=n*2.8, height=n*2.8)
ggsave(p2, file=file.path(fig_dir, "QC_pairs_background_CancerType.png"), width=n*2.8, height=n*2.8)

#p3 <- pairs_custom(iggs_df %>% dplyr::select(eval(iggs), Morphology), Morphology)
p3 <- pairs_custom(iggs_df %>% dplyr::select(eval(iggs), EpithelialStatus), EpithelialStatus)
#ggsave(p3, file=file.path(fig_dir, "QC_pairs_background_Morphology.svg"), width=n*2.8, height=n*2.8)
ggsave(p3, file=file.path(fig_dir, "QC_pairs_background_EpithelialStatus.png"), width=n*2.8, height=n*2.8)

###############################
# 2.4.2 Housekeeping Controls # 
###############################

# Run the same analysis for the housekeepers. This is the same overall process as for the isotypes.

hks_df <- log2(as.data.frame(dfs1[,hks]))
n <- ncol(hks_df)
#hks_df <- cbind(hks_df, dfs3 %>% dplyr::select(PanCK, TumorEdge, Morphology))
hks_df <- cbind(hks_df, dfs3 %>% dplyr::select(PanCK, CancerType, EpithelialStatus))
p4 <- ggpairs(data=hks_df, columns=1:n, progress=FALSE,
            ggplot2::aes(alpha=0.1)) + theme_bw()
#ggsave(p4,file=file.path(fig_dir, "QC_pairs_housekeepers.svg"), width=n*2.5, height=n*2.5)
ggsave(p4,file=file.path(fig_dir, "QC_pairs_housekeepers.png"), width=n*2.5, height=n*2.5)

p5 <- pairs_custom(hks_df %>% dplyr::select(eval(hks), PanCK), PanCK)
#ggsave(p5, file=file.path(fig_dir, "QC_pairs_housekeepers_CK.svg"), width=n*2.8, height=n*2.8)
ggsave(p5, file=file.path(fig_dir, "QC_pairs_housekeepers_PanCK.png"), width=n*2.8, height=n*2.8)

#p6 <- pairs_custom(hks_df %>% dplyr::select(eval(hks), TumorEdge), TumorEdge)
p6 <- pairs_custom(hks_df %>% dplyr::select(eval(hks), CancerType), CancerType)
#ggsave(p6, file=file.path(fig_dir, "QC_pairs_housekeepers_TumorEdge.svg"), width=n*2.8, height=n*2.8)
ggsave(p6, file=file.path(fig_dir, "QC_pairs_housekeepers_CancerType.png"), width=n*2.8, height=n*2.8)

#p7 <- pairs_custom(hks_df %>% dplyr::select(eval(hks), Morphology), Morphology)
p7 <- pairs_custom(hks_df %>% dplyr::select(eval(hks), EpithelialStatus), EpithelialStatus)
#ggsave(p7, file=file.path(fig_dir, "QC_pairs_housekeepers_Morphology.svg"), width=n*2.8, height=n*2.8)
ggsave(p7, file=file.path(fig_dir, "QC_pairs_housekeepers_EpithelialStatus.png"), width=n*2.8, height=n*2.8)

#####################################
# 2.4.3 Pairwise plot with geomeans #
#####################################
iggs_df <- (as.data.frame(dfs1[,iggs]))
hks_df <- (as.data.frame(dfs1[,hks]))
igg_geo_means <- as.numeric(apply(iggs_df, 1, EnvStats::geoMean))
hks_geo_means <- as.numeric(apply(hks_df, 1, EnvStats::geoMean))
geo_means <- as.data.frame(cbind(igg_geo_means, hks_geo_means))
geos_df <- cbind(geo_means, dfs3 %>% dplyr::select(Area, Nuclei))
n <- ncol(geos_df)
p8 <- GGally::ggpairs(data=geos_df, columns=1:4, progress=FALSE,
            ggplot2::aes(alpha=0.1)) + theme_bw()
#ggsave(p8, file=file.path(fig_dir, "QC_pairwise_geomeans_area_nuclei.svg"), width=n*2.8, height=n*2.8)
ggsave(p8, file=file.path(fig_dir, "QC_pairwise_geomeans_area_nuclei.png"), width=n*2.8, height=n*2.8)
p8

#######################
# 2.5 Signal to Noise #
#######################
backgrounds <- iggs

snr <- get_snr(exp=dfs1, backgrounds = backgrounds, n_processors = 5, 
               x=sd_multiplier, type="SBR")
snr2 <- reshape2::melt(snr, id.var=NA)
colnames(snr2) <- c("name", "protein", "snr")
#snr2 <- base::merge(snr2, dfs3 %>% mutate(name=rownames(dfs3)) %>% dplyr::select(name, PanCK, TumorEdge, Morphology), by="name")
snr2 <- base::merge(snr2, dfs3 %>% mutate(name=rownames(dfs3)) %>% dplyr::select(name, PanCK, CancerType, EpithelialStatus), by="name")
                            
# Order proteins based on median SNR. 
medians <- ddply(snr2, .(protein), summarize, Median=median(snr), 
            UQ=as.numeric(quantile(snr, probs=0.75))) %>% arrange(UQ)
snr2$protein <- factor(as.character(snr2$protein), 
              levels=as.character(medians$protein), order=TRUE)
snr2 <- snr2 %>% arrange(protein, name)

# Now we can plot the SNR for all probes. Within a given protein, it is valuable to assess the distribution of ROIs to determine if certain ROIs are above background while others are not; such differences can be the result of differing ROI type, tissue ID, or some other biological variable of interest. We will display these data colored by PanCK, TumorEdge, and Morphology, as well.
# Analyst Notes: The signal to noise ratio is above background for most proteins. We do not filter proteins based on SNR, but will use these results as context while investigating the data further.

###############
# All Samples #
###############
a <- ifelse(levels(snr2$protein) %in% iggs, "grey", "grey20")
p <- ggplot(data=snr2, aes(x=protein, y=log2(snr))) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(alpha=0.5, size=0.3) +
  geom_hline(yintercept = 0, lty="dashed") +
  theme_bw() + xlab("Protein") + ylab("log2(SNR)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1, colour=a),
        legend.background = element_rect(fill="white",
                                  size=0.5, linetype="solid", 
                                  colour ="grey20"),
        legend.position = c(0.05, 0.75)) +
  guides(colour=guide_legend(override.aes = list(size=2)))
ggsave(p, file= file.path(fig_dir, "snr_plot.png"), height=8, width=16)
p

#########
# PanCK #
#########
a <- ifelse(levels(snr2$protein) %in% iggs, "grey", "grey20")
p <- ggplot(data=snr2, aes(x=protein, y=log2(snr))) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(alpha=0.5, size=0.3) +
  geom_point(position=position_jitterdodge(0.1), aes(colour=PanCK), alpha=0.5, size=0.3) +
  geom_hline(yintercept = 0, lty="dashed") +
  #scale_colour_manual(values = c("#999999", "#FABFD2", "#318026" )) +
  scale_colour_manual(values = pal_main[["PanCK"]]) +
  theme_bw() + xlab("Protein") + ylab("log2(SNR)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1, colour=a),
        legend.background = element_rect(fill="white",
                                  size=0.5, linetype="solid", 
                                  colour ="grey20"),
        legend.position = c(0.05, 0.75)) +
  guides(colour=guide_legend(override.aes = list(size=2)))
ggsave(p, file= file.path(fig_dir, "snr_plot_PanCK.png"), height=8, width=16)
p

##############
# CancerType #
##############
a <- ifelse(levels(snr2$protein) %in% iggs, "grey", "grey20")
p <- ggplot(data=snr2, aes(x=protein, y=log2(snr))) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(alpha=0.5, size=0.3) +
  geom_point(position=position_jitterdodge(0.1), aes(colour=CancerType), alpha=0.8, size=0.3) +
  geom_hline(yintercept = 0, lty="dashed") +
  #scale_colour_manual(values =c("#3A6CA1", "#FFD861" )) +
 # scale_colour_manual(values =c("darkorange", "dodgerblue", "gold1", "lightgreen", "lightpink" )) +
  scale_colour_manual(values = pal_main[["CancerType"]]) +
  theme_bw() + xlab("Protein") + ylab("log2(SNR)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1, colour=a),
        legend.background = element_rect(fill="white",
                                  size=0.5, linetype="solid", 
                                  colour ="grey20"),
        legend.position = c(0.04, 0.75)) +
  guides(colour=guide_legend(override.aes = list(size=2)))
ggsave(p, file= file.path(fig_dir, "snr_plot_CancerType.png"), height=8, width=16)
p

##############
# EpithelialStatus #
##############
a <- ifelse(levels(snr2$protein) %in% iggs, "grey", "grey20")
p <- ggplot(data=snr2, aes(x=protein, y=log2(snr))) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(alpha=0.5, size=0.3) +
  geom_point(position=position_jitterdodge(0.1), aes(colour=EpithelialStatus), alpha=0.5, size=0.3) +
  geom_hline(yintercept = 0, lty="dashed") +
  #scale_colour_manual(values =c("#A66293", "#F28E2B", "#788D60")) +
  scale_colour_manual(values = pal_main[["EpithelialStatus"]]) +
  theme_bw() + xlab("Protein") + ylab("log2(SNR)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1, colour=a),
        legend.background = element_rect(fill="white",
                                  size=0.5, linetype="solid", 
                                  colour ="grey20"),
        legend.position = c(0.04, 0.75)) +
  guides(colour=guide_legend(override.aes = list(size=2)))
ggsave(p, file= file.path(fig_dir, "snr_plot_EpithelialStatus.png"), height=8, width=16)
p

#message("=== HERE ==")

#####################
# 2.6 Normalization #
#####################
# The purpose of normalization is to adjust for technical variables, such as ROI/AOI surface area and tissue quality, and enable meaningful biological and statistical discoveries. We normalize based on the two housekeeper proteins and check to see if there are AOIs that have very large normalization factors.

# Analyst notes: There are 3 AOI segments with a |log2 norm factor| of 5 or more. We removed 3 segments at this step. After normalization, we remove the housekeeper and IgG proteins.

dfs <- list(dfs1, NULL, dfs3)

##
## setdiff(hks, "GAPDH"): all elements of hks except "GAPDH"
##
## geo_means <- as.numeric(apply(dfs1[rownames(dfs[[1]]),iggs], 1, EnvStats::geoMean))
hks_geo_means <- as.numeric(apply(dfs1[rownames(dfs[[1]]),setdiff(hks, "GAPDH")], 1, EnvStats::geoMean))

norm_factors <- mean(hks_geo_means) / hks_geo_means
dfs[[2]] <- sweep(dfs[[1]], 1, norm_factors, "*")

###
# Compare histograms of normalization factors (log2 transformed).
###
dfs[[3]]$norm_factor <- norm_factors
p<-ggplot(dfs[[3]], aes(x=log2(norm_factor))) +
  geom_histogram(bins=20) + 
  geom_vline(xintercept = c(-nf_threshold, nf_threshold)) + 
  theme_bw() + 
  ylab("Number of AOIs") +
  scale_fill_manual(values =c("#3A6CA1", "#FFD861", "#CF4244", "#47BAB4"))
#ggsave(p, file= file.path(fig_dir, "normalization_factors.svg"), height=8, width=8)
ggsave(p, file= file.path(fig_dir, "normalization_factors.png"), height=8, width=8)
p

#################
# print message #
#################
message(paste0(
  "There are ", 
  nrow(filter(dfs[[3]], abs(log2(norm_factor))>nf_threshold)), 
  " (", 
  round(100*nrow(filter(dfs[[3]], abs(log2(norm_factor))>nf_threshold))/nrow(dfs[[3]]), 3), 
  "%) of AOIs with a |log2 norm factor| of ", 
  nf_threshold, 
  " or more."
))

#####################
# filtering samples #
#####################
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
to_keep <- row.names(filter(dfs[[3]], abs(log2(norm_factor))<=nf_threshold))
to_remove <- row.names(filter(dfs[[3]], abs(log2(norm_factor))>nf_threshold))
to_remove
## 
#to_remove
dfs <- filter_samples(DFS=dfs, samples_to_keep = to_keep)

############################################
# 2.6.1 Removing isotypes and housekeepers #
############################################

message("=== removed hks and igs to create the final dataset for downstream analysis. ===")

# Finally, we remove the background probes (IgGs) and housekeeper proteins, and create the final dataset for downstream analysis.

filter_features <- function(DFS, to_remove, proteins){
  if(nrow(DFS[[1]]) != nrow(DFS[[2]]) | nrow(DFS[[1]]) != nrow(DFS[[3]])){
    stop("Check dimensions. Error 1.")
  }
  before <- ncol(DFS[[1]])
  DFS[[1]] <- DFS[[1]][,-which(colnames(DFS[[1]]) %in% to_remove)]
  DFS[[2]] <- DFS[[2]][,-which(colnames(DFS[[2]]) %in% to_remove)]
  after <- ncol(DFS[[1]])
  message(paste0("Removed ", before - after, proteins, " features."))
  if(nrow(DFS[[1]]) != nrow(DFS[[2]]) | nrow(DFS[[1]]) != nrow(DFS[[3]])){
    stop("Check dimensions. Error 2.")
  }
  return(DFS)
}

# Remove the housekeeping proteins and iggs
dfs <- filter_features(dfs, to_remove=hks, " housekeeper proteins ")
dfs <- filter_features(dfs, to_remove=iggs, " IgGs ")

#####################################
# Lastly, we save the data to disk. #
#####################################
#Original dataset
openxlsx::write.xlsx(dfs[[1]], file.path(outdata_dir, "Protein_GeoMx_Raw_Data.xlsx"), rowNames=TRUE)
openxlsx::write.xlsx(dfs[[2]], file.path(outdata_dir, "Protein_GeoMx_Normalized_and_Filtered_Data.xlsx"), rowNames=TRUE)
openxlsx::write.xlsx(dfs[[3]], file.path(outdata_dir, "Protein_GeoMx_Normalized_and_Filtered_Annotations.xlsx"), rowNames=TRUE)

#Save R object
saveRDS(dfs, file.path(object_dir, "Geomx_Data_list.rds"))

## save image
save.image(file=file.path(object_dir,"step2.qc.norm.RData"))

##
sessionInfo()
