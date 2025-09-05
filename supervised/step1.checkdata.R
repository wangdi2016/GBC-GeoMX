###################################################
# GeoMx® Custom Report: NCI                       #
#                                                 #
# Difei Wang CGR/DCEG/NCI  2024-09-23             #
#                                                 #
# modified the following code                     #
# Felicia New, PhD, NanoString Technologies, Inc. #
# 2023-10-17                                      #
###################################################

# 1 Introduction
# We present here a data analysis report for GeoMx® Protein with nCounter readout data. Samples are from breast cancer patients.

#####################################
# 1.1 Slides and Profiling Strategy #
#####################################

# The dataset of interest includes 31 tumor samples, with one tissue per slide. Regions of Interest (ROIs) were placed to study the heterogeneity of the tumor microenvironment. ROIs were further segmented into CK+ and CK- compartments, known as Areas of Illumination (AOIs).

##############################
# 1.2 Experimental Questions #
##############################

## DW CK --> PanCK
## DW Proximal PanCK --> ROI neighborhood PanCK high vs. low
## DW Proximal CD45  --> ROI neighborhood CD45  high vs. low
## DW Luminal --> PAM50 (AIMS)

# This analysis report focuses on these key experimental comparisons (all comparisons are in the CK- compartment unless otherwise noted):

# Question 1: CK+ vs CK- AOIs
# Question 2: Predominantly tumor vs combined stroma
# Question 3: Tumor edge vs non-tumor edge
# Question 4: Distance to the centroid: high vs low
# Question 5: Proximal vs distal PanCK
# Question 6: Proximal vs distal CD45
# Question 7: Luminal vs non-luminal
# Question 8: Low vs. high grade
# Question 9: Parity status: high vs low
# Question 10: BMI: normal weight vs obese
# Question 11: Waist to hip ratio: high vs low
# Question 12: Younger vs older women

# For each comparison above, we look for differentially expressed proteins with various data visualizations - heatmaps, volcano plots, and violin plots of the top differentially expressed protein(s).

###################
# 1.3 Import Data #
###################

##
## load package, set parameters and load custom functions
##
source("./init.dw.R")

## read data table
#datadir <- "/Users/wangdi/work/DWX/089_sDAS_NCI_Lawrence_Report_DW2/data"
#datadir <- "/Users/wangdi/work/DWX/geomx/089_sDAS_NCI_Lawrence_Report_DW2/data"
#datadir <- "/Users/wangdi/work/DWX/geomx/2024-08-05/Rose-GeoMX/089_sDAS_NCI_Lawrence_Report_DW2/Data"
#datadir <- "/Users/wangdi/work/DWX/geomx/2024-09-23/089_sDAS_NCI_Lawrence_Report_DW3/Data"
#datadir <- "/Users/wangdi/work/DWX/geomx/2024-09-23/089_sDAS_NCI_Lawrence_Report_DW5-ALL-AOIs/Data"
datadir <- "/Users/wangdi/work/DWX/geomx/jill/2025-08-12/GBC/Data"

xlsxfile1 <- "Initial_Dataset_wholeslide_all_samples.dw.2025-08-12.xlsx"
xlsxfile2 <- "Annotation_file_wholeslide_all_sample.dw.2025-08-12.xlsx"

#raw <- openxlsx::read.xlsx(file.path(datadir, "Initial Dataset.xlsx"), sheet=2)
#raw <- openxlsx::read.xlsx(file.path(datadir, "Initial Dataset.xlsx"), sheet="Modified")
#raw <- openxlsx::read.xlsx(file.path(datadir, "InitialDataset.Hela-DW2.xlsx"), sheet="Sheet3-707Seg")
raw <- openxlsx::read.xlsx(file.path(datadir, xlsxfile1), sheet="Sheet1")
rownames(raw) <- raw$Custom.Segment.Name

#annotations <- openxlsx::read.xlsx(file.path(datadir, "Annotation_shared_091923.xlsx"), sheet=2)
#annotations <- openxlsx::read.xlsx(file.path(datadir, "Annotation_shared_091923.xlsx"), sheet="Sheet2")
#annotations <- openxlsx::read.xlsx(file.path(datadir, "Annotation_08182024.Hela-DW2.xlsx"), sheet="Sheet2CleanDW")
annotations <- openxlsx::read.xlsx(file.path(datadir, xlsxfile2), sheet="Sheet1")
annotations$Custom.Segment.Name <- paste0(annotations$Scan.Name, " | ", annotations$`ROI.(Label)`, " | ", annotations$`Segment(Name/Label)`)
rownames(annotations) <- annotations$Custom.Segment.Name

#Need to merge annotations with initial dataset
raw_annotated <- base::merge(annotations, raw, by=0)
rownames(raw_annotated) <- raw_annotated$Row.names
raw_annotated$Row.names <- NULL

#########################################################################
# Next, split the raw data into (raw) expression and GeoMx annotations. #
#########################################################################

#Specify the columns that are annotation vs protein counts. This can vary between projects.
#geomx_annots <- raw_annotated %>% dplyr::select(1:48)
geomx_annots <- raw_annotated %>% dplyr::select(1:34)   # from Scan.Name.x to Target.name.(display.name)
#raw <- raw_annotated %>% dplyr::select(49:100)
raw <- raw_annotated %>% dplyr::select(35:94)           # from CD20 to CD40

#Fix protein names
raw <- plyr::rename(raw, c("Ms.IgG1" = "Ms IgG1", "Ms.IgG2a"="Ms IgG2a", "Histone.H3"="Histone H3", "Rb.IgG" ="Rb IgG" ))
#change to numeric
for(i in c(1:ncol(raw))) {
    raw[,i] <- as.numeric(raw[,i])
}
raw <- data.matrix(raw)
dim(raw)

#Fix annotations
#geomx_annots <- geomx_annots[,c(3:21,32,33,36,42:47)] ## HKBC
geomx_annots <- geomx_annots[,c(3:6,17:18,21,27:33)] ## GBC JILL
geomx_annots <-  plyr::rename(geomx_annots, c("Segment(Name/Label)"="Segment", "AOI.surface.area"="Area", "AOI.nuclei.count"="Nuclei"))

#Modify Tumor edge annotation (LE vs the rest)
#geomx_annots$TumorEdge <- ifelse(geomx_annots$TumorEdge_LE_vs_rest == "LE", "Tumor edge", "Non tumor edge")
# TumorEdge

###########################
# 1.4 Overview of samples #
###########################

# We can summarize the study design in a Sankey diagram. This visualization is useful for showing the flow of samples and how they relate to biological annotations. The width of a cord in the figure represents how many segments are in common between the two annotations they connect.

annots <- geomx_annots

#factors_of_interest <- c("CK","Morphology","TumorEdge","DistanceToCentroid","PanCK_proximity","CD45_proximity","Luminal","Grade","Parity","BMI","Waist_to_hip_ratio","Age")
#factors_of_interest <- c("CK","Morphology","TumorEdge","PanCK_proximity","CD45_proximity","Luminal","Grade","Parity","BMI","Waist_to_hip_ratio","Age")
# factors_of_interest <- c("PanCK","Morphology","TumorEdge", "ROI_neighborhood_PanCK","ROI_neighborhood_CD45","PAM50", "Grade","Parity","BMI","Waist_to_hip_ratio","Age")  ## HKBC
factors_of_interest <- c("PanCK","CancerType","EpithelialStatus")  ## GBC JILL

pal_main <-
 getColorPalette(annots[, factors_of_interest], method="Map")
#Static plot
p <- easyalluvial::alluvial_wide(annots[factors_of_interest])
p <- p + theme(axis.text.x = element_text(face="plain", size=13, angle = 45, hjust = 1, vjust = 1))
ggsave(p, filename = file.path(fig_dir, "Main_Sankey.png") , width=16, height=8)
ggsave(p, filename = file.path(fig_dir, "Main_Sankey.svg") , width=16, height=8)
ggsave(p, filename = file.path(fig_dir, "Main_Sankey.pdf") , width=16, height=8)
#Interactive plot
pc <- parcats::parcats(p,
        marginal_histograms = FALSE,
        data_input = annots,
        labelfont = list(size = 12, color = "black"),
        hoveron="dimension",
        hoverinfo = "count",
        offset_marginal_histograms = 0.7)
#pc
#p

save.image(file=file.path(object_dir,"step1.checkdata.RData"))

##
sessionInfo()

