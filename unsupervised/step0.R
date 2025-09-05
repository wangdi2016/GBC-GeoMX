setwd("/Users/wangdi/work/DWX/geomx/jill/2025-03-26/GBC")

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
datadir <- "/Users/wangdi/work/DWX/geomx/jill/2025-03-26/GBC/Data"

#raw <- openxlsx::read.xlsx(file.path(datadir, "Initial Dataset.xlsx"), sheet=2)
#raw <- openxlsx::read.xlsx(file.path(datadir, "Initial Dataset.xlsx"), sheet="Modified")
#raw <- openxlsx::read.xlsx(file.path(datadir, "InitialDataset.Hela-DW2.xlsx"), sheet="Sheet3-707Seg")
raw <- openxlsx::read.xlsx(file.path(datadir, "Initial_Dataset.20250326.dw.xlsx"), sheet="Sheet1")
rownames(raw) <- raw$Custom.Segment.Name

#annotations <- openxlsx::read.xlsx(file.path(datadir, "Annotation_shared_091923.xlsx"), sheet=2)
#annotations <- openxlsx::read.xlsx(file.path(datadir, "Annotation_shared_091923.xlsx"), sheet="Sheet2")
#annotations <- openxlsx::read.xlsx(file.path(datadir, "Annotation_08182024.Hela-DW2.xlsx"), sheet="Sheet2CleanDW")
annotations <- openxlsx::read.xlsx(file.path(datadir, "Annotation_file.20250326.dw.xlsx"), sheet="Sheet1")
annotations$Custom.Segment.Name <- paste0(annotations$Scan.Name, " | ", annotations$`ROI.(Label)`, " | ", annotations$`Segment(Name/Label)`)
rownames(annotations) <- annotations$Custom.Segment.Name

#Need to merge annotations with initial dataset
raw_annotated <- base::merge(annotations, raw, by=0)
rownames(raw_annotated) <- raw_annotated$Row.names
raw_annotated$Row.names <- NULL

write.table(raw_annotated, file="raw_annotated.csv", sep=",")
raw_subset <- raw_annotated[,1:34]

