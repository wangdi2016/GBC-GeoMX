###
#source("./init.R")
source("./init.dw.R")

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
  fileType1 = "svg"
  fileType2 = "pdf"
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
  out_file_name1 <- paste0(out_dir, "/", prefix, '_volcano', '.', fileType1)
  out_file_name2 <- paste0(out_dir, "/", prefix, '_volcano', '.', fileType2)
  ggsave(plt, filename = out_file_name1, width=WIDTH, height=HEIGHT)
  ggsave(plt, filename = out_file_name2, width=WIDTH, height=HEIGHT)
  
  out_file_name2 <- paste0(out_dir, "/", prefix, '.tsv')
  write.table(dfx, file = out_file_name2, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
  
  return(list(plt, dfx))
  
} 
