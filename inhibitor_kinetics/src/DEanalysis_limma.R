### INHIBITOR KINETICS ###
# description:  differential expression analysis using MS intensities
# input:        output matrices from _parse-qiSPI-intensities.R
# output:       volcano plots, b5 specific peptides
# author:       HR, modified from YH

library(biobroom)
library(dplyr)
library(eulerr)
library(gprofiler2)
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(limma)
library(naniar)
library(stringr)
library(tidyr)
library(readr)
library(qvalue)

# source("src/data_utils.R")
source("../brainstorming/src/invitroSPI_utils.R")


### HYPERPARAMETERS ###
# Differential expression parameters
pval.threshold = 0.1
# fc.threshold = sqrt(2)
fc.threshold = log(2)

### INPUT ###
load("data/fold-changes.RData")
load("data/fold-changes_medpolish.RData")
load("data/fold-changes_double-log.RData")
load("data/intensities-4hrs-per-substrate.RData")

INT_4hrs = read.csv("data/intensity-table-4hrs.csv", stringsAsFactors = F)
intIDX = c(5:20)


### MAIN PART ###
# ----- data normalisation & prefiltering -----
# done in _parse-qi-SPI-intensities.R
# QUANT = FC %>% as.data.frame()
# QUANT = FCm %>% as.data.frame()
QUANT = FCl %>% as.data.frame()
# QUANT = FC_perSubs %>% as.data.frame()
# QUANT = FCsq %>% as.data.frame()


# ----- specify comparisons -----
names(QUANT)

# select peptides that are downregulated in the following comparisons
no_idx = c(1:4)
b1_idx = c(13:16)
b2_idx = c(9:12)
b5_idx = c(5:8)
all_idx = c(1:16)

comp = list(b5 = list(target = b5_idx,
                      b5_vs_all = c(b5_idx, all_idx[!all_idx %in% b5_idx]),
                      b5_vs_no = c(b5_idx, no_idx),
                      b5_vs_b1 = c(b5_idx, b1_idx),
                      b5_vs_b2 = c(b5_idx, b2_idx)),
            b1 = list(target = b1_idx,
                      b1_vs_all = c(b1_idx, all_idx[!all_idx %in% b1_idx]),
                      b1_vs_no = c(b1_idx, no_idx),
                      b1_vs_b5 = c(b1_idx, b5_idx),
                      b1_vs_b2 = c(b1_idx, b2_idx)),
            b2 = list(target = b2_idx,
                      b2_vs_all = c(b2_idx, all_idx[!all_idx %in% b2_idx]),
                      b2_vs_no = c(b2_idx, no_idx),
                      b2_vs_b5 = c(b2_idx, b5_idx),
                      b2_vs_b1 = c(b2_idx, b1_idx)))

suppressWarnings(sapply(paste0("results/",names(comp),"/"), dir.create))

################################################################################
# subunit = "b5"
# comparison = "b5_vs_no"

# function that does limma DE analysis for given subunit/comparison
DEanalysis = function(subunit, comparison) {
  
  suppressWarnings(dir.create(paste0("results/",subunit,"/",comparison,"/")))
  
  # ----- limma -----
  ### create design matrix ###
  # modify intensity table
  idx = comp[[subunit]][[comparison]]
  L = QUANT[, idx]
  names(L) = c(rep(paste0(subunit, "_inactive"), length(comp[[subunit]][["target"]])),
               rep(paste0(subunit, "_active"), length(idx)-length(comp[[subunit]][["target"]])))
  
  sample_info = data.frame(condition = names(L))
  design = model.matrix(~condition, sample_info)
  
  
  ### fit linear model ###
  data_fit = lmFit(object = L,
                   design = design,
                   method = "ls")
  
  eb = eBayes(data_fit)
  topTable(eb, adjust.method = "BH")
  res <- tidy(eb)
  
  ## check statistics
  qobj <- qvalue(p = res$p.value)
  plot(qobj)
  hist(qobj)
  
  res$padj <- p.adjust(res$p.value, "BH")
  res$qvalues <- qobj$qvalues
  
  # Plot p-value distribution per test
  pval <- ggplot(res, aes(p.value)) + 
    geom_histogram(bins = 100) +
    facet_wrap(~ term) +
    theme_bw() +
    ggtitle("p-values") 
  padj <- ggplot(res, aes(padj)) + 
    geom_histogram(bins = 100) +
    facet_wrap(~ term) +
    theme_bw() +
    ggtitle("Adjusted p-values") 
  
  gridExtra::grid.arrange(pval,padj, nrow=2)
  
  res <- res %>%
    mutate(threshold = ifelse(estimate > fc.threshold, "up", ifelse(estimate < -fc.threshold , "down", "small"))) %>%
    mutate(threshold = replace(threshold, padj >= pval.threshold,"not significant"))
  dim(res)
  
  # Export DE stats
  DE_stats <- res %>%
    dplyr::select(c("threshold","term")) %>%
    table() %>%
    as.data.frame.matrix() %>%
    t()
  DE_stats <- cbind(Contrast = rownames(DE_stats),
                    DE_stats) %>%
    as.data.frame()
  DE_stats
  
  
  # ----- volcano plots -----
  spliceTypes = factor(INT_4hrs$spliceType,
                       levels = c("cis", "revCis", "trans", "PCP", "type_multi-mapper"))
  
  Volcanoplots <- res %>%
    ggplot(aes(y=-log10(padj), x=estimate, color=spliceTypes)) + 
    geom_point(stat="identity", position="identity", alpha=0.5, size=1) + 
    theme_bw() + 
    scale_color_manual("type", values = c(plottingCols[["cis"]], plottingCols[["revCis"]],
                                          plottingCols[["trans"]], plottingCols[["PCP"]], "gray")) +
    theme(text=element_text(family="sans", face="plain", color="#000000", size=12, hjust=0.5, vjust=0.5)) + 
    xlab("log2 fold change") + 
    ylab("-log10 p.adj") +
    ggtitle(comparison,
            sub = paste0("differential peptide intensities - down: ",DE_stats$down,", up: ", DE_stats$up)) +
    geom_hline(aes(yintercept = -log10(pval.threshold)), colour = "black") +
    geom_vline(xintercept = c(-fc.threshold,fc.threshold), colour = "black", linetype = "dashed") +
    theme_classic()
  
  Volcanoplots
  
  
  # ggsave(filename = "~/Documents/Studium/Fachvertiefung+BA/Praktikumsbericht/plots/volcanoPlot.ps",
  #        plot = Volcanoplots, device = cairo_ps,
  #        dpi = "retina", height = 3*5, width = 4*5, units = "cm")
  # ggsave(filename = "~/Documents/Studium/Fachvertiefung+BA/Praktikumsbericht/plots/volcanoPlot.png",
  #        plot = Volcanoplots, device = "png",
  #        dpi = "retina", height = 3*5, width = 4*5, units = "cm")
  # 
  # ggsave(filename = "~/Documents/Studium/Fachvertiefung+BA/Praktikumsbericht/plots/qVals.ps",
  #        plot = plot(qobj), device = cairo_ps,
  #        dpi = "retina", height = 3*5, width = 4*5, units = "cm")
  
  # ----- extract peptides -----
  kk = which(res$estimate<(-fc.threshold) & res$padj<=pval.threshold)
  
  # ----- save plots -----
  if (length(kk) > 1) {
    
    pdf(paste0("results/",subunit,"/",comparison,"/DEanalysis_",comparison,".pdf"), height = 12, width = 6)
    par(mfrow = c(4,1))
    plotDensities(FC[kk,c(b1_idx,b2_idx,b5_idx,no_idx)],
                  group = c(rep("b1",length(b1_idx)),
                            rep("b2",length(b2_idx)),
                            rep("b5",length(b5_idx)),
                            rep("no",length(no_idx))),
                  col = c("black", "lightblue", "purple","gray"),
                  legend = F,
                  main = "LOG FOLD CHANGE - b5: purple, b1: black, b2: blue, no inhibitor: gray")
    boxplot(FC[kk,], main = "LOG FOLD CHANGE")
    
    plotDensities(log10(INT_4hrs[kk, 4+c(b1_idx,b2_idx,b5_idx,no_idx)]+1),  # !!!!!!
                  group = c(rep("b1",length(b1_idx)),
                            rep("b2",length(b2_idx)),
                            rep("b5",length(b5_idx)),
                            rep("no",length(no_idx))),
                  col = c("black", "lightblue", "purple","gray"),
                  legend = F,
                  main = "LOG10 INTENSITY @ 4HRS - b5: purple, b1: black, b2: blue, no inhibitor: gray")
    boxplot(log10(INT_4hrs[kk, intIDX]+1), main = "LOG10 INTENSITY @ 4HRS")
    
    dev.off()
  }
  
  
  # p-value diagnostics and volcanoplots
  pdf(paste0("results/",subunit,"/",comparison,"/DEanalysis_",comparison,"_stats.pdf"), height = 6, width = 10)
  
  print(Volcanoplots)
  print(plot(qobj))
  print(hist(qobj))
  print(gridExtra::grid.arrange(pval,padj, nrow=2))
  
  dev.off()
  
  if (length(kk) > 0) {
    ### OUTPUT ###
    write.csv(INT_4hrs[kk, ],
              paste0("results/",subunit,"/",comparison,"/DEanalysis_",comparison,"_intensities.csv"),
              row.names = F)
    write.csv(FC[kk,],
              paste0("results/",subunit,"/",comparison,"/DEanalysis_",comparison,"_fold-changes.csv"),
              row.names = T)
  }
  
  # save statistics
  write.csv(res,
            paste0("results/",subunit,"/",comparison,"/DEanalysis_",comparison,"_statistics.csv"),
            row.names = T)
  
}
################################################################################

# ----- APPLY -----


for (su in names(comp)) {
  
  cnt = comp[[su]]
  for (cn in 2:length(cnt)) {
    
    DEanalysis(subunit = su, comparison = names(cnt)[cn])
    
  }
  
}



