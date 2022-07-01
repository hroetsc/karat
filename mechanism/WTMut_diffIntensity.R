### karat projetc - PCPS mechanism ###
# description:  differential intensity analysis on WT/Mut pairs using limma
# input:        matrix with log-fold changes of peptides identified in WT/Mut,
#               table with WT/Mut pairs
# output:       up/downregulated peptides
# author:       HPR, modified from YH

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

suppressWarnings(dir.create("results/WTMut/DEanalysis/"))


### HYPERPARAMETERS ###
# Differential expression parameters
pval.threshold = 0.1
fc.threshold = 1.5


### INPUT ###
Y = read.csv("results/WTMut/WT-Mut-pairs+intensities.csv", stringsAsFactors = F)

### MAIN PART ###

DEanalysis = function(origSubs) {
  
  load(paste0("results/WTMut/fold-changes/",origSubs,".RData"))
  
  # ----- get comparisons -----
  cnames = str_split_fixed(str_split_fixed(colnames(FC), "-", Inf)[,1], "_", Inf)[,2]
  comp = combn(unique(cnames), 2)
  
  for (j in 1:ncol(comp)) {
    
    gr1 = which(cnames == comp[1,j])
    gr2 = which(cnames == comp[2,j])
    
    gr1Name = comp[1,j]
    gr2Name = comp[2,j]
    
    # ----- create contrast matrix -----
    
    L = FC[, c(gr1,gr2)] %>% as.data.frame()
    names(L) = c(rep(gr1Name, length(gr1)), rep(gr2Name, length(gr2)))
    sample_info = data.frame(condition = names(L))
    design = model.matrix(~condition, sample_info)
    
    # ----- fit linear model -----
    data_fit = lmFit(object = L,
                     design = design,
                     method = "ls")
    
    eb = eBayes(data_fit)
    topTable(eb, adjust.method = "BH")
    res <- tidy(eb)
    
    ## check statistics
    qobj <- try(qvalue(p = res$p.value))
    # error if dataset is too small to estimate pi0 from the data
    if ("try-error" %in% class(qobj)) {
      print("not enough data to estimate pi0 reliably")
      qobj <- qvalue(p = res$p.value, pi0 = 1)
    }
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
    
    # Export DE stats
    DE_stats <- res %>%
      dplyr::select(c("threshold","term")) %>%
      table() %>%
      as.data.frame.matrix() %>%
      t()
    DE_stats <- cbind(Contrast = rownames(DE_stats),
                      DE_stats) %>%
      as.data.frame()
    print(DE_stats)
    
    # ----- plot and export results -----
    
    Volcanoplots <- res %>%
      ggplot(aes(y=-log10(padj), x=estimate)) + 
      geom_point(stat="identity", position="identity", alpha=0.5, size=1) + 
      theme_bw() + 
      theme(text=element_text(family="sans", face="plain", color="#000000", size=12, hjust=0.5, vjust=0.5)) + 
      xlab("log2 fold change") + 
      ylab("-log10 p.adj") +
      ggtitle(paste0(origSubs,": ", gr1Name, " vs. ", gr2Name),
              sub = paste0("differential peptide intensities - down: ",DE_stats$down,", up: ", DE_stats$up)) +
      geom_hline(aes(yintercept = -log10(pval.threshold)), colour = "black") +
      geom_vline(xintercept = c(-fc.threshold,fc.threshold), colour = "black", linetype = "dashed") +
      theme_classic()
    
    
    # export significant results
    kk = which(abs(res$estimate)>fc.threshold & res$padj<=pval.threshold)
    signif = res[kk,] %>%
      rename(pepPos = gene)
    
    ZZ = left_join(signif,
                   Y %>% select(pepPos, identical, pepSeq, spliceType, positions,protein_name)) %>%
      mutate(origSubs = origSubs,
             group1 = gr1Name,
             group2 = gr2Name) %>%
      unique()
    write.csv(ZZ, file = paste0("results/WTMut/DEanalysis/SIGNIF_",origSubs,"_",gr1Name,"-vs-",gr2Name,".csv"), row.names = F)
    
    # export plots
    pdf(paste0("results/WTMut/DEanalysis/",origSubs,"_",gr1Name,"-vs-",gr2Name,"_stats.pdf"), height = 6, width = 10)
    
    print(Volcanoplots)
    print(plot(qobj))
    print(hist(qobj))
    print(gridExtra::grid.arrange(pval,padj, nrow=2))
    
    dev.off()
    
    
  }
  
}

DEanalysis("MPL")
DEanalysis("JAK2")
DEanalysis("IDH1")
DEanalysis("BRAF")
DEanalysis("NRAS")
DEanalysis("KRAS")


