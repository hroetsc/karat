### karat projetc - PCPS mechanism ###
# description:  fetch differentially produced peptides and link to mutation position
# input:        differentially produced WT/Mut pairs, table with WT/Mut pairs
# output:       mutation effects at different positions
# author:       HPR

library(dplyr)
library(stringr)
library(ggplot2)
source("src/_extract-aa.R")

theme_set(theme_classic())

### INPUT ###
Y = read.csv("results/WTMut/WT-Mut-pairs+intensities.csv", stringsAsFactors = F)
fs = list.files("results/WTMut/DEanalysis/", pattern = "SIGNIF", recursive = T, full.names = T)


### MAIN PART ###
suppressWarnings(dir.create("results/WTMut/PAIRS/"))

# get product types back
Y = Y %>%
  mutate(productType = ifelse(spliceType == "PCP", "PCP", "PSP"))

# ----- define amino acid group colours -----
AAnames = c("P","G","C",# special case
              "M", "A","V","I","L", # hydrophobic
              "F","Y","W",  # aromatic
              "H","R","K",  # basic
              "D","E",  # acidic
              "N","Q","S","T")  # nucleophilic

AAcolours = c(rep("lawngreen",3),  # special case
              rep("orange", 5),  # hydrophobic
              rep("mediumpurple", 4),  # aromatic
              rep("seagreen", 1), rep("seagreen2",1),  # basic
              rep("firebrick",2),  # acidic
              rep("dodgerblue",4))  # nucleophilic/polar
names(AAcolours) = AAnames

# ----- load all pairs -----

pairTBL = lapply(fs, read.csv, stringsAsFactors = F) %>%
  plyr::ldply()

origSubs = Y$origSubs %>% unique()

pdf("results/WTMut/PAIRS/foldChanges.pdf", height = 6, width = 15)
for (o in origSubs) {
  
  signifPos = pairTBL$pepPos[pairTBL$origSubs == o]
  
  # load fold changes
  load(paste0("results/WTMut/fold-changes/",o,".RData"))
  L = FC[rownames(FC) %in% signifPos, ] %>% as.data.frame()

  Lg = L %>% 
    tibble::rownames_to_column("pepPos")
  Lg = Lg %>% tidyr::gather(key, value, -pepPos)
  Lg$protein_name = str_split_fixed(Lg$key, "-", Inf)[,1]
  Lg$key = str_split_fixed(str_split_fixed(Lg$key, "-", Inf)[,1], "_", Inf)[,2]
  
  # get the pairs
  # info about position
  cntPairs = Y[Y$pepPos %in% signifPos, ] %>%
    select(productType,pepSeq, origSubs, protein_name, pepPos, identical) %>%
    tidyr::separate_rows(identical,sep = ";")
  
  # join both
  cntMaster = inner_join(cntPairs, Lg)
  cntMaster$identical = factor(cntMaster$identical, levels = SRnames)
  
  # get the mutated amino acids
  protNames = unique(cntMaster$protein_name)
  cntMaster = sapply(protNames, function(p){
    substr(Y$substrateSeq[Y$protein_name == p][1],
           Y$residue[Y$protein_name == p][1],
           Y$residue[Y$protein_name == p][1])
  }) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("protein_name") %>%
    right_join(cntMaster)
  names(cntMaster)[2] = "aa"
  cntMaster$aa = factor(cntMaster$aa, levels = AAnames)
  
  # add n-numbers
  cntMaster = cntMaster %>%
    group_by(productType, identical) %>%
    mutate(n = length(unique(pepSeq))) %>%
    ungroup()
  
  # plot
  psp = cntMaster %>%
    filter(productType == "PSP")
  if(nrow(psp) > 0) {
    psp = ggplot(psp, aes(x = identical, y = value, fill = aa)) +
      geom_boxplot(aes(alpha = .9), position = position_dodge(0.8)) +
      geom_vline(xintercept = c(8.5,16.5,24.5), lty = "dashed", col = "gray") +
      scale_x_discrete(limits = SRnames) +
      scale_fill_manual(values = AAcolours[unique(cntMaster$aa)][1:length(unique(cntMaster$aa[cntMaster$productType == "PSP"]))]) +
      ylim(c(-1.5,max(psp$value))) +
      geom_text(aes(label = factor(n)), y=-1, size = 3, stat = "count") +
      xlab("mutation position") +
      ylab("log fold change") +
      ggtitle("PSP", subtitle = paste(unique(cntMaster$protein_name), collapse = ", "))
  }
    
  pcp = cntMaster %>%
    filter(productType == "PCP")
  if (nrow(pcp) > 0) {
  pcp = ggplot(pcp, aes(x = identical, y = value, fill = aa)) +
      geom_boxplot(aes(alpha = .9), position = position_dodge(0.8)) +
      geom_vline(xintercept = c(8.5), lty = "dashed", col = "gray") +
      scale_x_discrete(limits = names(PCPpos)) +
      scale_fill_manual(values = AAcolours[unique(cntMaster$aa)][1:length(unique(cntMaster$aa[cntMaster$productType == "PCP"]))]) +
      ylim(c(-1.5,max(pcp$value))) +
      geom_text(aes(label = factor(n)), y=-1, size = 3, stat = "count") +
      xlab("mutation position") +
      ylab("log fold change") +
      ggtitle("PCP", subtitle = paste(unique(cntMaster$protein_name), collapse = ", "))
  }
    
  if (is.ggplot(psp) & is.ggplot(pcp)) {
    gridExtra::grid.arrange(psp,pcp, nrow = 2) %>% print()
  } else if (!is.ggplot(psp) & is.ggplot(pcp)) {
    print(pcp)
  } else {
    print(psp)
  }
  
}
dev.off()


### OUTPUT ####
