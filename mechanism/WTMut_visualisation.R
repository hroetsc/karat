### karat projetc - PCPS mechanism ###
# description:  fetch differentially produced peptides and link to mutation position
# input:        differentially produced WT/Mut pairs, table with WT/Mut pairs
# output:       mutation effects at different positions
# author:       HPR

library(dplyr)
library(stringr)
library(ggplot2)
library(ggseqlogo)
library(RColorBrewer)
library(lattice)
source("src/_extract-aa.R")
source("src/invitroSPI_utils.R")

# rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
rf <- colorRampPalette(rev(brewer.pal(11,'RdBu')))
rcol <- rf(100)

AAchar_here = c("P","G","C","M","A","V","I","L","F","Y","W","H","R","K","D","E","N","Q","S","T","X")
AAchar_here_sorted = sort(AAchar_here)

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

AAcolours = c(rep("deeppink",3),  # special case
              rep("orange", 5),  # hydrophobic
              rep("mediumpurple", 3),  # aromatic
              rep("seagreen3",1), rep("seagreen", 1), rep("seagreen2",1),  # basic
              rep("firebrick",2),  # acidic
              rep("dodgerblue",4))  # nucleophilic/polar
names(AAcolours) = AAnames


# create colour scheme
cs = make_col_scheme(chars = c("P","G","C",# special case
                               "M", "A","V","I","L", # hydrophobic
                               "F","Y","W",  # aromatic
                               "H","R","K",  # basic
                               "D","E",  # acidic
                               "N","Q","S","T"),  # nucleophilic/polar
                     cols = c(rep("deeppink",3),
                              rep("orange", 5),
                              rep("mediumpurple", 3),
                              rep("seagreen3",3),
                              rep("firebrick",2),
                              rep("dodgerblue",4)),
                     groups = c(rep("special",3),
                                rep("hydrophobic", 5),
                                rep("aromatic", 3),
                                rep("basic",3),
                                rep("acidic",2),
                                rep("polar",4)))

# ----- load all pairs -----

pairTBL = lapply(fs, read.csv, stringsAsFactors = F) %>%
  plyr::ldply()

origSubs = Y$origSubs %>% unique()
allMaster = list()
counter = 1

pdf("results/WTMut/PAIRS/foldChanges.pdf", height = 6, width = 15)
for (o in origSubs) {
  
  signifPos = pairTBL$pepPos[pairTBL$origSubs == o]
  
  # load fold changes
  load(paste0("results/WTMut/fold-changes/",o,".RData"))
  L = FC[rownames(FC) %in% signifPos, ] %>% as.data.frame()

  Lg = L %>% 
    tibble::rownames_to_column("pepPos")
  if (nrow(Lg) == 0) {
    next
  }
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
    mutate(n = length(unique(pepPos))) %>%
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
      ylab("log10 fold change") +
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
      ylab("log10 fold change") +
      ggtitle("PCP", subtitle = paste(unique(cntMaster$protein_name), collapse = ", "))
  }
    
  if (is.ggplot(psp) & is.ggplot(pcp)) {
    gridExtra::grid.arrange(psp,pcp, nrow = 2) %>% print()
  } else if (!is.ggplot(psp) & is.ggplot(pcp)) {
    print(pcp)
  } else {
    print(psp)
  }
  
  allMaster[[counter]] = cntMaster
  counter = counter +1
}
dev.off()


# ----- plot per amino acid / position -----

allMaster = plyr::ldply(allMaster)
# reorder amino acids
allMaster$aa = factor(allMaster$aa, levels = AAchar_here)


# plot per position
allP = list()
counter = 1
pos = c("P8","P7","P6","P5","P4", "P3", "P2", "P1", "P-8", "P-7", "P-6", "P-5", "P-4", "P-3","P-2", "P-1")
for (p in 1:length(pos)) {
  cntR = allMaster[allMaster$identical == as.character(pos[p]), ]
  
  if (nrow(cntR) > 0) {
    cntP = ggplot(cntR, aes(x = aa, y = value, fill = productType)) +
      geom_boxplot(outlier.shape = NA, coef = 0, position = position_dodge(0.5), aes(alpha = .8)) +
      scale_fill_manual(values = c(plottingCols[["PCP"]], plottingCols[["PSP"]])) +
      scale_x_discrete(limits = AAchar_here) +
      xlab("amino acid") +
      ylab("log10 fold change") +
      ggtitle(pos[p])
    
    allP[[counter]] = cntP
    counter = counter + 1
  }
}

ggsave(filename = "results/WTMut/PAIRS/fC_perPos.pdf", 
       plot = gridExtra::marrangeGrob(allP, nrow=8, ncol=2, byrow = F), 
       width = 15, height = 27, dpi = "retina")

# plot per amino acid
allPaa = list()
counter = 1
allMaster$identical = factor(allMaster$identical, levels = SRnames)
for (a in 1:length(AAchar_here)) {
  
  cntR = allMaster[allMaster$aa == as.character(AAchar_here[a]), ]
  
  if (nrow(cntR) > 0) {
    cntP = ggplot(cntR, aes(x = identical, y = value, fill = productType)) +
      geom_boxplot(outlier.shape = NA, coef = 0, position = position_dodge(0.5), aes(alpha = .8)) +
      scale_fill_manual(values = c(plottingCols[["PCP"]], plottingCols[["PSP"]])) +
      scale_x_discrete(limits = c("P8","P7","P6","P5","P4", "P3", "P2", "P1",
                                  "P-1","P-2","P-3","P-4","P-5","P-6","P-7","P-8")) +
      xlab("position") +
      ylab("log10 fold change") +
      ggtitle(AAchar_here[a])
    
    allPaa[[counter]] = cntP
    counter = counter + 1
  }
  
}

ggsave(filename = "results/WTMut/PAIRS/fC_perAA.pdf",
       plot = gridExtra::marrangeGrob(allPaa, nrow=3, ncol=4, byrow = T), 
       width = 28, height = 10, dpi = "retina")


# ----- get logoplots -----

getLogos = function(prodType, interesting_residues = SRnames) {
  
  # get frequency matrix
  freqMat = allMaster %>%
    filter(productType == prodType) %>%
    group_by(identical, aa) %>%
    summarise(med = median(value))
  
  c = log(median(exp(freqMat$med)))
  freqMat$med = freqMat$med - c
  
  DFt = freqMat %>%
    tidyr::spread(identical,med) %>%
    tibble::column_to_rownames("aa")
  
  
  # mising positions / aa
  lack = interesting_residues[!interesting_residues %in% colnames(DFt)]
  if (length(lack) > 0) {
    DFt = cbind(DFt, matrix(NA, nrow(DFt), length(lack)))
    colnames(DFt)[(ncol(DFt)-length(lack)+1):ncol(DFt)] = lack
  }
  
  DFt = as.matrix(DFt[,interesting_residues])
  
  # logoplots
  
  l1 = ggseqlogo(DFt, method = "custom", seq_type = "aa", font = "helvetica_bold",col_scheme = cs) +
    scale_x_discrete(limits = colnames(DFt)) +
    geom_vline(xintercept = c(8.5,16.5,24.5), lty = "dashed", col = "black") +
    xlab("position") +
    ylab("log10 fold change") +
    ggtitle(prodType)
  
  l2 = ggseqlogo(DFt-min(DFt,na.rm = T), method = "custom", seq_type = "aa", font = "helvetica_bold",col_scheme = cs) +
    scale_x_discrete(limits = colnames(DFt)) +
    geom_vline(xintercept = c(8.5,16.5,24.5), lty = "dashed", col = "black") +
    xlab("position") +
    ylab("entropy [bits]") +
    ggtitle(prodType)
  
  
  return(list(logo1 = l1, logo2 = l2))
}

pcpLogo = getLogos("PCP", interesting_residues = names(PCPpos))
pspLogo = getLogos("PSP", interesting_residues = SRnames)

ggsave(filename = "results/WTMut/PAIRS/logoPlots_1.png",
       plot = gridExtra::grid.arrange(pcpLogo$logo1, pspLogo$logo1),
       height = 10, width = 15, dpi = "retina")
ggsave(filename = "results/WTMut/PAIRS/logoPlots_2.png",
       plot = gridExtra::grid.arrange(pcpLogo$logo2, pspLogo$logo2),
       height = 10, width = 15, dpi = "retina")


# ----- fold change upon aa swap -----
interesting_residues = SRnames

protNames = unique(allMaster$protein_name)
protNames = sapply(protNames, function(p){
  substr(Y$substrateSeq[Y$protein_name == p][1],
         Y$residue[Y$protein_name == p][1],
         Y$residue[Y$protein_name == p][1])
}) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("protein_name")
names(protNames)[2] = "aa"
protNames$origSubs = str_split_fixed(protNames$protein_name,"_",Inf)[,1]
protNames$group = str_split_fixed(protNames$protein_name,"_",Inf)[,2]

SIGNIF = lapply(fs, function(x){
  print(basename(x))
  S = read.csv(x, stringsAsFactors = F) %>%
    select(origSubs, pepPos, spliceType, group1, group2,estimate, identical) %>%
    tidyr::separate_rows(identical, sep = ";") %>%
    unique()  # %>%
    # group_by(origSubs, group1, group2, identical) %>%
    # summarise(mean_fc = mean(estimate))
  
  if (nrow(S) > 0) {
    S$aa1 = protNames$aa[protNames$origSubs == S$origSubs[1] & protNames$group == S$group1[1]]
    S$aa2 = protNames$aa[protNames$origSubs == S$origSubs[1] & protNames$group == S$group2[1]]
    S$swap = paste0(S$aa2, "->", S$aa1)
    
    return(S)
  }
  
})

SIGNIF = plyr::ldply(SIGNIF)
SIGNIF$swap %>% unique() %>% sort()

# FC(B->A) = - FC(A->B)/(FC(A->B)+1)
# SIGNIF$estimate[SIGNIF$swap == "R->G"] = -1*(SIGNIF$estimate[SIGNIF$swap == "R->G"]/(SIGNIF$estimate[SIGNIF$swap == "R->G"]+1))
# SIGNIF$swap[SIGNIF$swap == "R->G"] = "G->R"

SIGNIF = SIGNIF[-which(SIGNIF$swap %in% c("D->D", "D->R")), ]
lim = c(min(SIGNIF$estimate), max(SIGNIF$estimate))

# swaps = c("L->W","F->V","E->V","D->R","D->G","R->G","R->Q","K->Q","K->R","H->R")
swaps = c("W->L","V->F","V->E","R->D","G->D","G->R","Q->R","Q->K","R->K","R->H")

getHeatmap = function(prodType, interesting_residues) {
  
  DFt = SIGNIF %>%
    select(-pepPos) %>%
    mutate(productType = ifelse(spliceType == "PCP", "PCP", "PSP")) %>%
    filter(productType == prodType) %>%
    group_by(swap, identical) %>%
    summarise(mean_fc = median(estimate),
              numPairs = n())
  N = DFt
  DFt = DFt %>%
    select(-numPairs) %>%
    tidyr::spread(identical,mean_fc) %>%
    tibble::column_to_rownames("swap")
  
  N = N %>%
    select(-mean_fc) %>%
    tidyr::spread(identical,numPairs) %>%
    tibble::column_to_rownames("swap")
  
  lack = interesting_residues[!interesting_residues %in% colnames(DFt)]
  if (length(lack) > 0) {
    DFt = cbind(DFt, matrix(NA, nrow(DFt), length(lack)))
    N = cbind(N, matrix(NA, nrow(N), length(lack)))
    colnames(DFt)[(ncol(DFt)-length(lack)+1):ncol(DFt)] = lack
    colnames(N)[(ncol(N)-length(lack)+1):ncol(N)] = lack
  }
  lack = swaps[!swaps %in% rownames(DFt)]
  if (length(lack) > 0) {
    DFt = rbind(as.matrix(DFt), matrix(NA, length(lack), ncol(DFt)))
    N = rbind(as.matrix(N), matrix(NA, length(lack), ncol(N)))
    rownames(DFt)[(nrow(DFt)-length(lack)+1):nrow(DFt)] = lack
    rownames(N)[(nrow(N)-length(lack)+1):nrow(N)] = lack
  }
  
  DFt = as.matrix(DFt[rev(swaps),interesting_residues])
  N = as.matrix(N[rev(swaps),interesting_residues])
  
  colnames(DFt) = str_replace_all(colnames(DFt), coll("_"), coll("'"))
  
  
  levelplot(DFt %>% t(), axes = F, col.regions = rcol,
            at = seq(lim[1], lim[2], length.out = 100),
            pretty = T, labels = T,
            main = "mean log fold change upon mutation",
            xlab = "position", ylab = "mutation",
            sub = prodType) %>% print()
  
  # levelplot(N %>% t(), axes = F, col.regions = colorRampPalette(rev(brewer.pal(11,'PiYG')))(100),
  #           pretty = T, labels = T,
  #           main = "number of peptide pairs",
  #           xlab = "position", ylab = "mutation",
  #           sub = prodType) %>% print()
  
  print(sum(N, na.rm = T))
  N[is.na(N)] = ""
  
  grid::grid.newpage()
  N[swaps, ] %>% gridExtra::grid.table()
  
  return(DFt)
}

pdf("results/WTMut/PAIRS/HEATMAP2.pdf", height = 6, width = 15)
DFtpcp = getHeatmap("PCP", interesting_residues = names(PCPpos))
DFtpsp = getHeatmap("PSP", interesting_residues = SRnames)
dev.off()

# TO-DO: nice plots for thesis
png("results/WTMut/PAIRS/_forThesis_heatmapPSP.png", height = 5, width = 15, units = "in", res = 600)
levelplot(DFtpsp %>% t(), axes = F, col.regions = rcol,
          at = seq(lim[1], lim[2], length.out = 100),
          panel = function(...) {
            panel.levelplot(...)
            panel.abline(v=8.5, lty = "dashed")
            panel.abline(v=16.5, lty = "dashed")
            panel.abline(v=24.5, lty = "dashed")
          },
          pretty = T, labels = T) %>% print()
dev.off()

png("results/WTMut/PAIRS/_forThesis_heatmapPCP.png", height = 5, width = 7.5, units = "in", res = 600)
levelplot(DFtpcp %>% t(), axes = F, col.regions = rcol,
          at = seq(lim[1], lim[2], length.out = 100),
          panel = function(...) {
            panel.levelplot(...)
            panel.abline(v=8.5, lty = "dashed")
          },
          pretty = T, labels = T) %>% print()
dev.off()


# ----- plot selected kinetics -----
# pick an example for thesis
krasPSP = SIGNIF[SIGNIF$swap %in% c("H->R","K->R") & SIGNIF$identical == "P1", ]
krasPSP = Y[Y$pepPos %in% c("IDH1-13_14_3_8","IDH1-13_14_4_8"), ]

krasPCP = SIGNIF[SIGNIF$identical == "P1" & SIGNIF$spliceType == "PCP", ]
krasPCP = Y[Y$pepPos %in% krasPCP$pepPos, ]
krasPCP = krasPCP[krasPCP$protein_name %in% c("KRAS_WT","KRAS_G12R","KRAS_G12D") & krasPCP$positions == "1_11", ]

# screening
pdf("results/WTMut/PAIRS/example_plots.pdf", height = 5, width = 5)
pos = krasPCP$positions %>% unique()
for (i in 1:length(pos)) {
  cnt = krasPCP[krasPCP$positions == pos[i], ]
  tmp = cnt %>%
    group_by(pepSeq,spliceType,positions,biological_replicate,digestTime, protein_name) %>%
    summarise(mean_int = if (all(intensity == 0) & all(digestTime != 0)) 0 else mean(intensity[intensity!=0 | digestTime == 0], na.rm=T),
              sd_int = if (all(intensity == 0) & all(digestTime != 0)) 0 else sd(intensity[intensity!=0 | digestTime == 0], na.rm=T)) %>%
    arrange(digestTime) %>%
    suppressMessages ()
  
  # Cols = c("darkorange","firebrick")
  # Cols = c("royalblue","lightblue")
  nCol = length(unique(tmp$protein_name))
  Cols = colorRampPalette(rev(brewer.pal(9,'Blues')))(nCol)
  # Cols = Cols[1:nCol]
  
  tmp = data.frame(protein_name = unique(tmp$protein_name), col = as.character(Cols)) %>%
    right_join(tmp) %>%
    suppressMessages()
  
  plot(x = tmp$digestTime, y = tmp$mean_int,
       pch = 16, cex = 1.5,
       xlab = "time [hrs]", ylab = "MS1 intensity",
       col = tmp$col,
       main = paste(unique(tmp$pepSeq), collapse = ", "), 
       sub = paste(unique(tmp$spliceType),cnt$identical[1], cnt$positions[1], sep = ", "))
  
  # ----- add lines (bioReps) -----
  tmp = tmp %>%
    mutate(combo = paste(biological_replicate, protein_name))
  for (gr in unique(tmp$combo)) {
    lines(x=tmp$digestTime[tmp$combo == gr],
          y=tmp$mean_int[tmp$combo == gr],
          col = tmp$col[tmp$combo == gr][1], lwd=1.5)
    
  }
  
  # ----- add legend -----
  legend(x = "topleft", legend = c(unique(tmp$protein_name)), 
         col = Cols,
         pch = 19, cex = 1.2, bty = "n", lty = "solid")
}
dev.off()


# plot the example kinetics
pos = krasPCP$positions %>% unique()
for (i in 1:length(pos)) {
  cnt = krasPCP[krasPCP$positions == pos[i], ]
  
  png(paste0("results/WTMut/PAIRS/_forThesis_",cnt$pepPos[1],".png"), height = 5, width = 5, units = "in", res = 600)
  
  tmp = cnt %>%
    group_by(pepSeq,spliceType,positions,biological_replicate,digestTime, protein_name) %>%
    summarise(mean_int = if (all(intensity == 0) & all(digestTime != 0)) 0 else mean(intensity[intensity!=0 | digestTime == 0], na.rm=T),
              sd_int = if (all(intensity == 0) & all(digestTime != 0)) 0 else sd(intensity[intensity!=0 | digestTime == 0], na.rm=T)) %>%
    arrange(digestTime) %>%
    suppressMessages ()
  
  tmp$mean_int[tmp$digestTime == 0] = 0
  
  
  Cols = c("darkorange","firebrick","goldenrod2")
  # Cols = c("royalblue","skyblue")
  nCol = length(unique(tmp$protein_name))
  # Cols = colorRampPalette(rev(brewer.pal(9,'Blues')))(nCol)
  Cols = Cols[1:nCol]
  
  tmp = data.frame(protein_name = unique(tmp$protein_name), col = as.character(Cols)) %>%
    right_join(tmp) %>%
    suppressMessages()
  
  plot(x = tmp$digestTime, y = tmp$mean_int,
       pch = 16, cex = 1.5,
       xlab = "time [hrs]", ylab = "MS1 intensity",
       col = tmp$col,
       main = paste(unique(tmp$pepSeq), collapse = ", "), 
       sub = paste(unique(tmp$spliceType),cnt$identical[1], cnt$positions[1], sep = ", "))
  
  # ----- add lines (bioReps) -----
  tmp = tmp %>%
    mutate(combo = paste(biological_replicate, protein_name))
  for (gr in unique(tmp$combo)) {
    lines(x=tmp$digestTime[tmp$combo == gr],
          y=tmp$mean_int[tmp$combo == gr],
          col = tmp$col[tmp$combo == gr][1], lwd=1.5)
    
  }
  
  # ----- add legend -----
  legend(x = "topleft", legend = c(unique(tmp$protein_name)), 
         col = Cols,
         pch = 19, cex = .9, bty = "n", lty = "solid")
  
  dev.off()
}

### OUTPUT ####

