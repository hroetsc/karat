### karat projetc - PCPS mechanism ###
# description:  compare precursor lengths between PCPs and PSPs
# input:        invitroSPI: Roetschke SciData, EGFR ery, WT substrates
#               number of peptides @ particular SR lengths
# output:       ?
# author:       HPR

library(dplyr)
library(stringr)
library(ggplot2)
library(twosamples)
library(dgof)
source("src/invitroSPI_utils.R")
source("src/plotting_utils.R")
source("../brainstorming/src/number-of-products.R")

theme_set(theme_classic())
Lext = 1
pepLen_prots = seq(7,30)
pepLen_polypeps = seq(5,40)


### INPUT ###
load("data/invitroSPI.RData")
load("../../proteinsPCPS/new/data/aSPIre.RData")

### MAIN PART ###
suppressWarnings(dir.create("results/precursor_lengths/"))

# ------ preprocessing -----
polypeps = ProteasomeDB %>%
  ILredundancy() %>%
  disentangleMultimappers.Type() %>%
  disentangleMultimappers.SRlen() %>%
  removeMultimappers.Type() %>%
  removeMultimappers.SRlen() %>%
  uniquePeptides()

proteins = Kinetics %>%
  ILredundancy() %>%
  disentangleMultimappers.Type() %>%
  disentangleMultimappers.SRlen() %>%
  removeMultimappers.Type() %>%
  removeMultimappers.SRlen() %>%
  uniquePeptides()

# ----- observed precursor length -----
PrecursorLength = function(trueDB) {
  
  pos = str_split_fixed(trueDB$positions, "_", Inf)[,c(1:4)]
  pos = apply(pos,2,as.numeric)
  
  trueDB$prec = trueDB$pepSeq
  trueDB$precLen = nchar(trueDB$pepSeq)
  trueDB$nucLen = nchar(trueDB$pepSeq)
  
  k = which(trueDB$productType == "PSP")
  trueDB$prec[k] = substr(trueDB$substrateSeq[k], pos[k,1], pos[k,2])
  trueDB$precLen[k] = pos[k,2]-pos[k,1]+1
  trueDB$nucLen[k] = pos[k,4]-pos[k,3]+1
  
  return(trueDB)
}

polypeps = PrecursorLength(polypeps)
proteins = PrecursorLength(proteins)

# ----- theoretical precursor length -----

getPrecursorComparison = function(DB, pepLen, dataset, reactant = "SR1") {
  
  OBSERVED = DB %>%
    rowwise() %>%
    mutate(L = nchar(substrateSeq),
           N = nchar(pepSeq),
           sr1 = ifelse(reactant == "SR1", precLen, nucLen)) %>%
    filter(sr1 %in% pepLen) %>%  # make sure to compare SR1s only to PCPs that would be detectable in MS
    group_by(substrateID,spliceType,L,N,sr1) %>%
    summarise(numPrec = dplyr::n())
  
  THEORETICAL = DB %>%
    mutate(L = nchar(substrateSeq)) %>%
    group_by(spliceType,substrateID,L) %>%
    summarise(N = paste(pepLen, collapse = ";")) %>%
    tidyr::separate_rows(N, sep = ";") %>%
    mutate(L = as.numeric(L),
           N = as.numeric(N)) %>%
    ungroup() %>% rowwise() %>%
    filter(N >= pepLen[1]+Lext) %>%
    mutate(sr1 = paste(seq(pepLen[1], N-Lext), collapse = ";")) %>%
    tidyr::separate_rows(sr1, sep = ";") %>%
    mutate(sr1 = as.numeric(sr1)) %>%
    filter(sr1 %in% pepLen) %>% 
    ungroup() %>% rowwise() %>%
    mutate(numPrec = ifelse(spliceType == "cis", numCis_fixedSR(L,N), numRevCis_fixedSR(L,N)),
           numPrec = ifelse(spliceType == "trans", numTrans_fixedSR(L,N,sr1), numPrec),
           numPrec = ifelse(spliceType == "PCP", numPCP(L,N), numPrec))
  
  # any(ALL$theonumPrec < 0)
  ALL = rbind(OBSERVED %>% mutate(df = "obsnumPrec"),
              THEORETICAL %>% mutate(df = "theonumPrec"))
  
  # scale between 1 and 100
  ALLdf = ALL %>%
    group_by(df) %>%
    mutate(normNumPrec = round((100-1)*((numPrec - min(numPrec))/(max(numPrec) - min(numPrec))) + 1))
  
  ALLdf_rep = data.frame(lapply(ALLdf, rep, ALLdf$normNumPrec))
  
  # ----- plot distributions -----
  # violin plot
  ALLdf_rep$spliceType = factor(ALLdf_rep$spliceType, levels = c("PCP","cis","revCis","trans"))
  # if (!reactant == "SR1") {
  #   ALLdf_rep = ALLdf_rep[-which(ALLdf_rep$spliceType == "PCP"), ]
  # }
  prl = ggplot(ALLdf_rep, aes(x = spliceType, y = sr1, fill = df)) +
    geom_split_violin(size = .25) +
    stat_summary(fun = median, fun.min = median, fun.max = median,
                 geom = "crossbar", 
                 width = .2,
                 position = position_dodge(width = .2)) +
    xlab("product type") +
    ylab("precursor length") +
    ggtitle("observed vs. calculated precursor length distributions", subtitle = dataset) +
    scale_fill_manual("database",
                      values = c(plottingCols[["ProteasomeDB"]], plottingCols[["randomDB"]]),
                      drop = F) +
    theme(axis.text.x = element_text(angle = 90))
  
  # all splice types
  precComp = ggplot(ALLdf_rep[ALLdf_rep$df == "obsnumPrec", ], aes(x = spliceType, y = sr1, fill = spliceType)) +
    geom_violin(draw_quantiles = 0.5) +
    scale_fill_manual(values = c(plottingCols[["PCP"]], plottingCols[["cis"]], plottingCols[["revCis"]], plottingCols[["trans"]])) +
    xlab("product type") +
    ylab("precursor length") +
    ggtitle("precursor length distributions", subtitle = dataset)
  
  sapply(c("cis","revCis","trans"), function(t){
    p = try(cvm.test(x = ALLdf_rep$sr1[ALLdf_rep$df == "obsnumPrec" & ALLdf_rep$spliceType == t],
                     y = ecdf(ALLdf_rep$sr1[ALLdf_rep$df == "obsnumPrec" & ALLdf_rep$spliceType == "PCP"]), type = "A2"))
    if ("try-error" %in% class(p)) {
      p = ad_test(a = ALLdf_rep$sr1[ALLdf_rep$df == "obsnumPrec" & ALLdf_rep$spliceType == t],
                  b = ALLdf_rep$sr1[ALLdf_rep$df == "obsnumPrec" & ALLdf_rep$spliceType == "PCP"])
    }
    return(p)
  }) %>%
    print()
  
  # distributions
  ALLsummarised = ALL %>%
    group_by(spliceType,df, sr1) %>%
    summarise(sumNum = sum(numPrec))
  
  if (reactant == "SR1") {
    types = c("PCP", "cis", "revCis", "trans")
  } else {
    types = c("cis", "revCis", "trans")
  }
  
  
  png(paste0("results/precursor_lengths/precDistr_",dataset,"_", reactant,".png"), height = 16, width = 8, units = "in", res = 300)
  par(mfrow = c(4,1))
  
  sapply(types, function(t){
    x1 = ALLsummarised$sr1[ALLsummarised$spliceType == t & ALLsummarised$df == "obsnumPrec"]
    y1 = ALLsummarised$sumNum[ALLsummarised$spliceType == t & ALLsummarised$df == "obsnumPrec"]
    plot(x = x1, y = y1,
         type = "l", col = plottingCols["ProteasomeDB"],
         xlab = "precursor length", ylab = "number of theoretical/detected sequences",
         sub = dataset, main = t,
         axes = F)
    polygon(c(x1, rev(x1)), c(y1, rep(0,length(y1))),
            col = add.alpha(plottingCols["ProteasomeDB"],0.5), lty = 0)
    axis(2, at = pretty(range(y1)))
    
    par(new = T)
    x2 = ALLsummarised$sr1[ALLsummarised$spliceType == t & ALLsummarised$df == "theonumPrec"]
    y2 = ALLsummarised$sumNum[ALLsummarised$spliceType == t & ALLsummarised$df == "theonumPrec"]
    plot(x = x2, y = y2, type = "l", col = plottingCols["randomDB"],
         axes = F, xlab = "", ylab = "")
    polygon(c(x2, rev(x2)), c(y2, rep(0,length(y2))),
            col = add.alpha(plottingCols["randomDB"],0.5), lty = 0)
    axis(4, at = pretty(range(y2)))
    axis(1)
    
    legend("topright", legend = c("identified", "theoretically"),
           col = c(plottingCols["ProteasomeDB"], plottingCols["randomDB"]), lty = c(1,1),
           cex = 0.8, bty = "n")
  })
  
  dev.off()
  
  # ----- statistical test -----
  
  print(dataset)
  print(reactant)
  sapply(types, function(t){
    p = try(cvm.test(x = ALLdf_rep$sr1[ALLdf_rep$df == "obsnumPrec" & ALLdf_rep$spliceType == t],
             y = ecdf(ALLdf_rep$sr1[ALLdf_rep$df == "theonumPrec" & ALLdf_rep$spliceType == t]), type = "A2"))
    if ("try-error" %in% class(p)) {
      p = ad_test(a = ALLdf_rep$sr1[ALLdf_rep$df == "obsnumPrec" & ALLdf_rep$spliceType == t],
                  b = ALLdf_rep$sr1[ALLdf_rep$df == "theonumPrec" & ALLdf_rep$spliceType == t])
    }
    return(p)
  }) %>%
    print()
  
  ggsave(filename = paste0("results/precursor_lengths/prec_",dataset,"_", reactant,".png"),
         plot = prl, dpi = "retina", height = 6, width = 6)
  ggsave(filename = paste0("results/precursor_lengths/TYPESprec_",dataset,"_", reactant,".png"),
         plot = precComp, dpi = "retina", height = 6, width = 6)
  
}


getPrecursorComparison(DB = polypeps, pepLen = pepLen_polypeps, dataset = "polypeps", reactant = "SR1")
getPrecursorComparison(DB = proteins, pepLen = pepLen_prots, dataset = "prots", reactant = "SR1")

getPrecursorComparison(DB = polypeps, pepLen = pepLen_polypeps, dataset = "polypeps", reactant = "SR2")
getPrecursorComparison(DB = proteins, pepLen = pepLen_prots, dataset = "prots", reactant = "SR2")
