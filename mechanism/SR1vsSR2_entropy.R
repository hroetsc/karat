### karat projetc - PCPS mechanism ###
# description:  sequence entropy for SR1 vs. SR2  (derived from random databases)
#               qualitative data set: Roetschke et al. SciData, EGFR, WT sequences of WT/Mut
#               random database for qualitative data set
# output:       entropy of different positions
# author:       HPR

library(dplyr)
library(stringr)
source("src/invitroSPI_utils.R")
source("src/plotting_utils.R")
source("src/_extract-aa.R")

theme_set(theme_classic())

### INPUT ###
load("data/invitroSPI.RData")
load("data/randomDB_smart.RData")

### MAIN PART ###
suppressWarnings(dir.create("results/SR1vsSR2/entropy/"))

# ----- preprocessing -----
Qual = ProteasomeDB %>%
  ILredundancy() %>%
  disentangleMultimappers.Type() %>%
  disentangleMultimappers.SRlen() %>%
  disentangleMultimappers.IVSeqLen() %>%
  removeMultimappers.SRlen() %>%
  removeMultimappers.AA() %>%
  uniquePeptides()

rndQual = rndDB %>%
  ILredundancy() %>%
  disentangleMultimappers.Type() %>%
  disentangleMultimappers.SRlen() %>%
  disentangleMultimappers.IVSeqLen() %>%
  removeMultimappers.SRlen() %>%
  removeMultimappers.AA() %>%
  uniquePeptides()

extractSRs = function(DB) {
  
  pos = str_split_fixed(DB$positions, "_", Inf)[,c(1:4)]
  pos = apply(pos,2,as.numeric)
  
  DB$sr1 = substr(DB$substrateSeq, pos[,1], pos[,2])
  DB$sr2 = substr(DB$substrateSeq, pos[,3], pos[,4])
  
  return(DB)
}

Qual = extractSRs(Qual)
rndQual = extractSRs(rndQual)

# ------ get entropy -----
# amino acids at relevant positions
DBQual = left_join(Qual, extract_aminoacids(Qual))
rndDBQual = left_join(rndQual, extract_aminoacids(rndQual))

same_pos = function(tbl) {
  
  if (! all(AA %in% names(tbl))) {
    k = which(! AA %in% names(tbl))
    
    for (i in 1:length(k)) {
      tbl = append(tbl, 0, k[i]-1)
    }
    names(tbl) = AA
  }
  
  return(tbl)
}

AAchar_here = c("P","G","A","V","L","M","F","Y","W","H","R","K","D","E","N","Q","S","T","C","X")

getEntropy = function(DB, rndDB) {
  
  getFreq = function(df) {
    freq = apply(df[,SRnames],2,function(x){
      y = table(x)
      names(y)[names(y) == ""] = "X"
      z = same_pos(y)
      z = z[match(AAchar_here,names(z))]
      # print(names(z))
      return(z)
    }) %>%
      as.matrix()
    rownames(freq) = AAchar_here
    freq[is.na(freq)] = 0
    
    Fs = apply(freq,2,function(x){
      x/sum(x)
    })
    # Fs = freq
    return(Fs)
  }
  
  htrue = getFreq(DB) %>% t()
  hrnd = getFreq(rndDB) %>% t()
  
  # Hs = htrue[SRnames,AAchar]
  Hs = htrue[SRnames,AAchar]/hrnd[SRnames,AAchar]
  # Hs = apply(Hs, 2, function(x){
  #     x/sum(x)
  # })
  
  H = apply(Hs,1,function(pos){
    pos = (pos-min(pos, na.rm = T))/(max(pos, na.rm = T) - min(pos, na.rm = T))
    return(-1*sum(pos*log2(pos), na.rm = T))
  })
  
  
  return(H)
}

# ----- iterate substrates and product types -----

types = c("cis","revCis","trans")
subU = DBQual$substrateID %>% unique()

Hall = lapply(types, function(t){
  Hsub = sapply(subU, function(s){
    
    xx = which(DBQual$substrateID == s & DBQual$spliceType == t)
    yy = which(rndDBQual$substrateID == s & rndDBQual$spliceType == t)
    # xx = which(DBQual$spliceType == t)
    if (length(xx) > 0) {
      getEntropy(DB = DBQual[xx,], rndDB = rndDBQual[yy,]) 
    }
  })
  HsubDF = plyr::ldply(Hsub)
  names(HsubDF)[1] = "substrateID"
  # HsubDF = t(Hsub) %>% as.data.frame()
  
  return(HsubDF)
})
names(Hall) = types
HallDF = plyr::ldply(Hall)
names(HallDF)[1] = "spliceType"

rem = which(rowSums(HallDF[,SRnames]) == 0)
HallDF = HallDF[-rem,]

HDF = HallDF %>% tidyr::gather(residue, entropy, -substrateID, -spliceType)
# HDF = HallDF %>% tidyr::gather(residue, entropy, -spliceType)
HDF = HDF %>%
  mutate(SR = ifelse(residue %in% c("P4","P3","P2","P1"), "SR1", "intervening"),
         SR = ifelse(residue %in% c("P4_","P3_","P2_","P1_"), "SR2", SR))


viop_all = ggplot(HDF, aes(x = factor(residue, levels = SRnames), y = entropy, fill = spliceType)) +
  geom_violin(draw_quantiles = c(0.5)) + 
  geom_vline(xintercept = c(4.5,8.5,12.5), col = "gray", lty = "dashed") +
  scale_fill_manual(values = c(plottingCols[["cis"]], plottingCols[["revCis"]], plottingCols[["trans"]])) +
  ggtitle("SR information content") +
  xlab("residue") + ylab("sequence entropy [bits]") +
  facet_wrap(~spliceType, scales = "free")

ggsave("results/SR1vsSR2/entropy/_residue-wise_entropy.png", plot = viop_all,
       height = 5, width = 12, dpi = "retina")

viop_srs = ggplot(HDF, aes(x = SR, y = entropy, fill = spliceType)) +
  geom_violin(draw_quantiles = c(0.5)) + 
  scale_fill_manual(values = c(plottingCols[["cis"]], plottingCols[["revCis"]], plottingCols[["trans"]])) +
  ggtitle("SR information content") +
    xlab("residue") + ylab("sequence entropy [bits]") +
  facet_wrap(~spliceType, scales = "free")

ggsave("results/SR1vsSR2/entropy/_SR-wise_entropy.png", plot = viop_srs,
       height = 5, width = 8, dpi = "retina")

