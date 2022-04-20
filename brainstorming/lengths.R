### PCPS brainstorming ###
# description:  SR and intervening sequence lengths
# input:        IDP protein data set
# output:       SR and intervening sequence length distributions
# author:       HPR

library(dplyr)
library(stringr)
library(ggplot2)
source("../../proteinsPCPS/new/src/invitroSPI_utils.R")
source("../../proteinsPCPS/new/src/plotting_utils.R")

theme_set(theme_classic())

### INPUT ###
load("../../proteinsPCPS/new/data/aSPIre.RData")
polypeps = read.csv("../../invitroSPI/revision/results/ProteasomeDB/ProteasomeDB.csv", stringsAsFactors = F)


### MAIN PART ###
# ----- data preprocessing -----
proteins = Kinetics

proteins = proteins %>%
  distinct(substrateID, pepSeq, .keep_all = T) %>%
  tidyr::separate_rows(spectralAngles, assignedScans, sep=";") %>%
  ILredundancy() %>%
  disentangleMultimappers.Type() %>%
  disentangleMultimappers.SRlen() %>%
  disentangleMultimappers.IVSeqLen()

polypeps = polypeps %>%
  ILredundancy() %>%
  filterPepLength() %>%
  disentangleMultimappers.Type() %>%
  disentangleMultimappers.SRlen() %>%
  disentangleMultimappers.IVSeqLen() %>%
  remSynthErrors() %>%
  filter20Sstandard() %>%
  uniquePeptides() %>%
  DBcosmetics()

# ----- splice-reactant length -----
SRLength = function(df, nm) {
  
  df = df[str_detect(df$productType, "PSP"), ] %>%
    removeMultimappers.Type() %>%
    removeMultimappers.SRlen()
  
  srlen = function(db) {
    pos = str_split_fixed(db$positions, pattern = "_", n = Inf)
    
    sr = data.frame(value_sr1 = as.numeric(pos[, 2]) - as.numeric(pos[, 1]) + 1,
                    value_sr2 = as.numeric(pos[, 4]) - as.numeric(pos[, 3]) + 1,
                    substrateID = db$substrateID,
                    type = db$spliceType)
    
    sr$type = factor(sr$type, levels = c("cis", "revCis", "trans"))
    return(sr)
  }
  
  sr = srlen(df)
  sr = sr %>% 
    tidyr::gather(value_sr1, value_sr2, -substrateID, -type) %>%
    rename(srtype = value_sr1,
           value = value_sr2)
  
  # plotting
  sr_plot = ggplot(data = sr, aes(x = type, y = value, fill = srtype)) +
    geom_split_violin(size = .25) +
    stat_summary(fun = median, fun.min = median, fun.max = median,
                 geom = "crossbar",width = .2,
                 position = position_dodge(width = .2)) +
    xlab("product type") +
    ylab("SR length (aa residues)") +
    ggtitle(nm) +
    scale_fill_manual("SR", values = c("gray","lightblue"), labels = c("SR1","SR2")) +
    theme(axis.text.x = element_text(angle = 90))
  
  ggsave(filename = paste0("results/SRlengths_",nm,".png"),
         plot = sr_plot, dpi = "retina", height = 4, width = 4)
}

SRLength(proteins, nm = "proteins")
SRLength(polypeps, nm = "polypeptides")


# ----- intervening sequence length -----
IVLength = function(df, nm) {
  
  df = df[str_detect(df$productType, "PSP"), ] %>%
    removeMultimappers.Type() %>%
    removeMultimappers.IVSeqLen()
  
  ivlen = function(db) {
    # intervening sequence length
    # calculated as in Specht et al., Paes et al. and SPI-ART
    pos = str_split_fixed(db$positions, pattern = "_", n = Inf)
    iv = rep(NA, nrow(db))
    
    cistrans.idx = which(db$spliceType %in% c("cis", "trans"))
    iv[cistrans.idx] = (abs(as.numeric(pos[cistrans.idx, 3]) - as.numeric(pos[cistrans.idx, 2])) - 1) %>%
      as.numeric()
    
    revcis.idx = which(db$spliceType == "revCis")
    iv[revcis.idx] = (abs(as.numeric(pos[revcis.idx, 1]) - as.numeric(pos[revcis.idx, 4])) - 1) %>%
      as.numeric()
    
    ivl = data.frame(value = iv,
                     type = db$spliceType,
                     substrateID = db$substrateID,
                     substrateSeq = db$substrateSeq)
    ivl$type = factor(ivl$type, levels = c("cis", "revCis", "trans"))
    return(ivl)
  }
  
  iv = ivlen(df)
  iv = iv %>% 
    mutate(rel = value/nchar(substrateSeq),
           placeholder = "x")
  
  # plotting
  iv_plot = ggplot(data = iv[iv$type != "trans", ], aes(x = placeholder,y = rel, fill = type)) +
    geom_split_violin(size = .25) +
    stat_summary(fun = median, fun.min = median, fun.max = median,
                 geom = "crossbar",width = .2,
                 position = position_dodge(width = .2)) +
    xlab("") +
    ylab("relative intervening sequence length") +
    ggtitle(nm) +
    scale_fill_manual("splice type", values = c(plottingCols["cis"],plottingCols["revCis"])) +
    theme(axis.text.x = element_text(angle = 90))
  
  ggsave(filename = paste0("results/IVlengths_",nm,".png"),
         plot = iv_plot, dpi = "retina", height = 4, width = 4)
}

IVLength(proteins, nm = "proteins")
IVLength(polypeps, nm = "polypeptides")


### OUTPUT ###

