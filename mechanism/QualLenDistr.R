### karat projetc - PCPS mechanism ###
# description:  get length distributions on qualitative level
# input:        invitroSPI: Roetschke SciData, EGFR ery, WT substrates
#               random database for qualitative DB 
# output:       peptide, SR and intervening sequence length distributions on qualitative level
# author:       HPR

library(dplyr)
library(stringr)
library(seqinr)
library(ggplot2)
library(twosamples)
library(dgof)
source("src/invitroSPI_utils.R")
source("src/plotting_utils.R")

theme_set(theme_classic())


### INPUT ###
# load("data/invitroSPI.RData")
# load("data/randomDB_smart.RData")

load("../../proteinsPCPS/new/data/aSPIre.RData")
load("../../proteinsPCPS/new/data/randomDB.RData")

### MAIN PART ###
suppressWarnings(dir.create("results/length/"))

# ------ preprocessing -----
DB = Kinetics %>%
  ILredundancy() %>%
  disentangleMultimappers.Type() %>%
  disentangleMultimappers.SRlen() %>%
  disentangleMultimappers.IVSeqLen() %>%
  uniquePeptides()

RANDOM = rndDB %>%
  ILredundancy() %>%
  disentangleMultimappers.Type() %>%
  disentangleMultimappers.SRlen() %>%
  disentangleMultimappers.IVSeqLen() %>%
  uniquePeptides()

save(RANDOM, file = "data/randomDB_disentangled+unique.RData")


# ----- peptide length -----
trueDB = DB
randomDB = RANDOM

PepLength = function(trueDB, randomDB) {
  
  trueDB = removeMultimappers.Type(trueDB)
  randomDB = removeMultimappers.Type(randomDB)
  
  # get peptide length
  pl = data.frame(value = nchar(trueDB$pepSeq),
                  substrateID = trueDB$substrateID,
                  type = trueDB$spliceType,
                  db = "ProteasomeDB")
  pl = rbind(pl,
             data.frame(value = nchar(randomDB$pepSeq),
                        substrateID = randomDB$substrateID,
                        type = randomDB$spliceType,
                        db = "randomDB"))
  
  # get counts
  pl = pl %>%
    group_by(type, db, value) %>%
    mutate(n = n()) %>%
    ungroup()
  
  pl$type = factor(pl$type, levels = c("PCP", "cis", "revCis", "trans"))
  
  # plot qualitative DB vs random
  pepLen = splittedViolinPlot(data = pl) +
    ylab("peptide product length (aa residues)") +
    ggtitle("peptide length distribution")
  
  # save plots
  ggsave(filename = "results/length/qual_PepLen.png", plot = pepLen, dpi = "retina", height = 6, width = 6)
}

PepLength(DB, RANDOM)


# ----- SR length ------
SRLength = function(trueDB, randomDB) {
  
  trueDB = trueDB[which(trueDB$spliceType != "PCP"), ] %>%
    removeMultimappers.Type() %>%
    removeMultimappers.SRlen()
  
  randomDB = randomDB[which(randomDB$spliceType != "PCP"), ] %>%
    removeMultimappers.Type() %>%
    removeMultimappers.SRlen()
  
  # get SR length
  srlen = function(db, db_name) {
    pos = str_split_fixed(db$positions, pattern = "_", n = Inf)
    
    sr = data.frame(value_sr1 = as.numeric(pos[, 2]) - as.numeric(pos[, 1]) + 1,
                    value_sr2 = as.numeric(pos[, 4]) - as.numeric(pos[, 3]) + 1,
                    substrateID = db$substrateID,
                    type = db$spliceType,
                    db = db_name)
    sr$type = factor(sr$type, levels = c("cis", "revCis", "trans"))
    return(sr)
  }
  
  sr = rbind(srlen(db = trueDB, db_name = "ProteasomeDB"),
             srlen(db = randomDB, db_name = "randomDB")) %>%
    na.omit()
  
  sr$type = factor(sr$type, levels = c("PCP", "cis", "revCis", "trans"))
  
  # distributions of both SRs
  sr1 = sr[, c("value_sr1", "type", "substrateID","db")]
  names(sr1)[1] = "value"
  sr2 = sr[, c("value_sr2", "type", "substrateID", "db")]
  names(sr2)[1] = "value"
  
  # plot qualitative DB vs random
  print("SR1")
  sr1Len = splittedViolinPlot(data = sr1) +
    ylab("splice reactant length (aa residues)") +
    ggtitle("", subtitle = "SR1 length distribution")
  print("SR2")
  sr2Len = splittedViolinPlot(data = sr2) +
    ylab("splice reactant length (aa residues)") +
    ggtitle("", subtitle = "SR2 length distribution")
  srLen = gridExtra::grid.arrange(sr1Len, sr2Len, ncol=2)
  
  # SR1 vs SR2
  sr = srlen(db = trueDB, db_name = "ProteasomeDB") %>% 
    select(-db) %>%
    tidyr::gather(value_sr1, value_sr2, -substrateID, -type) %>%
    rename(srtype = value_sr1,
           value = value_sr2)
  sr$type = factor(sr$type, levels = c("cis", "revCis", "trans"))
  
  sr_plot = ggplot(data = sr, aes(x = type, y = value, fill = srtype)) +
    geom_split_violin(size = .25) +
    stat_summary(fun = median, fun.min = median, fun.max = median,
                 geom = "crossbar",width = .2,
                 position = position_dodge(width = .2)) +
    xlab("product type") +
    ylab("SR length (aa residues)") +
    ggtitle("SR length distribution") +
    scale_fill_manual("SR", values = c("gray","lightblue"), labels = c("SR1","SR2")) +
    theme(axis.text.x = element_text(angle = 90))
  
  # save plots
  ggsave(filename = "results/length/qual_SRLen.png", plot = srLen, dpi = "retina", height = 6, width = 12)
  ggsave(filename = "results/length/qual_SR1vsSR2.png", plot = sr_plot, dpi = "retina", height = 6, width = 6)
}

SRLength(DB, RANDOM)

# ----- intervening sequence length -----
IVSeqLength = function(trueDB, randomDB) {
  
  trueDB = trueDB[which(trueDB$spliceType != "PCP"), ] %>%
    removeMultimappers.Type() %>%
    removeMultimappers.IVSeqLen()
  
  randomDB = randomDB[which(randomDB$spliceType != "PCP"), ] %>%
    removeMultimappers.Type() %>%
    removeMultimappers.IVSeqLen()
  
  # get intervening sequence length
  ivlen = function(db, db_name) {
    
    # intervening sequence length
    # calculated as in Specht et al., Paes et al. and SPI-ART
    pos = str_split_fixed(db$positions, pattern = "_", n = Inf)
    iv = (abs(as.numeric(pos[, 3]) - as.numeric(pos[, 2])) - 1) %>%
      as.numeric()
    
    ivl = data.frame(value = iv,
                     type = db$spliceType,
                     substrateID = db$substrateID,
                     db = db_name)
    return(ivl)
  }
  
  iv = rbind(ivlen(db = trueDB, db_name = "ProteasomeDB"),
             ivlen(db = randomDB, db_name = "randomDB")) %>%
    na.omit()
  iv$type = factor(iv$type, levels = c("cis", "revCis", "trans"))
  
  # plot qualitative DB vs random
  ivLen = splittedViolinPlot(data = iv) +
    ylab("intervening sequence length (aa residues)") +
    ggtitle("intervening sequence length distribution")
  
  # save plots
  ggsave(filename = "results/length/qual_IVLen.png", plot = ivLen, dpi = "retina", height = 6, width = 6)
  
}

IVSeqLength(DB, RANDOM)


