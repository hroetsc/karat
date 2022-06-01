### karat projetc - PCPS mechanism ###
# description:  get length distributions on quantitative level
# input:        aSPIre: EGFR ery, WT substrates
#               random database for quantitative 
# output:       peptide, SR and intervening sequence length distributions on quantitative level
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
load("data/aSPIre.RData")
load("data/randomDB_Quant_aSPIre.RData")


### MAIN PART ###
# ------ preprocessing -----
DB = Kinetics

DB = DB %>%
  tidyr::separate_rows(digestTimes, intensities, sep=";") %>%
  rename(digestTime = digestTimes,
         intensity = intensities) %>%
  mutate(intensity = as.numeric(intensity),
         digestTime = as.numeric(digestTime)) %>%
  filter(digestTime == 4) %>%
  ILredundancy() %>%
  disentangleMultimappers.Type() %>%
  disentangleMultimappers.SRlen() %>%
  disentangleMultimappers.IVSeqLen()

rndDB = randomQuant %>%
  ILredundancy() %>%
  disentangleMultimappers.Type() %>%
  disentangleMultimappers.SRlen() %>%
  disentangleMultimappers.IVSeqLen() %>%
  uniquePeptides()


# ----- functions for plotting -----
splittedViolinPlotQuant = function(data) {
  
  data$intensity[data$intensity == 0] = NA
  int = log10(data$intensity+1)
  scaled_intensity = round((100-1)*((int - min(int, na.rm = T))/(max(int, na.rm = T) - min(int, na.rm = T))) + 1)
  scaled_intensity[is.na(scaled_intensity)] = 0
  
  data$scaled_intensity = scaled_intensity
  data$db = "qualitative"
  
  y = data.frame(lapply(data, rep, data$scaled_intensity))
  y$db = "quantitative"
  
  DATA = rbind(data,y)
  DATA$db = factor(DATA$db, levels = c("qualitative","quantitative"))
  
  # print stats
  s = DATA %>%
    dplyr::group_by(type, db) %>%
    dplyr::summarise(n = dplyr::n(),
                     mean = mean(value),
                     median = median(value),
                     std = sd(value))
  print.data.frame(s)
  
  theme_set(theme_classic())
  
  # plotting
  sc = ggplot(data = DATA, aes(x = type, y = value, fill = db)) +
    geom_split_violin(size = .25) +
    stat_summary(fun = median, fun.min = median, fun.max = median,
                 geom = "crossbar", 
                 width = .2,
                 position = position_dodge(width = .2)) +
    xlab("product type") +
    scale_fill_manual("database",
                      values = c("olivedrab4", "darkolivegreen1"),
                      drop = F) +
    theme(axis.text.x = element_text(angle = 90))
  
  # add labs and title downstream
  
  return(sc)
  
}



# ----- peptide length -----
suppressWarnings(dir.create("results/length/"))

PepLength = function(trueDB, randomDB) {
  
  trueDB = removeMultimappers.Type(trueDB)
  randomDB = removeMultimappers.Type(randomDB)
  
  # get peptide length
  pl = data.frame(value = nchar(trueDB$pepSeq),
                  intensity = trueDB$intensity,
                  substrateID = trueDB$substrateID,
                  type = trueDB$spliceType,
                  db = "ProteasomeDB")
  pl = rbind(pl,
             data.frame(value = nchar(randomDB$pepSeq),
                        intensity = randomDB$intensity,
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
  
  # plot MS1 intensity distribution for each peptide length
  pepLenInt = ggplot(pl[which(pl$db == "ProteasomeDB"), ],
         aes(x = factor(value), y = log10(intensity+1), fill = type)) +
    geom_boxplot() + 
    ylim(c(-2, max(log10(pl$intensity+1)))) +
    xlab("peptide length") +
    ylab("intensity") +
    geom_text(aes(label = factor(n)), y=-1, size = 2, stat = "count") +
    scale_fill_manual(values = c(plottingCols[["PCP"]], plottingCols[["cis"]], plottingCols[["revCis"]], plottingCols[["trans"]])) +
    facet_wrap(~type, scales = "free", nrow = 4)
  
  # plot length weighted by intensity
  pepLenQuant = splittedViolinPlotQuant(data = pl[pl$db == "ProteasomeDB",]) +
    ylab("peptide product length (aa residues)") +
    ggtitle("peptide length distribution")
  
  # save plots
  ggsave(filename = "results/length/quant_PepLen.png", plot = pepLen, dpi = "retina", height = 6, width = 6)
  ggsave(filename = "results/length/quant_PepLenInt.png", plot = pepLenInt, dpi = "retina", height = 20, width = 10)
  ggsave(filename = "results/length/quant_PepLenQuant.png", plot = pepLenQuant, dpi = "retina", height = 6, width = 6)
}

PepLength(DB, rndDB)

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
                    intensity = db$intensity,
                    substrateID = db$substrateID,
                    type = db$spliceType,
                    db = db_name)
    sr$type = factor(sr$type, levels = c("cis", "revCis", "trans"))
    return(sr)
  }
  
  sr = rbind(srlen(db = trueDB, db_name = "ProteasomeDB"),
             srlen(db = randomDB, db_name = "randomDB"))
  
  sr$type = factor(sr$type, levels = c("PCP", "cis", "revCis", "trans"))
  
  # distributions of both SRs
  sr1 = sr[, c("value_sr1", "intensity", "type", "substrateID","db")]
  names(sr1)[1] = "value"
  sr2 = sr[, c("value_sr2", "intensity", "type", "substrateID", "db")]
  names(sr2)[1] = "value"
  
  # get counts
  sr1 = sr1 %>%
    group_by(type, db, value) %>%
    mutate(n = n()) %>%
    ungroup()
  sr2 = sr2 %>%
    group_by(type, db, value) %>%
    mutate(n = n()) %>%
    ungroup()
  
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
  
  # plot MS1 intensity distribution for each SR length
  sr1LenInt = ggplot(sr1[which(sr1$db == "ProteasomeDB"), ],
                     aes(x = factor(value), y = log10(intensity+1), fill = type)) +
    geom_boxplot() + 
    ylim(c(-2, max(log10(sr1$intensity+1)))) +
    xlab("SR1 length") +
    ylab("intensity") +
    geom_text(aes(label = factor(n)), y=-1, size = 2, stat = "count") +
    scale_fill_manual(values = c(plottingCols[["cis"]], plottingCols[["revCis"]], plottingCols[["trans"]])) +
    facet_wrap(~type, scales = "free", nrow = 4)
  
  sr2LenInt = ggplot(sr2[which(sr2$db == "ProteasomeDB"), ],
                     aes(x = factor(value), y = log10(intensity+1), fill = type)) +
    geom_boxplot() + 
    ylim(c(-2, max(log10(sr2$intensity+1)))) +
    xlab("SR2 length") +
    ylab("intensity") +
    geom_text(aes(label = factor(n)), y=-1, size = 2, stat = "count") +
    scale_fill_manual(values = c(plottingCols[["cis"]], plottingCols[["revCis"]], plottingCols[["trans"]])) +
    facet_wrap(~type, scales = "free", nrow = 4)
  
  srLenInt = gridExtra::grid.arrange(sr1LenInt,sr2LenInt,ncol=2)
  
  # plot length weighted by intensity
  print("SR1")
  SR1LenQuant = splittedViolinPlotQuant(data = sr1[sr1$db == "ProteasomeDB",]) +
    ylab("SR1 length (aa residues)") +
    ggtitle("SR1 length distribution")
  print("SR2")
  SR2LenQuant = splittedViolinPlotQuant(data = sr2[sr2$db == "ProteasomeDB",]) +
    ylab("SR2 length (aa residues)") +
    ggtitle("SR2 length distribution")
  
  srLenQuant = gridExtra::grid.arrange(SR1LenQuant,SR2LenQuant,ncol=2)
  
  # save plots
  ggsave(filename = "results/length/quant_SRLen.png", plot = srLen, dpi = "retina", height = 6, width = 12)
  ggsave(filename = "results/length/quant_SRLenInt.png", plot = srLenInt, dpi = "retina", height = 15, width = 20)
  ggsave(filename = "results/length/quant_SRLenQuant.png", plot = srLenQuant, dpi = "retina", height = 6, width = 12)
}

SRLength(DB, rndDB)

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
                     intensity = db$intensity,
                     type = db$spliceType,
                     substrateID = db$substrateID,
                     db = db_name)
    return(ivl)
  }
  
  iv = rbind(ivlen(db = trueDB, db_name = "ProteasomeDB"),
             ivlen(db = randomDB, db_name = "randomDB"))
  
  # get counts
  iv = iv %>%
    group_by(type, db, value) %>%
    mutate(n = n()) %>%
    ungroup()
  
  iv$type = factor(iv$type, levels = c("cis", "revCis", "trans"))
  
  # plot qualitative DB vs random
  ivLen = splittedViolinPlot(data = iv) +
    ylab("intervening sequence length (aa residues)") +
    ggtitle("intervening sequence length distribution")
  
  # plot MS1 intensity distribution for each intervening sequence length
  ivLenInt = ggplot(iv[which(iv$db == "ProteasomeDB" & iv$type %in% c("cis","revCis")), ],
                     aes(x = factor(value), y = log10(intensity), fill = type)) +
    geom_boxplot() + 
    ylim(c(-2, max(log10(iv$intensity+1)))) +
    xlab("intervening sequence length") +
    ylab("intensity") +
    geom_text(aes(label = factor(n)), y=-1, size = 2, stat = "count") +
    scale_fill_manual(values = c(plottingCols[["cis"]], plottingCols[["revCis"]])) +
    facet_wrap(~type, scales = "free", nrow = 2)
  
  # plot length weighted by intensity
  ivLenQuant = splittedViolinPlotQuant(data = iv[iv$db == "ProteasomeDB",]) +
    ylab("intervening sequence length (aa residues)") +
    ggtitle("intervening sequence length distribution")
  
  # save plots
  ggsave(filename = "results/length/quant_IVLen.png", plot = ivLen, dpi = "retina", height = 6, width = 6)
  ggsave(filename = "results/length/quant_IVLenInt.png", plot = ivLenInt, dpi = "retina", height = 10, width = 10)
  ggsave(filename = "results/length/quant_IVLenQuant.png", plot = ivLenQuant, dpi = "retina", height = 6, width = 6)
  
}

IVSeqLength(DB, rndDB)

