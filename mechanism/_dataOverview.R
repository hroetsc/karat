### karat projetc - PCPS mechanismS ###
# description:  peptide lengths and product type frequencies
# input:        all data sets
# output:       peptide length distributions and product type frequencies of
#               non-spliced, cis-spliced and trans-spliced peptides
# author:       HPR

library(stringr)
library(dplyr)
library(ggplot2)
source("src/invitroSPI_utils.R")

theme_set(theme_classic())


### INPUT ###
load("data/aSPIre.RData")
wt_polypeps = Kinetics
load("data/MutPairs_aSPIre.RData")
mut_polypeps = Kinetics
load("../../proteinsPCPS/new/data/aSPIre.RData")
proteins = Kinetics

qual_polypeps = read.csv("../../invitroSPI/revision/results/ProteasomeDB/ProteasomeDB.csv", stringsAsFactors = F) %>%
  ILredundancy() %>%
  remSynthErrors() %>%
  filterEarlyTimepoints() %>%
  filter20Sstandard() 

### MAIN PART ###
suppressWarnings(dir.create("results/dataOverview/"))

# ----- data preprocessing -----
preprocessing = function(DB) {
   DB = DB %>%
     ILredundancy() %>%
     filterPepLength() %>%
     disentangleMultimappers.Type() %>%
     uniquePeptides()
   
   DB$spliceType[is.na(DB$spliceType)] = "PCP"
   DB$spliceType[DB$spliceType %in% c("cis","revCis","type_multi-mapper")] = "allcis"
   
   print(table(DB$spliceType))
   print(table(DB$productType))
   
   return(DB)
}

qual_polypeps = preprocessing(qual_polypeps)
wt_polypeps = preprocessing(wt_polypeps)
mut_polypeps = preprocessing(mut_polypeps)
proteins = preprocessing(proteins)

# ----- sequence numbers -----


# ----- splice type frequencies -----

nm = intersect(names(proteins), intersect(names(qual_polypeps), intersect(names(mut_polypeps), names(wt_polypeps))))

splicetypes = function(DB, nm) {
  
  FREQ = DB %>%
    group_by(substrateID, spliceType) %>%
    summarise(n = n()) %>%
    tidyr::spread(spliceType,n, fill = 0) %>%
    tidyr::gather(spliceType,n, -substrateID) %>%
    ungroup() %>% group_by(substrateID) %>%
    mutate(freq = n/sum(n))
  FREQ$spliceType = factor(FREQ$spliceType, levels = c("PCP","allcis","trans"))
  
  p = ggplot(FREQ, aes(x = spliceType, y = freq*100, fill = spliceType)) +
    geom_boxplot() +
    scale_fill_manual(values = c(plottingCols[["PCP"]], plottingCols[["allcis"]], plottingCols[["trans"]])) +
    ylim(c(0,100)) +
    xlab("") + ylab("")
  
  ggsave(filename = paste0("results/dataOverview/typeFreq_",nm,".png"), height = 4, width = 3, dpi = "retina")
}

splicetypes(DB = qual_polypeps, nm = "qual_polypeps")
splicetypes(DB = wt_polypeps, nm = "wt_polypeps")
splicetypes(DB = mut_polypeps, nm = "mut_polypeps")
splicetypes(DB = proteins, nm = "proteins")

# ----- peptide length distributions -----

peplengths = function(DB, nm) {
  
  DB$pepLen = nchar(DB$pepSeq)
  DB$spliceType = factor(DB$spliceType, levels = c("PCP","allcis","trans"))
  
  p = ggplot(DB, aes(x = spliceType, y = pepLen, fill = spliceType)) +
    geom_violin(draw_quantiles = 0.5) +
    scale_fill_manual(values = c(plottingCols[["PCP"]], plottingCols[["allcis"]], plottingCols[["trans"]])) +
    xlab("") + ylab("")
  
  ggsave(filename = paste0("results/dataOverview/pepLen_",nm,".png"), height = 4, width = 3, dpi = "retina")
}

peplengths(DB = qual_polypeps, nm = "qual_polypeps")
peplengths(DB = wt_polypeps, nm = "wt_polypeps")
peplengths(DB = mut_polypeps, nm = "mut_polypeps")
peplengths(DB = proteins, nm = "proteins")


