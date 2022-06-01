### karat projetc - PCPS mechanism ###
# description:  compare number of theoretical vs. number of observed peptides over product lengths
# input:        invitroSPI: Roetschke SciData, EGFR ery, WT substrates --> spliced peptides only
# output:       fraction of observed peptides over length
# author:       HPR

library(dplyr)
library(stringr)
source("src/invitroSPI_utils.R")
source("src/plotting_utils.R")
source("../brainstorming/src/number-of-products.R")

### INPUT ###
load("data/invitroSPI.RData")

### MAIN PART ###
suppressWarnings(dir.create("results/SRspecificity/"))

# ----- preprocessing -----
DB = ProteasomeDB %>%
  ILredundancy() %>%
  disentangleMultimappers.Type() %>%
  removeMultimappers.Type() %>%
  uniquePeptides()

# ----- get lengths -----

getLen = function(DB) {
  
  DB$pepLen = nchar(DB$pepSeq)
  
  pos = str_split_fixed(DB$positions, pattern = "_", n = Inf)
  pspidx = which(DB$spliceType != "PCP")
  DB$SR1Len = NA
  DB$SR2Len = NA
  DB$SRshortLen = NA
  DB$SR1Len[pspidx] = as.numeric(pos[pspidx, 2]) - as.numeric(pos[pspidx, 1]) + 1
  DB$SR2Len[pspidx] = as.numeric(pos[pspidx, 4]) - as.numeric(pos[pspidx, 3]) + 1
  DB$SRshortLen[pspidx] = do.call(pmin, DB[pspidx,c("SR1Len", "SR2Len")])
  
  DB$IVlen = NA
  cistrans = which(DB$spliceType %in% c("cis","trans"))
  revcis = which(DB$spliceType == "revCis")
  DB$IVlen[c(cistrans,revcis)] = (abs(as.numeric(pos[c(cistrans,revcis), 3]) - as.numeric(pos[c(cistrans,revcis), 2])) - 1) %>%
    as.numeric()
  # DB$IVlen[revcis] = (abs(as.numeric(pos[revcis, 1]) - as.numeric(pos[revcis, 4])) - 1) %>%
  #   as.numeric()
  
  return(DB)
}

DB = getLen(DB)

