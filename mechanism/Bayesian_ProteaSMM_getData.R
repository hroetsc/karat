### karat projetc - PCPS mechanism ###
# description:  get input data for Bayesian ProteaSMM
# input:        aSPIre: Roetschke SciData, EGFR ery, WT substrates (quantitative DB)
# output:       substrate matrices for each residue in each substrate + SCS/PSP-P1(')
# author:       HPR

library(dplyr)
library(stringr)
library(pheatmap)
source("src/invitroSPI_utils.R")
source("src/_extract-aa.R")
source("src/SCS+PSP-P1.R")

suppressWarnings(dir.create("data/ProteaSMM/"))

AAchar_here = c("P","G","C","M","A","V","I","L","F","Y","W","H","R","K","D","E","N","Q","S","T","X")
AAchar_here_sorted = sort(AAchar_here)

### INPUT ###
load("data/aSPIre.RData")
polypeps = Kinetics
load("../../proteinsPCPS/new/data/aSPIre.RData")
proteins = Kinetics


### MAIN PART ###
# ----- preprocessing -----
nm = intersect(names(polypeps), names(proteins))
# rbind(polypeps[,nm] %>% mutate(dataset = "polypeps"),
#       proteins[,nm] %>% mutate(dataset = "proteins"))
Quant = polypeps %>% mutate(dataset = "polypeps") %>%
  # ILredundancy() %>%
  disentangleMultimappers.Type() %>%
  tidyr::separate_rows(digestTimes, intensities, sep=";") %>%
  rename(digestTime = digestTimes,
         intensity = intensities) %>%
  mutate(intensity = as.numeric(intensity),
         digestTime = as.numeric(digestTime)) %>%
  filter((dataset == "polypeps" & digestTime == 4) | (dataset == "proteins" & digestTime %in% c(3,4))) %>%
  resolve_multimapper() %>%
  tidyr::separate_rows(positions, sep = ";")

# ----- get binary matrices for substrates -----
# for only one SR
getSubstrateCounts = function(DB, SRpos, allCombos) {
  subSeqs = DB$substrateID %>% unique()
  
  allFeatures = lapply(c(1:length(subSeqs)), function(i){
    
    S = DB$substrateSeq[DB$substrateID == subSeqs[i]][1]
    k = which(DB$substrateID == subSeqs[i])
    
    pos = data.frame(substrateSeq = S, residue = c(1:nchar(S)))
    SRTBL = sapply(SRpos, function(x){
      substr(pos$substrateSeq, start = pos$residue+x, stop = pos$residue+x)
    }) %>%
      as.data.frame()
    SRTBL[SRTBL == ""] = "X"
    
    master = matrix(0, nrow = nrow(SRTBL), ncol = length(allCombos))
    colnames(master) = allCombos
    
    for (j in 1:nrow(SRTBL)) {
      cntN = paste(colnames(SRTBL),SRTBL[j,],sep = ";")
      master[j,cntN] = 1
    }
    
    return(cbind(c(1:nchar(S)),master) %>% as.data.frame())
  })
  names(allFeatures) = subSeqs
  
  return(allFeatures)
}



# for combinations of SR1 and SR2
getSubstrateCountsForSTS = function(DB, SRpos, reactant, allCombos) {
  
  subSeqs = DB$substrateID %>% unique()
  
  allFeatures = lapply(c(1:length(subSeqs)), function(i){
    
    S = DB$substrateSeq[DB$substrateID == subSeqs[i]][1]
    
    residues = tidyr::crossing(c(1:nchar(S)), c(1:nchar(S)), .name_repair = "minimal") %>%
      as.data.frame()
    names(residues) = c("P1","P1_")
    
    pos = data.frame(substrateSeq = S, P1 = residues$P1, P1_ = residues$P1_)
    
    SR1TBL = sapply(SRpos[reactant == "SR1"], function(x){
      substr(pos$substrateSeq, start = pos$P1+x, stop = pos$P1+x)
    }) %>%
      as.data.frame()
    SR1TBL[SR1TBL == ""] = "X"
    
    SR2TBL = sapply(SRpos[reactant == "SR2"], function(x){
      substr(pos$substrateSeq, start = pos$P1_+x, stop = pos$P1_+x)
    }) %>%
      as.data.frame()
    SR2TBL[SR2TBL == ""] = "X"
    
    SRTBL = cbind(SR1TBL, SR2TBL)
    
    master = matrix(0, nrow = nrow(SRTBL), ncol = length(allCombos))
    colnames(master) = allCombos
    
    for (j in 1:nrow(SRTBL)) {
      cntN = paste(colnames(SRTBL),SRTBL[j,],sep = ";")
      master[j,cntN] = 1
    }
    
    return(cbind(residues,master) %>% as.data.frame())
  })
  names(allFeatures) = subSeqs
  
  return(allFeatures)
}


# ----- calculate STS -----

getDataSTS = function(nm) {
  fname = paste0("data/ProteaSMM/",nm,"_STS/")
  print(fname)
  suppressWarnings(dir.create(fname))
  
  # splice-reactants
  SRpos = c("P4"=-3, "P3"=-2, "P2"=-1, "P1"=0,
            "P-1"=1, "P-2"=2, "P-3"=3, "P-4"=4,
            "P-4_"=-4, "P-3_"=-3, "P-2_"=-2, "P-1_"=-1,
            "P1_"=0, "P2_"=1, "P3_"=2, "P4_"=3)
  # SRpos = c("P6"=-5,"P5"=-4,"P4"=-3, "P3"=-2, "P2"=-1, "P1"=0,
  #           "P-1"=1, "P-2"=2, "P-3"=3, "P-4"=4, "P-5"=5, "P-6"=6,
  #           "P-6_"=-6,"P-5_"=-5,"P-4_"=-4, "P-3_"=-3, "P-2_"=-2, "P-1_"=-1,
  #           "P1_"=0, "P2_"=1, "P3_"=2, "P4_"=3, "P5_"=4, "P6_"=5)
  # SRpos = c("P8"=-7,"P7"=-6,"P6"=-5,"P5"=-4,"P4"=-3, "P3"=-2, "P2"=-1, "P1"=0,
  #           "P-1"=1, "P-2"=2, "P-3"=3, "P-4"=4, "P-5"=5, "P-6"=6, "P-7"=7, "P-8"=8,
  #           "P-8_"=-8,"P-7_"=-7, "P-6_"=-6,"P-5_"=-5,"P-4_"=-4, "P-3_"=-3, "P-2_"=-2, "P-1_"=-1,
  #           "P1_"=0, "P2_"=1, "P3_"=2, "P4_"=3, "P5_"=4, "P6_"=5, "P7_"=6, "P8_"=7)
  reactant = c(rep("SR1", 8), rep("SR2", 8))
  
  # get colnames
  interesting_residues = names(SRpos)
  allCombos = tidyr::crossing(interesting_residues,AAchar_here_sorted)
  allCombos = do.call(paste, c(allCombos[c("interesting_residues","AAchar_here_sorted")], sep=";"))
  
  # get SCS and PSP for P1 or P1' for each bio rep
  out = STS_allSubs(Quant, meanOverBio = F)
  yPredDF = plyr::ldply(out) %>%
    as.data.frame()
  names(yPredDF)[1] = "substrateID"
  
  yPredDF_filter = yPredDF
  paste0(nrow(yPredDF_filter), " residues passed filtering") %>% print()
  
  # get binary substrate counts
  substrateMatrix = getSubstrateCountsForSTS(DB = Quant, SRpos, reactant, allCombos)
  substrateCounts = plyr::ldply(substrateMatrix)
  names(substrateCounts)[1:3] = c("substrateID","P1","P1_")
  
  # merge with SCS and PSP-P1
  # add prediction target
  ALL = left_join(substrateCounts, yPredDF_filter %>% dplyr::select(substrateID, P1,P1_, sts_mean)) %>%
    rename(y = sts_mean) %>%
    na.omit()
  
  # separate by bio reps
  t = str_split_fixed(ALL$y, ";", Inf) %>% as.matrix()
  t = apply(t,2,as.numeric)
  
  X = ALL[,allCombos] %>% as.matrix()
  X = X[,!colnames(X) %in% c("P1;X", "P1_;X")]
  
  # merge and return
  DATA = list(X = X,
              t = t,
              substrateIDs = ALL$substrateID)  # !!!!!
  
  save(DATA, file = paste0(fname, "DATA.RData"))
}

getDataSTS(nm = "PSP")


# ----- calculate SCS/PSP-P1 -----

getData = function(target, SRpos, col, nCol, nm, features_from) {
  fname = paste0("data/ProteaSMM/",nm,"_",features_from,"feat_",target,"/")
  print(fname)
  suppressWarnings(dir.create(fname))
  
  # get colnames
  interesting_residues = names(SRpos)
  allCombos = tidyr::crossing(interesting_residues,AAchar_here_sorted)
  allCombos = do.call(paste, c(allCombos[c("interesting_residues","AAchar_here_sorted")], sep=";"))
  
  # get SCS and PSP for P1 or P1' for each bio rep
  out = SCSandPSP_allSubs(Quant, target, meanOverBio = F, Zscale = F, SR2forSCS = T)
  yPredDF = plyr::ldply(out) %>%
    as.data.frame()
  yPredDF$target = target
  names(yPredDF)[1] = "substrateID"
  
  # keep only residues with sufficient number of peptides detected
  # yPredDF_filter = yPredDF %>%
  #   rename(nnum = nCol) %>%
  #   filter(nnum > 0) %>%
  #   group_by(substrateID) %>%
  #   filter(nnum > quantile(nnum, .05)) %>% ungroup()
  yPredDF_filter = yPredDF
  
  paste0(nrow(yPredDF_filter), " residues passed filtering") %>% print()
  
  # get binary substrate counts
  substrateMatrix = getSubstrateCounts(Quant, SRpos, allCombos)
  substrateCounts = plyr::ldply(substrateMatrix)
  names(substrateCounts)[1:2] = c("substrateID","residue")
  
  # merge with SCS and PSP-P1
  # add prediction target
  ALL = left_join(substrateCounts, yPredDF_filter %>% dplyr::select(substrateID, residue, col)) %>%
    rename(y = col) %>%
    na.omit()
  
  # separate by bio reps
  t = str_split_fixed(ALL$y, ";", Inf) %>% as.matrix()
  t = apply(t,2,as.numeric)
  
  X = ALL[,allCombos] %>% as.matrix()
  if (col == "scs_mean") {
    X = X[,!colnames(X) %in% c("P1;X", "P1_;X", "P-1;X")]
  } else {
    X = X[,!colnames(X) %in% c("P1;X", "P1_;X")]
  }
  
  
  # keep only the residues that have peptides detected
  # rem = which(rowSums(t, na.rm = T) == 0)
  # if (length(rem) > 0) {
  #   paste0("removing ", length(rem), " empty rows") %>% print()
  #   t = t[-rem,]
  #   X = X[-rem,]
  # }
  
  pdf(file = paste0(fname, "DATA.pdf"), width = 50, height = 50)
  pheatmap(cbind(X, t/100), cluster_cols = F) %>% print()
  dev.off()
  
  # t[-grep("MM",ALL$substrateID),] %>% as.vector() %>% log10() %>% na.omit() %>% density() %>% plot(ylim = c(0,.5))
  # t[grep("MM",ALL$substrateID),] %>% as.vector() %>% log10() %>% na.omit() %>% density() %>% lines(col = "red")
  
  # merge and return
  DATA = list(X = X,
              t = t,
              substrateIDs = ALL$substrateID)
  
  save(DATA, file = paste0(fname, "DATA.RData"))
}

# ----- just 4 positions -----
# PSP: SR1
getData(target = "P1",
        features_from = "SR1",
        SRpos = c("P4"=-3, "P3"=-2, "P2"=-1, "P1"=0,
                  "P-1"=1, "P-2"=2, "P-3"=3, "P-4"=4),
        col = "psp_mean", nCol = "psp_n",
        nm = "PSP")

# PSP: SR2
getData(target = "P1_",
        features_from = "SR2",
        SRpos = c("P-4_"=-4, "P-3_"=-3, "P-2_"=-2, "P-1_"=-1,
                  "P1_"=0, "P2_"=1, "P3_"=2, "P4_"=3),
        col = "psp_mean", nCol = "psp_n",
        nm = "PSP")

# PCP: precursor
getData(target = "P1",
        features_from = "SR1",
        SRpos = c("P4"=-3, "P3"=-2, "P2"=-1, "P1"=0,
                  "P-1"=1, "P-2"=2, "P-3"=3, "P-4"=4),
        col = "scs_mean", nCol = "scs_n",
        nm = "PCP")


# PSP-P1: SR1 + SR2
getData(target = "P1",
        features_from = "SR1+2",
        SRpos = c("P4"=-3, "P3"=-2, "P2"=-1, "P1"=0,
                  "P-1"=1, "P-2"=2, "P-3"=3, "P-4"=4,
                  "P-4_"=-4, "P-3_"=-3, "P-2_"=-2, "P-1_"=-1,
                  "P1_"=0, "P2_"=1, "P3_"=2, "P4_"=3),
        col = "psp_mean", nCol = "psp_n",
        nm = "PSP")

# ----- extended positions -----
getData(target = "P1",
        features_from = "SR1ext",
        SRpos = c("P6"=-5,"P5"=-4,"P4"=-3, "P3"=-2, "P2"=-1, "P1"=0,
                  "P-1"=1, "P-2"=2, "P-3"=3, "P-4"=4, "P-5"=5, "P-6"=6),
        col = "psp_mean", nCol = "psp_n",
        nm = "PSP")

# PSP: SR2
getData(target = "P1_",
        features_from = "SR2ext",
        SRpos = c("P-6_"=-6,"P-5_"=-5,"P-4_"=-4, "P-3_"=-3, "P-2_"=-2, "P-1_"=-1,
                  "P1_"=0, "P2_"=1, "P3_"=2, "P4_"=3, "P5_"=4, "P6_"=5),
        col = "psp_mean", nCol = "psp_n",
        nm = "PSP")

# PCP: precursor
getData(target = "P1",
        features_from = "SR1ext",
        SRpos = c("P6"=-5,"P5"=-4,"P4"=-3, "P3"=-2, "P2"=-1, "P1"=0,
                  "P-1"=1, "P-2"=2, "P-3"=3, "P-4"=4, "P-5"=5, "P-6"=6),
        col = "scs_mean", nCol = "scs_n",
        nm = "PCP")


# ----- generate in silico data set -----

# SCS
rtruncnorm <- function(n, mu, sigma, low, high) {
  # find quantiles that correspond the the given low and high levels.
  p_low <- pnorm(low, mu, sigma)
  p_high <- pnorm(high, mu, sigma)
  
  # draw quantiles uniformly between the limits and pass these
  # to the relevant quantile function.
  qnorm(runif(n, p_low, p_high), mu, sigma)
}

load("data/ProteaSMM/PCP_SR1extfeat_P1/DATA.RData")
X = DATA$X

set.seed(42)
w = rtruncnorm(n = ncol(X), mu = 0.6, sigma = 1, low = 0, high = 2)
tsim = X%*%log(w)

DATA = list(X = X,
            t = tsim,
            w = w,
            substrateIDs = DATA$substrateIDs)

dir.create("data/ProteaSMM/insilico/")
save(DATA, file = "data/ProteaSMM/insilico/DATA.RData")



# STS
load("data/ProteaSMM/PSP_STS/DATA.RData")
X = DATA$X

set.seed(42)
w = rtruncnorm(n = ncol(X), mu = 0.3, sigma = 1, low = 0, high = 1.3)
tsim = X%*%log(w)

DATA = list(X = X,
            t = tsim,
            w = w,
            substrateIDs = DATA$substrateIDs)

dir.create("data/ProteaSMM/STS_insilico/")
save(DATA, file = "data/ProteaSMM/STS_insilico/DATA.RData")

