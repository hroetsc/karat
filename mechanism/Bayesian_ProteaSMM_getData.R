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
  out = SCSandPSP_allSubs(Quant, target, meanOverBio = F, Zscale = F)
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
  X = X[,!colnames(X) %in% c("P1;X", "P1_;X")]
  
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
