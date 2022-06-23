### karat projetc - PCPS mechanism ###
# description:  prediction of SCS/PSP-P1 based on posterior parameter distribution
# input:        BayesianTools: posterior.RData
#               data of prediction target
# output:       performance of model on validation data set
# author:       HPR

library(dplyr)
library(stringr)
source("src/invitroSPI_utils.R")
source("src/_extract-aa.R")
source("src/SCS+PSP-P1.R")

AAchar_here = c("P","G","C","M","A","V","I","L","F","Y","W","H","R","K","D","E","N","Q","S","T","X")
AAchar_here_sorted = sort(AAchar_here)

### INPUT ###
# chains
load("results/Bayesian_ProteaSMM/CHAINS/LOVpcp_chain.RData")
load("results/Bayesian_ProteaSMM/CHAINS/LOVsr1_chain.RData")
load("results/Bayesian_ProteaSMM/CHAINS/LOVsr2_chain.RData")

load("results/Analytical_ProteaSMM/PCP_SR1extfeat_P1/SMManalytical_PCP_LOV_P1.RData")
pcp_chain = out$allParams
load("results/Analytical_ProteaSMM/PSP_SR1extfeat_P1/SMManalytical_PSP_LOV_P1.RData")
sr1_chain = out$allParams
load("results/Analytical_ProteaSMM/PSP_SR2extfeat_P1_/SMManalytical_PSP_LOV_P1_.RData")
sr2_chain = out$allParams

# prediction target
# load("data/MutPairs_aSPIre.RData")
# Target = Kinetics[Kinetics$substrateID == "MM500", ]

# load("../../proteinsPCPS/new/data/aSPIre.RData")
# Target = Kinetics[Kinetics$substrateID == "P37840aSyn", ]
# Target = Kinetics[Kinetics$substrateID == "hTau", ]

load("data/aSPIre.RData")
Target = Kinetics[Kinetics$substrateID == "MM136", ]

### MAIN PART ###
# ----- preprocessing of prediction target -----
Target = Target %>%
  # ILredundancy() %>%  # !!!!!
  disentangleMultimappers.Type() %>%
  tidyr::separate_rows(digestTimes, intensities, sep=";") %>%
  rename(digestTime = digestTimes,
         intensity = intensities) %>%
  mutate(intensity = as.numeric(intensity),
         digestTime = as.numeric(digestTime)) %>%
  filter(digestTime %in% c(3,4)) %>%  # !!!!!
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

getData = function(Quant, target, SRpos, col, nCol, features_from) {
  
  # get colnames
  interesting_residues = names(SRpos)
  allCombos = tidyr::crossing(interesting_residues,AAchar_here_sorted)
  allCombos = do.call(paste, c(allCombos[c("interesting_residues","AAchar_here_sorted")], sep=";"))
  
  # get SCS and PSP for P1 or P1' for each bio rep
  out = SCSandPSP_allSubs(Quant, target, meanOverBio = F, Zscale = F, rawVals = F)
  yPredDF = plyr::ldply(out) %>%
    as.data.frame()
  yPredDF$target = target
  names(yPredDF)[1] = "substrateID"
  
  # plot distributions
  p = data.frame(residue = rep(yPredDF$residue, 2),
                 n = c(yPredDF$scs_n, yPredDF$psp_n),
                 p1 = c(rowMeans(apply(str_split_fixed(yPredDF$scs_mean,";",Inf),2,as.numeric)),
                        -1*rowMeans(apply(str_split_fixed(yPredDF$psp_mean,";",Inf),2,as.numeric))),
                 stdev = c(yPredDF$scs_sd, yPredDF$psp_sd),
                 col = c(rep(plottingCols["PCP"], nrow(yPredDF)),
                         rep(plottingCols["PSP"], nrow(yPredDF))))
  plot(p1~residue, data=p,
       type = "h", lwd = 2,
       col = col,
       xlab = "", ylab = "cleavage/splicing strength (%)", axes=T,
       main = Quant$substrateID[1])
  arrows(p$residue, p$p1-p$stdev, p$residue, p$p1+p$stdev,
         length=0.03, angle=90, code=3,
         lty = "solid", lwd = 1, col = p$col) %>%
    suppressWarnings()
  
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
  ALL = left_join(substrateCounts, yPredDF_filter %>% dplyr::select(substrateID, residue, col, nCol)) %>%
    rename(y = col,
           n = nCol) %>%
    na.omit()
  
  # separate by bio reps
  t = str_split_fixed(ALL$y, ";", Inf) %>% as.matrix()
  t = apply(t,2,as.numeric)
  
  X = ALL[,allCombos] %>% as.matrix()
  X = X[,!colnames(X) %in% c("P1;X", "P1_;X")]
  
  # merge and return
  DATA = list(X = X,
              t = t,
              n = ALL$n,
              substrateIDs = ALL$substrateID)
  
  return(DATA)
}

# PCP: SR1
DATA_P1pcp = getData(Quant = Target,
                  target = "P1", 
                  features_from = "SR1",
                  SRpos = c("P6"=-5,"P5"=-4,"P4"=-3, "P3"=-2, "P2"=-1, "P1"=0,
                            "P-1"=1, "P-2"=2, "P-3"=3, "P-4"=4, "P-5"=5, "P-6"=6),
                  col = "scs_mean", nCol = "scs_n")

DATA_P1psp = getData(Quant = Target,
                     target = "P1", 
                     features_from = "SR1",
                     SRpos = c("P6"=-5,"P5"=-4,"P4"=-3, "P3"=-2, "P2"=-1, "P1"=0,
                               "P-1"=1, "P-2"=2, "P-3"=3, "P-4"=4, "P-5"=5, "P-6"=6),
                     col = "psp_mean", nCol = "psp_n")

DATA_P1_psp = getData(Quant = Target,
                     target = "P1_",
                     features_from = "SR2",
                     SRpos = c("P-6_"=-6,"P-5_"=-5,"P-4_"=-4, "P-3_"=-3, "P-2_"=-2, "P-1_"=-1,
                               "P1_"=0, "P2_"=1, "P3_"=2, "P4_"=3, "P5_"=4, "P6_"=5),
                     col = "psp_mean", nCol = "psp_n")

# ----- apply parameters to create prediction -----
# tmp

rbPal = colorRampPalette(c('blue','red'))

getPrediction = function(DATA, params, nm) {
  
  X = DATA$X
  t = DATA$t
  
  dataNames = colnames(X)
  paramNames = colnames(params)
  
  paramsF = params[,paramNames %in% dataNames]
  
  
  tsims = apply(paramsF,1,function(p){
    t = X%*%log(p)
    return(t)
  })
  
  tsims_mean = apply(tsims,1,mean,na.rm = T)
  tsims_sd = apply(tsims,1,sd,na.rm = T)
  
  ttrue = log(t/100)
  ttrue[!is.finite(ttrue)] = NA
  ttrue_mean = apply(ttrue,1,mean,na.rm = T)
  ttrue_sd = apply(ttrue,1,sd,na.rm = T)
  
  keep = which(!(is.na(ttrue_sd) | is.na(ttrue_mean)))
  
  
  boxplot(tsims[keep,] %>% t(), outline=FALSE,
          lty = "blank", whisklty="blank", medlty = "solid")
  points(ttrue_mean[keep],col="red", cex = .6)
  
  mini = min(ttrue_mean[keep], ttrue_mean[keep], na.rm = T)
  maxi = max(ttrue_mean[keep], ttrue_mean[keep], na.rm = T)
  
  pcc = cor(ttrue_mean[keep], tsims_mean[keep])
  plot(x = ttrue_mean[keep], y = tsims_mean[keep],
       pch = 16, xlab = "true", ylab = "predicted",
       col = rbPal(50)[as.numeric(cut(log(DATA$n+1),breaks = 50))],
       # ylim = c(mini,maxi), xlim = c(mini,maxi),
       main = paste0(nm, ": ", DATA$substrateIDs[1],", PCC = ", round(pcc,4)))
  abline(a = 0, b = 1, col = "black")
  
  arrows(x0 = ttrue_mean[keep]+ttrue_sd[keep],
         x1 = ttrue_mean[keep]-ttrue_sd[keep],
         y0 = tsims_mean[keep], y1 = tsims_mean[keep],
         code = 3, angle = 90, length = 0.03, lwd = .5,
         col = rbPal(50)[as.numeric(cut(log(DATA$n+1),breaks = 50))])
  
  arrows(x0 = ttrue_mean[keep], x1 = ttrue_mean[keep],
         y0 = tsims_mean[keep]+tsims_sd[keep], y1 = tsims_mean[keep]-tsims_sd[keep],
         code = 3, angle = 90, length = 0.03, lwd = .5,
         col = rbPal(50)[as.numeric(cut(log(DATA$n+1),breaks = 50))])
  
  # text(labels = c(1:nrow(t))[keep],x = ttrue_mean[keep], y = tsims_mean[keep]+0.01)
  
}

suppressWarnings(dir.create("results/Bayesian_ProteaSMM/PREDICTION/"))
suppressWarnings(dir.create("results/Bayesian_ProteaSMM/PREDICTION/POLYPEPTIDES/"))
suppressWarnings(dir.create("results/Bayesian_ProteaSMM/PREDICTION/POLYPEPTIDES_ANALYTICAL/"))

png("results/Bayesian_ProteaSMM/PREDICTION/POLYPEPTIDES_ANALYTICAL/TSN5_220616.png", height = 9, width = 14, units = "in", res = 600)
par(mfrow = c(3,2))
getPrediction(DATA = DATA_P1pcp, params = pcp_chain, nm = "PCP")
getPrediction(DATA = DATA_P1psp, params = sr1_chain, nm = "SR1")
getPrediction(DATA = DATA_P1_psp, params = sr2_chain, nm = "SR2")
dev.off()

### OUTPUT ###


