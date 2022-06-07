### karat projetc - PCPS mechanism ###
# description:  ProteaSMM implementation to determine relecance of each AA for SCS/PSP-P1
# input:        aSPIre: Roetschke SciData, EGFR ery, WT substrates (quantitative DB)
#               random database for quantitative DB (sanity check)
# output:       importance of each aa for SCS-PSP-P1
# author:       HPR


library(caret)
library(dplyr)
library(stringr)
library(reshape2)
library(matlib)
library(RColorBrewer)
library(dbscan)
library(MASS)
library(pracma)
source("src/invitroSPI_utils.R")
source("src/_extract-aa.R")
source("src/SCS+PSP-P1.R")
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(100)


target = "P1"
features_from = "SR1+SR2"
# features_from = "SR2"

# SRpos = c("P4"=-3, "P3"=-2, "P2"=-1, "P1"=0,
#           "P-1"=1, "P-2"=2)
# SRpos = c("P-2_"=-2, "P-1_"=-1,"P1_"=0, "P2_"=1, "P3_"=2, "P4_"=3)

SRpos = c("P4"=-3, "P3"=-2, "P2"=-1, "P1"=0,
          "P-1"=1, "P-2"=2,
          "P1_"=0, "P2_"=1, "P3_"=2, "P4_"=3)

interesting_residues = names(SRpos)
AAchar_here = c("P","G","C","M","A","V","L","F","Y","W","H","R","K","D","E","N","Q","S","T","X")
AAchar_here_sorted = sort(AAchar_here)

allCombos = tidyr::crossing(interesting_residues,AAchar_here_sorted)
allCombos = do.call(paste, c(allCombos[c("interesting_residues","AAchar_here_sorted")], sep=";"))


### INPUT ###
load("data/aSPIre.RData")
# load("data/randomDB_Quant_aSPIre.RData")

### MAIN PART ###
suppressWarnings(dir.create("results/ProteaSMM/"))

# ----- preprocessing -----
Quant = Kinetics %>%
  ILredundancy() %>%
  disentangleMultimappers.Type() %>%
  disentangleMultimappers.SRlen() %>%
  disentangleMultimappers.IVSeqLen() %>%
  removeMultimappers.SRlen() %>%
  removeMultimappers.AA() %>%
  tidyr::separate_rows(digestTimes, intensities, sep=";") %>%
  rename(digestTime = digestTimes,
         intensity = intensities) %>%
  mutate(intensity = as.numeric(intensity),
         digestTime = as.numeric(digestTime)) %>%
  filter(digestTime == 4)


extractSRs = function(DB) {
  
  pos = str_split_fixed(DB$positions, "_", Inf)[,c(1:4)]
  pos = apply(pos,2,as.numeric)
  
  DB$sr1 = substr(DB$substrateSeq, pos[,1], pos[,2])
  DB$sr2 = substr(DB$substrateSeq, pos[,3], pos[,4])
  
  DB$posP1 = pos[,2]
  DB$posP1_ = pos[,3]
  return(DB)
}

Quant = extractSRs(Quant)

# ----- calculate SCS/PSP-P1 -----
out = SCSandPSP_allSubs(Quant, target)
yPredDF = plyr::ldply(out) %>%
  as.data.frame()
yPredDF$target = target
names(yPredDF)[1] = "substrateID"

# ----- extract amino acid counts @ each position of the substrate -----
# get all peptide positions
DBQuant = left_join(Quant, extract_aminoacids(Quant, onlyValidSeq = T)) %>%
  filter(biological_replicate == 1) %>%
  unique()

# get P-1 and P-2
pos = str_split_fixed(DBQuant$positions,"_",Inf)[,c(1:4)]
if (features_from %in% c("SR1","SR1+SR2")) {
  DBQuant$`P-1` = substr(DBQuant$substrateSeq, as.numeric(pos[,2])+1, as.numeric(pos[,2])+1)
  DBQuant$`P-2` = substr(DBQuant$substrateSeq, as.numeric(pos[,2])+2, as.numeric(pos[,2])+2)
} else if (features_from == "SR2") {
  DBQuant$`P-1_` = substr(DBQuant$substrateSeq, as.numeric(pos[,3])-1, as.numeric(pos[,3])-1)
  DBQuant$`P-2_` = substr(DBQuant$substrateSeq, as.numeric(pos[,3])-2, as.numeric(pos[,3])-2)
}

DBQuant[,SRnames][DBQuant[,SRnames] == ""] = "X"

# get all positions for all SR1s in each substrate
# input: either PSPs or PCPs
# DB = DBQuant[DBQuant$productType == "PSP", ]  # tmp!
getSubstrateCounts = function(DB, substrateOnly = T) {
  
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
    
    
    if (!substrateOnly) {
      
      # get counts of peptides for each position
      cnt = DB[k,]
      
      allRES = lapply(pos$residue, function(r){
        if (target == "P1"){
          kk = which(cnt$posP1 == r)
        } else if (target == "P1_") {
          kk = which(cnt$posP1_ == r)
        }
        
        if (length(kk) > 0) {
          dd = cnt[kk,interesting_residues]
          
          res = lapply(1:ncol(dd),function(p){
            length(dd[p][dd[p]==SRTBL[r,p]])
          })
          res = as.numeric(unlist(res))
          names(res) = interesting_residues
          return(res)
        }
        
      }) %>%
        plyr::ldply() %>%
        t() %>%
        melt() %>%
        rename(residue = Var2, position = Var1, count = value) %>%
        mutate(residue = as.numeric(residue))
      
      master = melt(as.matrix(SRTBL)) %>%
        rename(position = Var2, residue = Var1, aa = value) %>%
        mutate(position = as.character(position)) %>%
        left_join(allRES)
      
      master$count[is.na(master$count)] = 0
      # normalise aa-wise
      # master = master %>%
      #   group_by(aa) %>%
      #   mutate(count = count/sum(count)) %>%
      #   ungroup()
      
      master$pos = do.call(paste, c(master[c("position","aa")], sep=";"))
      
      master2 = master %>%
        select(-aa, -position) %>%
        tidyr::spread(pos, count)
      master2[is.na(master2)] = 0
      missing = allCombos[! allCombos %in% names(master2)]
      
      if (length(missing) > 0) {
        app = matrix(0, nrow = nrow(master2), ncol = length(missing))
        colnames(app) = missing
        master2 = cbind(master2, app) %>% as.data.frame()
      }
      master2 = master2[,c("residue", allCombos)]
      
      return(master2)
      
    } else {
      
      master = matrix(0, nrow = nrow(SRTBL), ncol = length(allCombos))
      colnames(master) = allCombos
      
      for (j in 1:nrow(SRTBL)) {
        cntN = paste(colnames(SRTBL),SRTBL[j,],sep = ";")
        master[j,cntN] = 1
      }
      
      return(cbind(c(1:nchar(S)),master) %>% as.data.frame())
    }
    
  })
  names(allFeatures) = subSeqs
  
  return(allFeatures)
}

PSPsubstrateCounts = getSubstrateCounts(DBQuant[DBQuant$productType == "PSP", ])
PCPsubstrateCounts = getSubstrateCounts(DBQuant[DBQuant$productType == "PCP", ])

# ----- match arguments and prediction targets -----

PSPfeatures = plyr::ldply(PSPsubstrateCounts)
names(PSPfeatures)[1:2] = c("substrateID","residue")
PCPfeatures = plyr::ldply(PCPsubstrateCounts)
names(PCPfeatures)[1:2] = c("substrateID","residue")

# check which combinations never occur and remove them to reduce the matrix sparsity
kpsp = which(colSums(PSPfeatures[,allCombos]) == 0)
kpcp = which(colSums(PCPfeatures[,allCombos]) == 0)
rem = intersect(kpsp, kpcp)
# P1; X --> do not remove!
# all other combinations occur


# add prediction target
PSPall = left_join(PSPfeatures, yPredDF %>% dplyr::select(substrateID, residue, psp_mean)) %>%
  rename(y = psp_mean)
PCPall = left_join(PCPfeatures, yPredDF %>% dplyr::select(substrateID, residue, scs_mean)) %>%
  rename(y = scs_mean)


# ----- estimate lambda via cross-validation ------

# lambdas = c(seq(0.1,1,0.1), seq(2,10,1), seq(20,100,10), seq(200,1000,100))
lambdas = c(seq(0.05,3,0.05))
fold = 10
noBags = 5e02
split = 0.9

determineWeights = function(ALL, nm) {
  
  png(paste0("results/ProteaSMM/SMManalytical_",nm,"_",target,"_",features_from,"feat.png"), height = 20, width = 13, units = "in", res = 300)
  par(mfrow = c(3,2))
  
  # get data
  X = ALL[,allCombos] %>% as.matrix()
  X[is.na(X)] = 0
  t = ALL$y
  
  # k-fold cross-validation to determine optimal lambda
  flds = createFolds(c(1:nrow(X)), k = fold)
  
  Lcv = sapply(c(1:fold), function(i){
    # split into train/test data set in current fold
    k = flds[[i]]
    
    Xtest = X[k,]
    t_test = t[k]
    
    Xtrain = X[-k,]
    t_train = t[-k]
    
    MSE_lambda = sapply(lambdas, function(a){
      
      # determine weights using the train set
      A = diag(a, nrow = ncol(Xtrain), ncol = nrow(Xtrain))
      w = (ginv(Xtrain)+A)%*%t_train
      
      # get the loss on the test data set
      t_pred = Xtest%*%w
      MSE = sum((t_pred - t_test)^2)
      return(MSE)
    })
    
    bestL = lambdas[which.min(MSE_lambda)]
    return(c(bestL = bestL, MSE = min(MSE_lambda)))
    
  })
  
  bestL = median(Lcv[1,])
  plot(x = c(1:fold), y = Lcv[2,], pch = 16,
       ylim = c(min(Lcv[2,]), max(Lcv[2,])+500),
       main = "MSE across folds", xlab = "fold", ylab = "best MSE",
       sub = paste0("λ = ", bestL))
  text(Lcv[1,], x = c(1:fold), y = Lcv[2,]+100)
  
  # determine weights using bagging strategy
  MSEs = rep(NA, noBags)
  allParams = matrix(NA, nrow = noBags, ncol = ncol(X))
  bestL = 1
  
  pb = txtProgressBar(min = 0, max = noBags, style = 3)
  for (b in 1:noBags) {
    setTxtProgressBar(pb, b)
    
    # sample subset of the data
    kk = sample(nrow(X), ceiling(nrow(X)*split))
    Xcnt = X[kk,]
    t_cnt = t[kk]
    
    # get weights
    A = diag(bestL, nrow = ncol(Xcnt), ncol = nrow(Xcnt))
    # w = ((inv(t(Xcnt)%*%Xcnt + A))%*%t(Xcnt))%*%t_cnt
    w = (ginv(Xcnt)+A)%*%t_cnt
    
    # make prediction
    t_pred = X[-kk,]%*%w
    MSEs[b] = sum((t_pred - t[-kk])^2)
    
    allParams[b,] = w
  }
  
  # keep best 5 % and average over the parameters
  cutoff = quantile(MSEs, 1)
  ii = which(MSEs <= cutoff)
  
  # plot the bags
  plot(x = c(1:noBags), y = MSEs, type = "l",
       main = "bootstrap aggregation",
       xlab = "bag", ylab = "MSE")
  abline(h = cutoff, lty = "dashed", col = "red")
  
  params = colMeans(allParams[ii,])
  names(params) = allCombos
  
  
  # visualise parameters
  DF = cbind(params, str_split_fixed(allCombos,";",Inf)) %>%
    as.data.frame()
  names(DF) = c("mean_weight","position","aa")
  DFt = DF %>% tidyr::spread(position,mean_weight)
  DFt = as.matrix(DFt[,interesting_residues])
  DFt = apply(DFt,2,as.numeric)
  rownames(DFt) = AAchar_here_sorted
  DFt = DFt[AAchar_here, ]
  
  image(DFt %>% t(), axes = F, col = r,
        main = "estimated regression parameters",
        sub = "low: blue, high: red")
  axis(2, at = seq(0,1,1/(length(AAchar_here)-1)), labels = AAchar_here)
  axis(1, at = seq(0,1,1/(length(interesting_residues)-1)), labels = interesting_residues)
  
  # get parameter distribution
  boxplot(DFt, main = "multiple linear regression weights", sub = "distribution over positions")
  
  # cl = hdbscan(DFt[,colnames(DFt) %in% c("P1","P2","P-1")], minPts = 2)
  cl = hdbscan(DFt, minPts = 2)
  plot(cl$hc,
       main = "HDBSCAN* Hierarchy",
       # xlab = "on linear regression weights @ P2 - P-1",
       labels = gsub("L","I/L",rownames(DFt)))
  
  dev.off()
  
  # output values
  save(DFt, file = paste0("results/ProteaSMM/SMManalytical_",nm,"_",target,"_",features_from,"feat.RData"))
}

determineWeights(PSPall, nm = "PSP")

# only meaningful for SCS-P1 and SR1
determineWeights(PCPall, nm = "PCP")

