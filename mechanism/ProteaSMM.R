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
library(ggplot2)
library(BayesianTools)
source("src/invitroSPI_utils.R")
source("src/_extract-aa.R")
source("src/SCS+PSP-P1.R")

rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
rcol <- rf(100)
theme_set(theme_classic())

AAchar_here = c("P","G","C","M","A","V","I","L","F","Y","W","H","R","K","D","E","N","Q","S","T","X")
AAchar_here_sorted = sort(AAchar_here)

# hyperparameters
suppressWarnings(dir.create("results/Analytical_ProteaSMM/"))

### INPUT ###


### MAIN PART ###
suppressWarnings(dir.create("results/ProteaSMM/"))

# ----- estimate lambda via cross-validation ------

noBags = 5e02
split = 0.8

determineWeights = function(inpFolder, nm, target, interesting_residues) {
  
  load(paste0(inpFolder,"DATA.RData"))
  X = DATA$X
  
  t = rowMeans(DATA$t, na.rm = T) %>% as.matrix()
  t = log(t/100+1e-08)
  t[!is.finite(t)] = NA
  
  
  # tmp !
  k = which(DATA$substrateIDs != "MM136")
  t = t[k,] %>% as.matrix()
  X = X[k,]
  
  folderN = str_replace_all(inpFolder, "data/ProteaSMM/", "results/Analytical_ProteaSMM/")
  suppressWarnings(dir.create(folderN))
  
  
  png(paste0(folderN,"SMManalytical_",nm,"_",target,".png"), height = 13, width = 13, units = "in", res = 300)
  par(mfrow = c(2,2))
  
  # get data
  paramNames = colnames(X)
  numRep = ncol(t)
  numParam = ncol(X)
  
  # determine weights using bagging strategy
  CORREL = rep(NA, noBags)
  allParams = matrix(NA, nrow = noBags, ncol = ncol(X))
  
  pb = txtProgressBar(min = 0, max = noBags, style = 3)
  for (b in 1:noBags) {
    setTxtProgressBar(pb, b)
    
    # sample subset of the data
    kk = sample(nrow(X), ceiling(nrow(X)*split))
    Xcnt = X[kk,]
    t_cnt = t[kk,] %>% as.matrix()
    
    # get weights
    # for each bio rep
    # lambda = diag(0.5, nrow = numParam, ncol = nrow(Xcnt))
    w = sapply(1:numRep,function(n){
      # nn = which(!is.na(t_cnt[,n]))
      w = exp((ginv(Xcnt))%*%t_cnt[,n])
    })
    allParams[b,] = rowMeans(w)
    
    # make prediction
    t_pred = sapply(1:numRep,function(n){
      X[-kk,]%*%log(w[,n])
    })
    # MSEs[b] = sum((rowMeans(t_pred) - apply(t[-kk,] %>% as.matrix(),1,mean,na.rm = T))^2)
    CORREL[b] = cor(rowMeans(t_pred), apply(t[-kk,] %>% as.matrix(),1,mean,na.rm = T))
  }
  colnames(allParams) = paramNames
  # keep best 5 % and average over the parameters
  # cutoff = quantile(MSEs, .8, na.rm = T)
  # ii = which(MSEs < cutoff)
  
  cutoff = quantile(CORREL, .2, na.rm = T)
  ii = which(CORREL > cutoff)
  
  xx = allParams[ii, ]
  save(xx, file = paste0(folderN,"SMManalytical_",nm,"_",target,".RData"))
  
  # plot the bags
  plot(x = c(1:noBags), y = CORREL, type = "l",
       main = "bootstrap aggregation",
       xlab = "bag", ylab = "PCC")
  abline(h = log(cutoff), lty = "dashed", col = "red")
  
  
  params = apply(allParams[ii,],2,function(x){
    return(c(mu = mean(x), confLower = quantile(x, .05), confUpper = quantile(x, .95)))
  }) %>% as.data.frame()
  names(params) = paramNames
  
  
  # visualise parameters
  DF = cbind(params[1,] %>% t(), str_split_fixed(paramNames,";",Inf)) %>%
    as.data.frame()
  names(DF) = c("mean_weight","position","aa")
  DFt = DF %>% tidyr::spread(position,mean_weight)
  
  DFt = as.matrix(DFt[,interesting_residues])
  DFt = apply(DFt,2,as.numeric)
  rownames(DFt) = AAchar_here_sorted
  DFt = DFt[AAchar_here, ]
  
  image(DFt %>% t(), axes = F, col = rcol,
        main = "estimated regression parameters",
        sub = "low: blue, high: red")
  axis(2, at = seq(0,1,1/(length(AAchar_here)-1)), labels = AAchar_here)
  axis(1, at = seq(0,1,1/(length(interesting_residues)-1)), labels = interesting_residues)
  
  # get parameter distribution
  boxplot(DFt, main = "multiple linear regression weights", sub = "distribution over positions")
  
  # cluster parameters
  cl = hdbscan(DFt[-which(rownames(DFt) == "X"),], minPts = 2)
  # plot(cl$hc,
  #      main = "HDBSCAN* Hierarchy",
  #      labels = gsub("L","I/L",rownames(DFt)[-which(rownames(DFt) == "X")]))
  
  # prediction
  tsim = X%*%log(apply(allParams[ii,],2,median))
  tmean = rowMeans(t, na.rm = T)
  pcc = cor(tmean, tsim)
  
  plot(x = tmean, y = tsim,
       pch = 16, xlab = "true", ylab = "predicted",
       main = paste0("PCC = ", round(pcc,4)))
  abline(a = 0, b = 1, col = "red")
  
  dev.off()
  
  
  pdf(paste0(folderN,"SMManalytical_",nm,"_",target,".pdf"), height = 21, width = 16)
  par(mfcol = c(4,2))
  
  sapply(unique(DF$position), function(r){
    print(r)
    k = which(colnames(allParams) %in% paste(r,AAchar_here,sep = ";"))
    cnt = allParams[ii,k[match(paste(r,AAchar_here,sep = ";"),colnames(allParams[,k]))]]
    colnames(cnt) = AAchar_here
    
    boxplot(cnt, main = r, outline = F)
    # if(!all(is.na(scores))) {
    #   points(log(scores[k[match(paste(r,AAchar_here,sep = ";"),colnames(allParams[,k]))],]),
    #          col = "red")
    # }
  })
  
  dev.off()
  
  # output values
  out = list(allParams = allParams[ii,],
             DFt = DFt)
  save(out, file = paste0(folderN,"SMManalytical_",nm,"_",target,".RData"))
}

determineWeights(inpFolder = "data/ProteaSMM/PSP_STS/", nm = "PSP_STS", target = "P1+1_",
                 interesting_residues = c("P4", "P3", "P2", "P1","P-1", "P-2", "P-3", "P-4", "P-4_", "P-3_", "P-2_", "P-1_", "P1_", "P2_", "P3_", "P4_"))


determineWeights(inpFolder = "data/ProteaSMM/PCP_SR1extfeat_P1/", nm = "PCP_LOV", target = "P1",
                 interesting_residues = c("P6","P5","P4", "P3", "P2", "P1", "P-1", "P-2", "P-3", "P-4","P-5","P-6"))

determineWeights(inpFolder = "data/ProteaSMM/PSP_SR1extfeat_P1/", nm = "PSP_LOV", target = "P1",
                 interesting_residues = c("P6","P5","P4", "P3", "P2", "P1", "P-1", "P-2", "P-3", "P-4","P-5","P-6"))

determineWeights(inpFolder = "data/ProteaSMM/PSP_SR2extfeat_P1_/", nm = "PSP_LOV", target = "P1_",
                 interesting_residues = c("P-6_","P-5_","P-4_", "P-3_", "P-2_", "P-1_", "P1_", "P2_", "P3_", "P4_", "P5_", "P6_"))


# ----- create prior -----
load("results/Analytical_ProteaSMM/PSP_SR1extfeat_P1/SMManalytical_PSP_P1.RData")
sr1_chain = out$allParams
load("results/Analytical_ProteaSMM/PSP_SR2extfeat_P1_/SMManalytical_PSP_P1_.RData")
sr2_chain = out$allParams
load("results/Analytical_ProteaSMM/PCP_SR1extfeat_P1/SMManalytical_PCP_P1.RData")
pcp_chain = out$allParams


sr1_chain = cbind(sr1_chain,runif(n = nrow(sr1_chain), min = 0, max = 10))
colnames(sr1_chain)[ncol(sr1_chain)] = "sigma"
prior = createTruncatedNormalPrior(mean = apply(sr1_chain,2,mean), sd = apply(sr1_chain,2,sd), upper = apply(sr1_chain,2,max), lower = apply(sr1_chain,2,min))
save(prior, file = "data/ProteaSMM/PSP_SR1extfeat_P1/analytical_prior.RData")

sr2_chain = cbind(sr2_chain,runif(n = nrow(sr2_chain), min = 0, max = 10))
colnames(sr2_chain)[ncol(sr2_chain)] = "sigma"
prior = createTruncatedNormalPrior(mean = apply(sr2_chain,2,mean), sd = apply(sr2_chain,2,sd), upper = apply(sr2_chain,2,max), lower = apply(sr2_chain,2,min))
save(prior, file = "data/ProteaSMM/PSP_SR2extfeat_P1_/analytical_prior.RData")

pcp_chain = cbind(pcp_chain,runif(n = nrow(pcp_chain), min = 0, max = 10))
colnames(pcp_chain)[ncol(pcp_chain)] = "sigma"
prior = createTruncatedNormalPrior(mean = apply(pcp_chain,2,mean), sd = apply(pcp_chain,2,sd), upper = apply(pcp_chain,2,max), lower = apply(pcp_chain,2,min))
save(prior, file = "data/ProteaSMM/PCP_SR1extfeat_P1/analytical_prior.RData")


# for STS
load("results/Analytical_ProteaSMM/PSP_STS/SMManalytical_PSP_STS_P1+1_.RData")
sr_chain = out$allParams
sr_chain = cbind(sr_chain,runif(n = nrow(sr_chain), min = 0, max = 10))
colnames(sr_chain)[ncol(sr_chain)] = "sigma"
prior = createTruncatedNormalPrior(mean = apply(sr_chain,2,mean), sd = 2*apply(sr_chain,2,sd), lower = rep(0, ncol(sr_chain)), upper = 2*apply(sr_chain,2,max))
save(prior, file = "data/ProteaSMM/PSP_STS/analytical_prior.RData")



# ----- compare P1 and P1' -----
suppressWarnings(dir.create("results/Analytical_ProteaSMM/PLOTS/"))

load("results/Analytical_ProteaSMM/PSP_SR1extfeat_P1/SMManalytical_PSP_P1.RData")
sr1 = out
load("results/Analytical_ProteaSMM/PSP_SR2extfeat_P1_/SMManalytical_PSP_P1_.RData")
sr2 = out
# load("results/Analytical_ProteaSMM/PSP_SR1+2feat_P1/SMManalytical_PSP_P1.RData")
# sr = out
load("results/Analytical_ProteaSMM/PCP_SR1extfeat_P1/SMManalytical_PCP_P1.RData")
pcp = out


SR1param = tidyr::gather(sr1$allParams %>% as.data.frame()) %>%
  mutate(reactant = "SR1")
SR2param = tidyr::gather(sr2$allParams %>% as.data.frame()) %>%
  mutate(reactant = "SR2")
# SRparam = tidyr::gather(sr$allParams %>% as.data.frame()) %>%
#   mutate(reactant = "SR")
PCPparam = tidyr::gather(pcp$allParams %>% as.data.frame()) %>%
  mutate(reactant = "PCP")


bothR = rbind(SR1param,SR2param,PCPparam)
# bothR = rbind(SRparam,PCPparam)
# bothR$value = exp(bothR$value)

bothR$position = str_split_fixed(bothR$key,";",Inf)[,1]

# tmp!
# bothR$reactant[bothR$position %in% c("P6","P5","P4", "P3", "P2", "P1", "P-1", "P-2", "P-3", "P-4") & bothR$reactant != "PCP"] = "SR1"
# bothR$reactant[bothR$position %in% c("P-4_", "P-3_", "P-2_", "P-1_", "P1_", "P2_", "P3_", "P4_")] = "SR2"

bothR$aa = str_split_fixed(bothR$key,";",Inf)[,2]

# compare positions
{
  bothR$commonPos = NA
  bothR$commonPos[bothR$position %in% c("P1","P1_")] = "P1/P1_"
  bothR$commonPos[bothR$position %in% c("P2","P2_")] = "P2/P2_"
  bothR$commonPos[bothR$position %in% c("P3","P3_")] = "P3/P3_"
  bothR$commonPos[bothR$position %in% c("P4","P4_")] = "P4/P4_"
  bothR$commonPos[bothR$position %in% c("P5","P5_")] = "P5/P5_"
  bothR$commonPos[bothR$position %in% c("P6","P6_")] = "P6/P6_"
  bothR$commonPos[bothR$position %in% c("P-1","P-1_")] = "P-1/P-1_"
  bothR$commonPos[bothR$position %in% c("P-2","P-2_")] = "P-2/P-2_"
  bothR$commonPos[bothR$position %in% c("P-3","P-3_")] = "P-3/P-3_"
  bothR$commonPos[bothR$position %in% c("P-4","P-4_")] = "P-4/P-4_"
  bothR$commonPos[bothR$position %in% c("P-5","P-5_")] = "P-5/P-5_"
  bothR$commonPos[bothR$position %in% c("P-6","P-6_")] = "P-6/P-6_"
}


# reorder amino acids
bothR$aa = factor(bothR$aa, levels = AAchar_here)
# reorder positions
# bothR$commonPos = factor(bothR$commonPos, levels = c("P4/P4_","P3/P3_","P2/P2_","P1/P1_","P-1/P-1_","P-2/P-2_","P-3/P-3_","P-4/P-4_"))


# plot
pos = c("P6/P6_","P5/P5_","P4/P4_","P3/P3_","P2/P2_","P1/P1_",
        "P-6/P-6_","P-5/P-5_","P-4/P-4_","P-3/P-3_","P-2/P-2_","P-1/P-1_")
allP = list()

for (p in 1:length(pos)) {
  cntR = bothR[bothR$commonPos == as.character(pos[p]), ]
  
  cntP = ggplot(cntR, aes(x = aa, y = value, fill = reactant)) +
    geom_boxplot(outlier.shape = NA, coef = 0, position = position_dodge(0.5), aes(alpha = .8)) +
    geom_hline(yintercept = c(mean(cntR$value[cntR$reactant == "PCP" & cntR$value <= 5]),
                              mean(cntR$value[cntR$reactant == "SR1" & cntR$value <= 5]),
                              mean(cntR$value[cntR$reactant == "SR2" & cntR$value <= 5])),
               col = c(plottingCols[["PCP"]],"gray","lightblue"), lty = "dashed", lwd = 0.8) +
    ylim(c(0,5)) +
    scale_fill_manual(values = c(plottingCols[["PCP"]],"gray","lightblue")) +
    xlab("amino acid") +
    ylab("regression weight") +
    ggtitle(pos[p])
  
  allP[[p]] = cntP
  
}


ggsave(filename = "results/Analytical_ProteaSMM/PLOTS/_paramDistributions.pdf", 
       plot = gridExtra::marrangeGrob(allP, nrow=6, ncol=2, byrow = F), 
       width = 15, height = 27, dpi = "retina")

# sum over positions

# plot per amino acid
bothR$position = factor(bothR$position, levels = c("P6","P5","P4", "P3", "P2", "P1",
                                                   "P-1","P-2","P-3","P-4","P-5","P-6",
                                                   "P-6_","P-5_","P-4_", "P-3_", "P-2_", "P-1_",
                                                   "P1_", "P2_", "P3_", "P4_", "P5_", "P6_"))
allPaa = list()
for (a in 1:length(AAchar_here)) {
  
  cntR = bothR[bothR$aa == as.character(AAchar_here[a]), ]
  cntP = ggplot(cntR, aes(x = position, y = value, fill = reactant)) +
    geom_boxplot(outlier.shape = NA, coef = 0, position = position_dodge(0.5), aes(alpha = .8)) +
    geom_hline(yintercept = c(mean(cntR$value[cntR$reactant == "PCP" & cntR$value <= 5]),
                              mean(cntR$value[cntR$reactant == "SR1" & cntR$value <= 5]),
                              mean(cntR$value[cntR$reactant == "SR2" & cntR$value <= 5])),
               col = c(plottingCols[["PCP"]],"gray","lightblue"), lty = "dashed", lwd = 0.8) +
    geom_vline(xintercept = c(6.5,12.5,18.5), lty = "dotted") +
    ylim(c(0,5)) +
    scale_fill_manual(values = c(plottingCols[["PCP"]],"gray","lightblue")) +
    xlab("position") +
    ylab("regression weight") +
    ggtitle(AAchar_here[a])
  
  allPaa[[a]] = cntP
}

ggsave(filename = "results/Analytical_ProteaSMM/PLOTS/_paramDistributions_aa-wise.pdf", 
       plot = gridExtra::marrangeGrob(allPaa, nrow=7, ncol=3, byrow = T), 
       width = 35, height = 35, dpi = "retina")


