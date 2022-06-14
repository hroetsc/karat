### karat projetc - PCPS mechanism ###
# description:  analyse posterior distributions of Bayesian ProteaSMM
# input:        BayesianTools: out.RData, posterior.RData
# output:       importance of each aa for SCS-PSP-P1
# author:       HPR

library(dplyr)
library(stringr)
library(ggplot2)
library(BayesianTools)
library(numDeriv)
library(corrplot)
library(RColorBrewer)
source("src/invitroSPI_utils.R")
source("src/plotting_utils.R")

theme_set(theme_classic())
AAchar_here = c("P","G","C","M","A","V","L","F","Y","W","H","R","K","D","E","N","Q","S","T","X")

### INPUT ###
load("Bayesian_ProteaSMM/server/PCP_SR1feat_P1_logModel_server6_polyRepshighNoise/posterior.RData")
pcp_post = posterior
load("Bayesian_ProteaSMM/server/PSP_SR1feat_P1_logModel_server6_polyRepshighNoise/posterior.RData")
sr1_post = posterior
load("Bayesian_ProteaSMM/server/PSP_SR2feat_P1__logModel_server6_polyRepshighNoise/posterior.RData")
sr2_post = posterior
load("Bayesian_ProteaSMM/server/PSP_SR1+1feat_logModel_server6_polyRepshighNoise/posterior.RData")
sr_post = posterior

### MAIN PART ###
suppressWarnings(dir.create("results/Bayesian_ProteaSMM/PLOTS/"))
suppressWarnings(dir.create("results/Bayesian_ProteaSMM/PLOTS/PCP/"))
suppressWarnings(dir.create("results/Bayesian_ProteaSMM/PLOTS/ALL/"))


# ----- correlation between true and predicted -----
plotCorrelation = function(inpFolder, posterior) {
  
  # ---- load data
  if(length(inpFolder) == 1) {
    
    load(paste0(inpFolder,"DATA.RData"))
    X = DATA$X
    
    t = DATA$t
    # t = rowMeans(DATA$t, na.rm = T) %>% as.matrix()
    t = log(t/100)
    t[!is.finite(t)] = NA
    
    k = grep("MM", DATA$substrateIDs)
    X = X[k,]
    t = t[k,] %>% as.matrix()
    
  } else {
    
    load(paste0(inpFolder[1],"DATA.RData"))
    X1 = DATA$X
    t1 = DATA$t
    # t1 = rowMeans(DATA$t, na.rm = T) %>% as.matrix()
    # t1 = log(t1/100)
    t1 = log(t1)
    t1[!is.finite(t1)] = NA
    # only polypeptides
    k = grep("MM", DATA$substrateIDs)
    X1 = X1[k,]
    t1 = t1[k,] %>% as.matrix()
    
    
    
    load(paste0(inpFolder[2],"DATA.RData"))
    X2 = DATA$X
    # t2 = rowMeans(DATA$t, na.rm = T) %>% as.matrix()
    t2 = DATA$t
    # t2 = log(t2/100)
    t2 = log(t2)
    t2[!is.finite(t2)] = NA
    # only polypeptides
    k = grep("MM", DATA$substrateIDs)
    X2 = X2[k,]
    t2 = t2[k,] %>% as.matrix()
    
    
    # concatenate matrices and targets
    X = rbind(X1,X2)
    t = rbind(t1,t2)
    
    i1 = c(1:ncol(X1))
    i2 = c((ncol(X1)+1):(ncol(X1)+ncol(X2)))
    
  }
  
  N = dim(posterior)[1]
  burnIn = round(0.9*N)
  chain = posterior[-(1:burnIn),]
  if(length(inpFolder) == 1) {
    colnames(chain) = c(colnames(X),"sigma")
  } else {
    colnames(chain) = c(colnames(X1),colnames(X2),"sigma")
  }
  
  numRep = ncol(t)
  
  # ---- get chain and tsim
  NN = min(c(dim(chain)[1],100))
  
  sampledIndex = sample(c(1:dim(chain)[1]),size=NN,replace=FALSE)
  param = chain[sampledIndex,]
  p = param[,-dim(param)[2]]
  
  if (length(inpFolder) == 1) {
    tsim = matrix(NA,dim(param)[1],nrow(t))
    for(i in 1:dim(param)[1]){
      tsim[i,] = X%*%log(p[i,])
    }
    
  } else {
    
    tsim1 = matrix(NA,dim(param)[1],nrow(X1))
    tsim2 = matrix(NA,dim(param)[1],nrow(X2))
    for(i in 1:dim(param)[1]){
      tsim1[i,] = X1%*%log(p[i,i1])
      tsim2[i,] = X2%*%log(p[i,i2])
    }
    
    tsim = cbind(tsim1, tsim2)
  }
  
  tsim_mean = apply(tsim,2,mean)
  
  # ----- plot correlation
  # par(mfrow = c(1,numRep))
  sapply(1:numRep, function(r){
    rr = which(!is.na(t[,r]))
    pcc = cor(x = t[rr,r], y = tsim_mean[rr], method = "pearson")
    
    plot(x=t[rr,r], y=tsim_mean[rr], pch = 16,
         # xlim = c(-13,0), ylim = c(-13,0),
         xlab = "log SCS-P1", ylab = "simulated log SCS-P1",
         main = paste0("PCC = ", round(pcc,4)), sub = inpFolder[1])
    abline(a=0,b=1,col = "red")
  })
  
  # return chain
  return(chain)
}

pdf("results/Bayesian_ProteaSMM/PLOTS/ALL/_SR1+2tog_CORRELATION.pdf",width = 10, height = 10)
par(mfrow = c(2,3))
pcp_chain = plotCorrelation(inpFolder = "Bayesian_ProteaSMM/ProteaSMM/PCP_SR1feat_P1/", posterior = pcp_post)
sr_chain = plotCorrelation(inpFolder = c("Bayesian_ProteaSMM/ProteaSMM/PSP_SR1feat_P1/", "Bayesian_ProteaSMM/ProteaSMM/PSP_SR2feat_P1_/"), posterior = sr_post)
# sr1_chain = plotCorrelation(inpFolder = "data/ProteaSMM/PSP_SR1feat_P1/", posterior = sr1_post)
# sr2_chain = plotCorrelation(inpFolder = "data/ProteaSMM/PSP_SR2feat_P1_/", posterior = sr2_post)
dev.off()

# ----- boxplots -----
PCPparam = tidyr::gather(pcp_chain %>% as.data.frame()) %>%
  mutate(reactant = "PCP")
# SR1param = tidyr::gather(sr1_chain %>% as.data.frame()) %>%
#   mutate(reactant = "SR1")
# SR2param = tidyr::gather(sr2_chain %>% as.data.frame()) %>%
#   mutate(reactant = "SR2")
SRparam = tidyr::gather(sr_chain %>% as.data.frame()) %>%
  mutate(reactant = "SR")

# bothR = rbind(SR1param,SR2param,PCPparam)
bothR = rbind(SRparam,PCPparam)
bothR$position = str_split_fixed(bothR$key,";",Inf)[,1]
# tmp!
bothR$reactant[bothR$position %in% c("P4", "P3", "P2", "P1", "P-1", "P-2", "P-3", "P-4") & bothR$reactant != "PCP"] = "SR1"
bothR$reactant[bothR$position %in% c("P-4_", "P-3_", "P-2_", "P-1_", "P1_", "P2_", "P3_", "P4_")] = "SR2"

bothR$aa = str_split_fixed(bothR$key,";",Inf)[,2]

bothR = bothR[-which(bothR$position == "sigma"), ]
# reorder amino acids
bothR$aa = factor(bothR$aa, levels = AAchar_here)

# compare positions
# bothR$commonPos = NA
# bothR$commonPos[bothR$position %in% c("P1","P1_")] = "P1/P1_"
# bothR$commonPos[bothR$position %in% c("P2","P2_")] = "P2/P2_"
# bothR$commonPos[bothR$position %in% c("P3","P3_")] = "P3/P3_"
# bothR$commonPos[bothR$position %in% c("P4","P4_")] = "P4/P4_"
# bothR$commonPos[bothR$position %in% c("P-1","P-1_")] = "P-1/P-1_"
# bothR$commonPos[bothR$position %in% c("P-2","P-2_")] = "P-2/P-2_"
# bothR$commonPos[bothR$position %in% c("P-3","P-3_")] = "P-3/P-3_"
# bothR$commonPos[bothR$position %in% c("P-4","P-4_")] = "P-4/P-4_"


# plot per position
# pos = c("P4", "P3", "P2", "P1","P-4","P-3","P-2","P-1")
pos = c("P4/P4_","P3/P3_","P2/P2_","P1/P1_","P-4/P-4_","P-3/P-3_","P-2/P-2_","P-1/P-1_")
allP = list()

for (p in 1:length(pos)) {
  cntR = bothR[bothR$commonPos == as.character(pos[p]), ]
  
  cntP = ggplot(cntR, aes(x = aa, y = value, fill = reactant)) +
    geom_boxplot(outlier.shape = NA, coef = 0, position = position_dodge(0.5), aes(alpha = .8)) +
    geom_hline(yintercept = c(mean(cntR$value[cntR$reactant == "PCP"]),
                              mean(cntR$value[cntR$reactant == "SR1"]),
                              mean(cntR$value[cntR$reactant == "SR2"])),
               col = c(plottingCols[["PCP"]],"gray","lightblue"), lty = "dashed", lwd = 0.8) +
    # ylim(c(-6,6)) +
    scale_fill_manual(values = c(plottingCols[["PCP"]],"gray","lightblue")) +
    xlab("amino acid") +
    ylab("regression weight") +
    ggtitle(pos[p])
  
  allP[[p]] = cntP
  
}

ggsave(filename = "results/Bayesian_ProteaSMM/PLOTS/ALL/_SR1+2tog_paramDistributions.pdf", 
       plot = gridExtra::marrangeGrob(allP, nrow=4, ncol=2, byrow = F), 
       width = 15, height = 18, dpi = "retina")



# plot per amino acid
bothR$commonPos = factor(bothR$commonPos, levels = c("P4/P4_","P3/P3_","P2/P2_","P1/P1_","P-1/P-1_","P-2/P-2_","P-3/P-3_","P-4/P-4_"))
allPaa = list()
for (a in 1:length(AAchar_here)) {
  
  cntR = bothR[bothR$aa == as.character(AAchar_here[a]), ]
  cntP = ggplot(cntR, aes(x = commonPos, y = value, fill = reactant)) +
    geom_boxplot(outlier.shape = NA, coef = 0, position = position_dodge(0.5), aes(alpha = .8)) +
    geom_hline(yintercept = c(mean(cntR$value[cntR$reactant == "PCP"]),
                              mean(cntR$value[cntR$reactant == "SR1"]),
                              mean(cntR$value[cntR$reactant == "SR2"])),
               col = c(plottingCols[["PCP"]],"gray","lightblue"), lty = "dashed", lwd = 0.8) +
    # ylim(c(-6,6)) +
    scale_fill_manual(values = c(plottingCols[["PCP"]],"gray","lightblue")) +
    xlab("position") +
    ylab("regression weight") +
    ggtitle(AAchar_here[a])
  
  allPaa[[a]] = cntP
}

ggsave(filename = "results/Bayesian_ProteaSMM/PLOTS/ALL/_SR1+2tog_paramDistributions_aa-wise.pdf", 
       plot = gridExtra::marrangeGrob(allPaa, nrow=7, ncol=3, byrow = T), 
       width = 22, height = 35, dpi = "retina")


# ----- parameter correlation -----
myCols = c(plottingCols[["PCP"]],"gray","lightblue")
bothR = bothR %>%
  mutate(col = ifelse(reactant == "SR1", "gray", "lightblue"),
         col = ifelse(reactant == "PCP", plottingCols["PCP"], col))

pdf("results/Bayesian_ProteaSMM/PLOTS/ALL/_SR1+2tog_PARAMETER_CORRELATION.pdf", height = 7, width = 21)
par(mfrow = c(1,3))

poi = c("P4/P4_","P3/P3_","P2/P2_","P1/P1_","P-1/P-1_","P-2/P-2_","P-3/P-3_","P-4/P-4_")
for (p in 1:length(poi)) {
  print(poi[p])
  cntR = bothR[bothR$commonPos == as.character(poi[p]), ] %>%
    mutate(aa = as.character(aa)) %>%
    as.data.frame()
  
  sapply(c("SR1","SR2","PCP"), function(t) {
    cnt = sapply(AAchar_here[-which(AAchar_here == "X")], function(a){
      cntR$value[cntR$aa == a & cntR$reactant == t]
    })
    
    M = cor(cnt, method = "pearson")
    corrplot(M, type="upper", order = "hclust",
             col=brewer.pal(n=8, name="RdYlBu"), tl.col = "black",
             diag = T, addCoefasPercent = F, title = paste0("\n",t," - ",poi[p]))
  })
  
}
dev.off()

pdf("results/Bayesian_ProteaSMM/PLOTS/ALL/_SR1+2tog_PARAMETER_CORRELATION_allPos.pdf", height = 25, width = 25)
sapply(c("SR1","SR2","PCP"), function(t) {
  print(t)
  
  cntR = bothR[bothR$reactant == t, ]
  allVal = sapply(unique(cntR$key), function(k){
    cntR$value[cntR$key == k & cntR$reactant == t]
  })
  M = cor(allVal)
  corrplot(M, type="upper", order = "hclust",
           col=brewer.pal(n=8, name="RdYlBu"),
           tl.cex	= .4, tl.col = "black", 
           diag = T, addCoefasPercent = F, title = paste0("\n",t))
})
dev.off()


# ----- analyse correlation matrix ------
cntR = bothR[bothR$reactant == "PCP", ]
allVal = sapply(unique(cntR$key), function(k){
  cntR$value[cntR$key == k & cntR$reactant == "PCP"]
})
M = cor(allVal)

ndx = sort(M[upper.tri(M)], decreasing = T)
res = sapply(ndx, function(i){
  return(c(i,rownames(which(M == i, arr.ind = T, useNames = T))))
}) %>% 
  t() %>%
  as.data.frame()
names(res) = c("correlation", "param1", "param2")
res$correlation = as.numeric(res$correlation)

res$param1_val = sapply(res$param1, function(p){
  median(bothR$value[bothR$reactant == "PCP" & bothR$key == p])
})

res$param2_val = sapply(res$param2, function(p){
  median(bothR$value[bothR$reactant == "PCP" & bothR$key == p])
})



# ----- curvature matrix -----

w = colMeans(chain)
w = w[-length(w)]

getPred = function(x0) {
  X <- matrix(x0[1:159*1297],ncol = 159,nrow = 1297); w <- x0[(159*1297+1):length(x0)]
  t <- X%*%log(w)
  t
}

numDeriv::hessian(func=getPred, x= c(X,w))
pracma::hessian(f = getPred, x0 = c(X,w))

numDeriv::jacobian(func = getPred, x = w)

