### karat projetc - PCPS mechanism ###
# description:  predictive power of in silico datastes
# input:        BayesianTools: posterior.RData of in silico dataset
# output:       prediction performance of simulated dataset based on training data set size
# author:       HPR

library(dplyr)
library(stringr)
library(RColorBrewer)
library(ggplot2)
library(DescTools)
source("src/invitroSPI_utils.R")
source("src/plotting_utils.R")
source("src/_ROC-curve.R")


theme_set(theme_classic())
AAchar_here = c("P","G","C","M","A","V","I","L","F","Y","W","H","R","K","D","E","N","Q","S","T","X")
AAchar_here_sorted = sort(AAchar_here)

subIDs = c("MM499", "MM577", "MM578", "MM580", "MM835", "MM836", "MM537", "MM136", "MM539", "MM581", "MM336", "MM579", "TSN5","MM582")


### INPUT ###
fs = list.files("Bayesian_ProteaSMM/server/_insilico_singleSubs_PCP_0622/", pattern = "posterior.RData", recursive = T, full.names = T)
load("data/ProteaSMM/insilico/DATA.RData")


### MAIN PART ###
suppressWarnings(dir.create("results/Bayesian_ProteaSMM/PLOTS/insilico/"))

# ----- load posterior and data -----
pseudo = 1e-05

X = DATA$X
t = DATA$t

substrates = DATA$substrateIDs
paramNames = colnames(X)


# --- data
# make sure to get the same number of iterations in all chains
minLen = sapply(fs, function(f){
  load(f)
  return(dim(posterior)[1])
})
N = min(minLen)

CHAINS = lapply(fs, function(f){
  load(f)
  
  posterior = posterior[1:N, ]
  # N = dim(posterior)[1]
  burnIn = round(0.3*N)
  chain = posterior[-(1:burnIn),]
  colnames(chain) = c(colnames(X),"sigma")
  
  # sample 1 k particles
  chain = chain[sample(1:nrow(chain), 1e03, replace = F), ]
  
  return(chain)
})

names(CHAINS) = fs
CHAINSdf = plyr::ldply(CHAINS)

save(CHAINS, file = "results/Bayesian_ProteaSMM/PLOTS/insilico/chains.RData")

# ----- prediction on MM582 -----
numberOfSubstrates = str_extract_all(fs, pattern = "[:alnum:]+(?=/posterior.RData)", simplify = T) %>% as.vector() %>% as.numeric()

pdf("results/Bayesian_ProteaSMM/PLOTS/insilico/performance.pdf", height = 5*4, width = 3*4)
par(mfrow = c(5,3))

PCC = rep(NA, length(CHAINS))
allROCs = list()
for (i in order(numberOfSubstrates)) {
  
  keep = which(substrates == "MM582")
  
  # get prediction
  cntPosterior = CHAINS[[i]]
  cntPosterior = cntPosterior[,colnames(X)]
  tsims = apply(cntPosterior,1,function(p){
    tsim = X[keep,]%*%log(p)
    return(tsim)
  })
  
  tsims_mean = apply(tsims,1,mean,na.rm = T)
  tsims_sd = apply(tsims,1,sd,na.rm = T)
  
  ttrue_mean = t[keep, ]
  ttrue_sd = rep(0, length(t[keep, ]))
  
  
  # correlation coefficient
  pcc = cor(ttrue_mean, tsims_mean)
  
  plot(x = ttrue_mean, y = tsims_mean,
       pch = 16, xlab = "true", ylab = "predicted",
       main = paste0(numberOfSubstrates[i], " substrates, PCC = ", round(pcc,4)))
  
  abline(a = 0, b = 1, col = "red")
  abline(v = log(0.01+pseudo), lty = "dashed", col = "blue")
  abline(h = log(0.01+pseudo), lty = "dashed", col = "blue")
  
  arrows(x0 = ttrue_mean+ttrue_sd,
         x1 = ttrue_mean-ttrue_sd,
         y0 = tsims_mean, y1 = tsims_mean,
         code = 3, angle = 90, length = 0.03, lwd = .5) %>% suppressWarnings()
  
  arrows(x0 = ttrue_mean, x1 = ttrue_mean,
         y0 = tsims_mean+tsims_sd, y1 = tsims_mean-tsims_sd,
         code = 3, angle = 90, length = 0.03, lwd = .5) %>% suppressWarnings()
  
  PCC[i] = pcc
  allROCs[[i]] = getROCcurve(ttrue_mean, tsims_mean, substrate = numberOfSubstrates[i], threshPerc = 0.01, retValsOnly = T)
}

barplot(PCC)

dev.off()


names(PCC) = numberOfSubstrates
names(allROCs) = numberOfSubstrates

# ----- summarised performance of in silico dataset -----

png("results/Bayesian_ProteaSMM/PLOTS/insilico/pcc.png", height = 4, width = 6, units = "in", res = 300)
barplot(PCC[order(as.numeric(names(PCC)))],
        xlab = "number of substrates in training dataset",
        ylab = "Pearson Correlation Coefficient",
        ylim = c(0,1))
dev.off()

allROCsDF = allROCs %>% plyr::ldply()
allROCsDF$numSubs = factor(allROCsDF$.id, levels = sort(numberOfSubstrates))
allROCsDF = allROCsDF[order(as.numeric(allROCsDF$.id)), ]

roc = ggplot(allROCsDF, aes(col = numSubs, alpha = .8)) +
  geom_path(aes(1 - specificity, sensitivity), lwd = 1) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
  scale_color_manual("number of substrates",
                     values = rainbow(max(numberOfSubstrates)),
                     labels = paste0(unique(allROCsDF$numSubs), " - AUC: ", round(unique(allROCsDF$roc_auc), 2))) +
  xlim(c(0,1)) +
  ylim(c(0,1)) +
  ggtitle("ROC curve")


trueBinary = ifelse(t > log(0.01+pseudo), 1, 0)
pr = ggplot(allROCsDF, aes(col = numSubs, alpha = .8)) +
  geom_path(aes(recall, precision), lwd = 1) +
  geom_abline(intercept = length(which(trueBinary == 1))/length(trueBinary),
              slope = 0, linetype = "dotted") +
  xlim(c(0,1)) +
  ylim(c(0,1)) +
  scale_color_manual("number of substrates",
                     values = rainbow(max(numberOfSubstrates)),
                     labels = paste0(unique(allROCsDF$numSubs), " - AUC: ", round(unique(allROCsDF$pr_auc), 2))) +
  ggtitle("PR curve")

ggsave(filename = "results/Bayesian_ProteaSMM/PLOTS/insilico/roc_auc.png", plot = gridExtra::grid.arrange(roc, pr, ncol = 2),
       height = 4, width = 10, dpi = "retina")


### OUTPUT ###

