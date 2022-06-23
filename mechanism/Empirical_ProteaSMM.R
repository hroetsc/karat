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
library(Rlinsolve)
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

WeightsOptimisation = function(inpFolder, nm, target, interesting_residues) {
  
  # ----- load data + preprocessing -----
  load(paste0(inpFolder,"DATA.RData"))
  X = DATA$X
  
  tpure = rowMeans(DATA$t, na.rm = T) %>% as.matrix()
  t = log(tpure/100+1e-08)
  t[!is.finite(t)] = NA
  
  
  # tmp !
  k = which(DATA$substrateIDs != "MM136")
  tpure = tpure[k, ] %>% as.matrix()
  t = t[k,] %>% as.matrix()
  X = X[k,]
  
  folderN = str_replace_all(inpFolder, "data/ProteaSMM/", "results/Empirical_ProteaSMM/")
  suppressWarnings(dir.create(folderN, recursive = T))
  
  
  # ----- define optimisation function -----
  numParam = ncol(X)+2
  w0 = runif(n = numParam, min = 0, max = 2)
  
  E = function(w) {
    t0 = w[length(w)-1]
    lambda = w[length(w)]
    p = w[1:(length(w)-2)]
    p[p<=0] = 1e-08
    
    loss = (X%*%log(p)+t0 - t)^2 + lambda*(X%*%log(p))^2
    return(sum(loss))
  }
  
  opt = optim(par = w0, E, method = "SANN", control = list(maxit = 1e03), hessian = T)
  w = opt$par
  mse = opt$value
  
  tsim = X%*%log(w[1:(length(w)-2)])
  cor(t, tsim)
  
  kk = which(tpure == 0)
  tpure[kk, ] = rnorm(n = length(kk), mean = 0, sd = 1e-03)
  
  opt2 = lsolve.bicgstab(A=X[-kk,], B=tpure[-kk,], reltol = 1e-12)
  tsim = X[-kk,]%*%opt2$x
  cor(t[-kk,],tsim)
  
  density(t) %>% plot()
  density(tsim) %>% lines()
  
  
}

WeightsOptimisation(inpFolder = "data/ProteaSMM/PSP_STS/", nm = "PSP_STS", target = "P1+1_",
                 interesting_residues = c("P4", "P3", "P2", "P1","P-1", "P-2", "P-3", "P-4", "P-4_", "P-3_", "P-2_", "P-1_", "P1_", "P2_", "P3_", "P4_"))

