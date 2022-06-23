### karat projetc - PCPS mechanism ###
# description:  Bayesian inference for determining importance of aa for positions
# input:        substrate matrices for each residue in each substrate
# output:       importance of each aa for SCS-PSP-P1
# author:       JL, modified by HPR

library(BayesianTools)
library(RColorBrewer)
library(grDevices)
library(dplyr)
library(stringr)
library(MASS)
library(pracma)
library(optparse)
suppressWarnings(dir.create("results/Bayesian_ProteaSMM/"))


# ----- parse command line arguments -----
option_list = list(
  make_option(c('-f', '--folderN'), type = 'character', default = "",
              help = 'output folder'),
  make_option(c('-i', '--inpFolder'), type = 'character', default = "",
              help = 'input folder'),
  make_option(c('-n', '--substrate'), type = 'character', default = "",
              help = 'substrate used for inference'),
  make_option(c('-o', '--leaveOut'), type = 'character', default = "MM136",
              help = 'substrate left completely out'),
  make_option(c('-c', '--cntleaveOut'), type = 'character', default = "",
              help = 'substrate left out in current CV'))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

folderN = opt$folderN
inpFolder = opt$inpFolder
leaveOut = opt$leaveOut
cntleaveOut = opt$cntleaveOut
cntSubstrate = opt$substrate


pseudo = 1e-05

### INPUT ###
inpFolder = "data/ProteaSMM/PSP_STS/"
folderN = "results/Bayesian_ProteaSMM/STS_0621/"

load(paste0(inpFolder,"DATA.RData"))
# load(paste0(inpFolder, "analytical_prior.RData"))
X = DATA$X
t = DATA$t

# keep = which(rowSums(t) > 0)
# t = t[keep, ]
# X = X[keep, ]
# subIDs = DATA$substrateIDs[keep]

t = log(t/100+pseudo)
# t = log(t/100)
t[!is.finite(t)] = NA


k = which(!DATA$substrateIDs == "MM582")
X = X[k, ]
t = t[k, ]

suppressWarnings(dir.create(folderN, recursive = T))

numRep = ncol(t)
numParam = ncol(X)+1

### MAIN PART ###
# ---- settings -----
Niter = 5*10**5

mini = rep(0,numParam)  # 1 parameter more than required by model --> sigma
maxi = rep(1.3,numParam)

mini[length(mini)] = 0
maxi[length(maxi)] = 10 # initial sigma

# ----- plotting analytics -----

plotChain <- function(chain){
  
  # define parameters
  paramNames <- c(colnames(X),"sigma")
  
  N = dim(chain)[1]
  burnIn = round(0.9*N)
  
  
  # plot histograms and chains
  pdf(paste(folderN,"/chain_",N,".pdf",sep=""), width=9, height=12)
  par(mfrow = c(5,2))
  
  for(i in 1:dim(chain)[2]){
    hist(chain[-(1:burnIn),i], main=paramNames[i], xlab="",breaks=70,col="grey" )
    abline(v = mean(chain[-(1:burnIn),i]), col="red")
    plot(chain[,i],xlab="iteration",ylab=paramNames[i],main="",axes=FALSE,type="l")
    axis(1)
    axis(2)
  }
  
  
  dev.off()
  
  
  # get boxplots
  
  # compute and plot time series based on posterior sample
  post = chain[-(1:burnIn),]
  NN = min(c(dim(post)[1],100))
  
  sampledIndex = sample(c(1:dim(post)[1]),size=NN,replace=FALSE)
  param = post[sampledIndex,]
  
  p = param[,-dim(param)[2]]
  
  # compute likelihood
  tsim = matrix(NA,dim(param)[1],nrow(t))
  for(i in 1:dim(param)[1]){
    tsim[i,] = X%*%log(p[i,])
  }
  
  
  
  pdf(paste(folderN,"/residuals_",N,".pdf",sep=""), width=25, height=10)
  
  maxi = max(c(as.vector(t),as.vector(tsim)), na.rm = T)
  mini = min(c(as.vector(t),as.vector(tsim)), na.rm = T)
  # boxplot(tsim,outline=FALSE,ylim=c(mini,maxi),
  #         lty = "blank", whisklty="blank", medlty = "solid")
  # 
  # sapply(c(1:numRep), function(r){
  #   points(c(1:nrow(t)),t[,r],col="red", cex = .4)
  # })
  
  PP = chain[-(1:burnIn),]
  colnames(PP) = paramNames
  boxplot(PP)
  # sapply(c(1:numRep), function(r){
  #   points(c(1:nrow(scores)),scores[,r],col="red", cex = .4)
  # })
  
  par(mfrow = c(1,numRep))
  sapply(c(1:numRep), function(r){
    rr = which(!is.na(t[,r]))
    pcc = cor(x = t[rr,r], y = apply(tsim,2,mean,na.rm = T)[rr], method = "pearson")
    plot(t[,r],apply(tsim,2,mean),pch = 16,
         xlab="cleavage/splicing strength",ylab="simulation",
         main = paste0("bio rep: ", r, " - PCC = ", round(pcc,4)),
         ylim = c(mini,maxi),
         xlim = c(mini,maxi),
         axes=FALSE)
    axis(1)
    axis(2,las=2)
    abline(a=0,b=1)
    
  })
  
  dev.off()
  
}

# ----- likelihood function -----
# overwrites the function in functionsInference.R
likelihoodFun <- function(param){
  
  sigma = param[length(param)]  # sd of prior
  p = param[-length(param)]  # parameters to infer
  
  tsim = X%*%log(p)
  
  
  SD = sigma
  likelihood = sapply(1:numRep, function(r){
    rr = which(!is.na(t[,r]))
    return(sum(dnorm(x=tsim[rr,],mean=t[rr,r],sd=SD,log=TRUE)))
  }) %>%
    sum(na.rm = T)
  
  if(is.na(likelihood) | !is.finite(likelihood)){
    likelihood = -10**11
  }
  
  
  k = which(tsim>1)
  if(length(k)>0){
    likelihood = likelihood-length(k)*10000
  }
  
  return(likelihood)
}

# ----- initiate MCMC -----

# define prior and likelihood
prior <- createUniformPrior(mini, maxi)
bayesianSetup <- createBayesianSetup(likelihood = likelihoodFun, prior = prior)


# initialize and run sampler
settings <- list(iterations = Niter,
                 consoleUpdates = 5000,
                 nrChains = 1,
                 Z = NULL,  # starting population
                 startValue=5,  # number of chains
                 pSnooker = 1e-06,  # probability of Snooker update
                 burnin = 0,  # number of iterations as burn in (not recorded)
                 thin = 10,  # thinning parameter
                 f = 2.38,  # scaling factor gamma
                 eps = 1e-06,  # small number to avoid singularity
                 pGamma1 = 0.1,  # probability determining the frequency with which the scaling is set to 1 
                 eps.mult = 2,  # random term (multiplicative error)
                 eps.add = 0.0,  # random term
                 zUpdateFrequency = 1,
                 currentChain = 1,
                 message = TRUE)


print("START")
out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DEzs", settings = settings)
print("END")


posterior = getSample(out,end = NULL, thin = 10,parametersOnly=TRUE, whichParameters = 1:numParam)
plotChain(posterior)


for(ii in 1:100){
  
  print(ii)
  print("RESTART")
  out <- runMCMC(out, sampler = "DEzs", settings = settings)
  print("END")
  posterior = getSample(out,end = NULL, thin = 10,parametersOnly=TRUE, whichParameters = 1:numParam)
  plotChain(posterior)
  save(posterior,file=paste(folderN,"/posterior.RData",sep=""))
  
}

save(out, file=paste(folderN,"/out.RData",sep=""))

