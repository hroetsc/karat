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

# hyperparameters

suppressWarnings(dir.create("results/Bayesian_ProteaSMM/"))
folderN = "results/Bayesian_ProteaSMM/"

### INPUT ###
inpFolder = "data/ProteaSMM/PCP_SR1feat_P1_zscaled/"
load(paste0(inpFolder,"DATA.RData"))
X = DATA$X

# t = DATA$t
t = rowMeans(DATA$t, na.rm = T) %>% as.matrix()
# t = log10(DATA$t)
# t[!is.finite(t)] = NA
ptheor = ginv(X)%*%t
ttheor = X%*%ptheor

# sample t based on actual distribution
S = list()
MSEs = rep(NA, 1e05)
for (iter in 1:1e05) {
  # scores = matrix(runif(ncol(X), min = -20, max = 20), nrow = ncol(X), ncol = 1)
  scores = ptheor+rnorm(n=ncol(X),mean = 0, sd = 1e-2)
  MSEs[iter] = mean((X%*%scores - t)^2)
  S[[iter]] = scores
}
k = which.min(MSEs)
MSEs[k]
S[[k]]

tsample = X%*%S[[k]]
tsample[tsample<0] = 0

# tround = round(t,digits = 1) %>% as.data.frame() %>% group_by(V1) %>% summarise(n = n())
# tsample = rep(0, tround$n[1])
# for (j in 2:nrow(tround)) {
#   cntsample = sample(x = seq(tround$V1[j-1]+1e-05,tround$V1[j],0.001), size = tround$n[j], replace = T)
#   tsample = c(tsample, cntsample+1e-05)
# }
# tsample = tsample + runif(length(tsample), min = 1e-06,1e-05)

density(t) %>% plot()
density(tsample) %>% lines()
# t = tsample[sample(1:length(tsample))] %>% as.matrix()
t = tsample

scores = S[[k]]

folderN = str_replace_all(inpFolder, "data/ProteaSMM/", "results/Bayesian_ProteaSMM/")
folderN = "results/Bayesian_ProteaSMM/PCP_SR1feat_P1_zscaled_proofRealDistr/"
suppressWarnings(dir.create(folderN))

save(t, file = paste0(folderN, "sampledTarget.RData"))
save(scores, file = paste0(folderN, "sampledScores.RData"))

numRep = ncol(t)
numParam = ncol(X)+1

### MAIN PART ###
# ---- settings -----
Niter = 3*10**5

# change prior according to models
# set prior: initial conditions followed by parameters
mini = rep(-10,numParam)  # 1 parameter more than required by model --> sigma
maxi = rep(20,numParam)
mini[length(mini)] = 0
maxi[length(maxi)] = 20  # initial sigma

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
  
  # par(mfrow = c(1,1))
  # xxx = chain[-(1:burnIn),]
  # colnames(xxx) = paramNames
  
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
    tsim[i,] = X%*%p[i,]
    tsim[i,which(tsim[i,]<0)] = 0
  }
  
  
  
  pdf(paste(folderN,"/residuals_",N,".pdf",sep=""), width=25, height=10)
  
  maxi = max(c(as.vector(t),as.vector(tsim)), na.rm = T)
  mini = min(c(as.vector(t),as.vector(tsim)), na.rm = T)
  boxplot(tsim,outline=FALSE,ylim=c(mini,maxi),
          lty = "blank", whisklty="blank", medlty = "solid")
  
  sapply(c(1:numRep), function(r){
    points(c(1:nrow(t)),t[,r],col="red", cex = .4)
  })
  
  PP = chain[-(1:burnIn),]
  colnames(PP) = paramNames
  boxplot(PP)
  sapply(c(1:numRep), function(r){
    points(c(1:nrow(scores)),scores[,r],col="red", cex = .4)
  })
  
  par(mfrow = c(1,numRep))
  sapply(c(1:numRep), function(r){
    plot(t[,r],apply(tsim,2,mean),pch = 16,
         xlab="cleavage/splicing strength",ylab="simulation",
         main = paste0("bio rep: ", r),
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
  
  # noise = matrix(rnorm(nrow(X)*ncol(X), mean = 0, sd = 1e-05), nrow = nrow(X), ncol = ncol(X))
  tsim = X%*%p
  tsim[which(tsim<0)] = 0
  
  # # compute likelihood
  # likelihoods = rep(NA, numRep)
  # for (r in 1:numRep) {
  #   rr = which(!is.na(t[,r]))
  #   likelihood = sum(dnorm(x=tsim[rr,],mean=t[rr,r],sd=sigma,log=TRUE))
  # 
  #   # replace missing values with strongly negative ones
  #   # print(likelihood)
  # 
  #   if(is.na(likelihood) | !is.finite(likelihood)){
  #     likelihood = -10**11
  #   }
  # 
  #   likelihoods[r] = likelihood
  # }
  # likelihood = sum(likelihoods, na.rm = T)
  
  
  SD = sigma
  likelihood = sum(dnorm(x=tsim,mean=t,sd=SD,log=TRUE))

  if(is.na(likelihood) | !is.finite(likelihood)){
    likelihood = -10**11
  }


  k = which(tsim<0 | tsim>100)
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
                 consoleUpdates = 10000,
                 nrChains = 1,
                 Z = NULL,  # starting population
                 startValue=3,  # number of chains
                 pSnooker = 0.01,  # probability of Snooker update
                 burnin = 0,  # number of iterations as burn in (not recorded)
                 thin = 10,  # thinning parameter
                 f = 2.38,  # scaling factor gamma
                 eps = 0.01,  # small number to avoid singularity
                 parallel = NULL,
                 pGamma1 = 0.1,  # probability determining the frequency with which the scaling is set to 1 
                 eps.mult = 2,  # random term (multiplicative error)
                 eps.add = 0.5,  # random term
                 zUpdateFrequency = 1,
                 currentChain = 1,
                 message = TRUE)


print("START")
out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DEzs", settings = settings)
print("END")


posterior = getSample(out,end = NULL, thin = 10,parametersOnly=TRUE, whichParameters = 1:length(mini))
plotChain(posterior)


for(ii in 1:100){
  
  print(ii)
  print("RESTART")
  out <- runMCMC(out, sampler = "DEzs", settings = settings)
  print("END")
  posterior = getSample(out,end = NULL, thin = 10,parametersOnly=TRUE, whichParameters = 1:length(mini))
  plotChain(posterior)
  save(posterior,file=paste(folderN,"/posterior.RData",sep=""))
  
}


