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
source("../_from-others/functionsInference.R")

# hyperparameters
suppressWarnings(dir.create("results/Bayesian_ProteaSMM/"))

### INPUT ###
inpFolder = "data/ProteaSMM/PSP_SR1feat_P1/"
load(paste0(inpFolder,"DATA.RData"))
X1 = DATA$X

t1 = rowMeans(DATA$t, na.rm = T) %>% as.matrix()
t1 = log(t1)
t1[!is.finite(t1)] = NA


inpFolder = "data/ProteaSMM/PSP_SR2feat_P1_/"
load(paste0(inpFolder,"DATA.RData"))
X2 = DATA$X

t2 = rowMeans(DATA$t, na.rm = T) %>% as.matrix()
t2 = log(t2)
t2[!is.finite(t2)] = NA

X = rbind(X1,X2)
t = rbind(t1,t2)
i1 = c(1:ncol(X1))
i2 = c((ncol(X1)+1):(ncol(X1)+ncol(X2)))

folderN = str_replace_all(inpFolder, "data/ProteaSMM/PSP_SR2feat_P1_/",
                          "results/Bayesian_ProteaSMM/PSP_SR1+2feat/")
suppressWarnings(dir.create(folderN))


numRep = ncol(t)
numParam = ncol(X1)+ncol(X2)+1

### MAIN PART ###
# ---- settings -----
Niter = 5*10**4

# change prior according to models
# set prior: initial conditions followed by parameters
mini = rep(-100,numParam)  # 1 parameter more than required by model --> sigma
maxi = rep(100,numParam)
maxi[length(maxi)] = 50  # initial sigma

# ----- plotting analytics -----

plotChain <- function(chain){
  
  # define parameters
  paramNames <- c(colnames(X1),colnames(X2),"sigma")
  
  N = dim(chain)[1]
  burnIn = round(0.7*N)
  
  
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
  tsim1 = matrix(NA,dim(param)[1],nrow(X1))
  tsim2 = matrix(NA,dim(param)[1],nrow(X2))
  for(i in 1:dim(param)[1]){
    tsim1[i,] = X1%*%p[i,i1]
    tsim2[i,] = X2%*%p[i,i2]
  }
  
  pdf(paste(folderN,"/residuals_",N,".pdf",sep=""), width=25, height=10)
  
  maxi = max(c(as.vector(t),as.vector(tsim1),as.vector(tsim2)), na.rm = T)
  mini = min(c(as.vector(t),as.vector(tsim1),as.vector(tsim2)), na.rm = T)
  
  boxplot(tsim1,outline=FALSE,ylim=c(mini,maxi),
          lty = "blank", whisklty="dashed", medlty = "solid")
  sapply(c(1:numRep), function(r){
    points(c(1:nrow(t1)),t1[,r],col="red", cex = .4)
  })
  
  boxplot(tsim2,outline=FALSE,ylim=c(mini,maxi),
          lty = "blank", whisklty="dashed", medlty = "solid")
  sapply(c(1:numRep), function(r){
    points(c(1:nrow(t2)),t2[,r],col="red", cex = .4)
  })
  
  PP = chain[-(1:burnIn),]
  colnames(PP) = paramNames
  boxplot(PP)
  
  par(mfrow = c(1,numRep))
  sapply(c(1:numRep), function(r){
    plot(t1[,r],apply(tsim1,2,mean),pch = 16,
         xlab="cleavage/splicing strength",ylab="simulation",
         main = paste0("bio rep: ", r),
         ylim = c(mini,maxi),
         xlim = c(mini,maxi),
         col = "gray",
         axes=FALSE)
    points(t2[,r],apply(tsim2,2,mean),pch = 16,col = "lightblue")
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
  
  tsim1 = X1%*%p[i1]
  tsim2 = X2%*%p[i2]
  
  # compute likelihood
  likelihoods = rep(NA, numRep)
  for (r in 1:numRep) {
    rr1 = which(!is.na(t1[,r]))
    likelihood1 = sum(dnorm(x=tsim1[rr1,],mean=t1[rr1,r],sd=sigma,log=TRUE))
    
    rr2 = which(!is.na(t2[,r]))
    likelihood2 = sum(dnorm(x=tsim2[rr2,],mean=t2[rr2,r],sd=sigma,log=TRUE))
    
    likelihoods[r] = likelihood1+likelihood2
  }
  # likelihoods = sapply(c(1:numRep), function(r){
  #   rr1 = which(!is.na(t1[,r]))
  #   likelihood1 = sum(dnorm(x=tsim1[rr1,],mean=t1[rr1,r],sd=sigma,log=TRUE))
  #   
  #   rr2 = which(!is.na(t2[,r]))
  #   likelihood2 = sum(dnorm(x=tsim2[rr2,],mean=t2[rr2,r],sd=sigma,log=TRUE))
  #   
  #   return(likelihood1+likelihood2)
  # })
  # 
  likelihood = sum(likelihoods, na.rm = T)
  return(likelihood)
}

# ----- initiate MCMC -----

# define prior and likelihood
prior <- createUniformPrior(mini, maxi)
bayesianSetup <- createBayesianSetup(likelihood = likelihoodFun, prior = prior)

# initialize and run sampler
settings <- list(iterations = Niter,
                 consoleUpdates = 1000,
                 nrChains = 1,
                 Z = NULL,  # starting population
                 startValue=5,  # number of chains
                 pSnooker = 0.001,  # probability of Snooker update
                 burnin = 0,  # number of iterations as burn in (not recorded)
                 thin = 10,  # thinning parameter
                 f = 2.38,  # scaling factor gamma
                 eps = 0.01,  # small number to avoid singularity
                 parallel = NULL,
                 pGamma1 = 0.01,  # probability determining the frequency with which the scaling is set to 1 
                 eps.mult = 0.5,  # random term (multiplicative error)
                 eps.add = 0.0,  # random term
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


