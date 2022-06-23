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
# ----- SR1
load(paste0(inpFolder,"DATA.RData"))
X = DATA$X
t1 = DATA$t1
t1 = log(t1/100+pseudo)
t1[!is.finite(t1)] = NA

k = which(DATA$substrateIDs != cntleaveOut)
X = X[k, ]
t1 = t1[k, ]

# ----- SR1
t2 = DATA$t2
t2 = log(t2/100+pseudo)
t2[!is.finite(t2)] = NA

t2 = t2[k, ]


# concatenate matrices and targets
i1 = c(1:ncol(X1))
i2 = c((ncol(X1)+1):(ncol(X1)+ncol(X2)))

numRep = ncol(t)
numParam = ncol(X1)+ncol(X2)+1

### MAIN PART ###
# ---- settings -----
Niter = 4*10**5

# change prior according to models
# set prior: initial conditions followed by parameters
mini = rep(0,numParam)  # 1 parameter more than required by model --> sigma
maxi = rep(2,numParam)
mini[length(mini)] = 0
maxi[length(maxi)] = 10 # initial sigma

# ----- plotting analytics -----

plotChain <- function(chain){
  
  # define parameters
  paramNames <- c(colnames(X1),colnames(X2),"sigma")
  
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
  tsim1 = matrix(NA,dim(param)[1],nrow(X1))
  tsim2 = matrix(NA,dim(param)[1],nrow(X2))
  for(i in 1:dim(param)[1]){
    tsim1[i,] = X1%*%log(p[i,i1])
    tsim2[i,] = X2%*%log(p[i,i2])
  }
  
  pdf(paste(folderN,"/residuals_",N,".pdf",sep=""), width=25, height=10)
  
  maxi = max(c(as.vector(t),as.vector(tsim1),as.vector(tsim2)), na.rm = T)
  mini = min(c(as.vector(t),as.vector(tsim1),as.vector(tsim2)), na.rm = T)
  
  boxplot(tsim1,outline=FALSE,ylim=c(mini,maxi),
          lty = "blank", whisklty="blank", medlty = "solid")
  sapply(c(1:numRep), function(r){
    points(c(1:nrow(t1)),t1[,r],col="red", cex = .4)
  })
  
  boxplot(tsim2,outline=FALSE,ylim=c(mini,maxi),
          lty = "blank", whisklty="blank", medlty = "solid")
  sapply(c(1:numRep), function(r){
    points(c(1:nrow(t2)),t2[,r],col="red", cex = .4)
  })
  
  PP = chain[-(1:burnIn),]
  colnames(PP) = paramNames
  boxplot(PP)
  
  par(mfrow = c(1,numRep))
  sapply(c(1:numRep), function(r){
    rr1 = which(!is.na(t1[,r]))
    pcc1 = cor(x = t1[rr1,r], y = apply(tsim1,2,mean,na.rm = T)[rr1], method = "pearson")
    
    rr2 = which(!is.na(t2[,r]))
    pcc2 = cor(x = t2[rr2,r], y = apply(tsim2,2,mean,na.rm = T)[rr2], method = "pearson")
    
    plot(t1[,r],apply(tsim1,2,mean),pch = 16,
         xlab="cleavage/splicing strength",ylab="simulation",
         main = paste0("bio rep: ", r, " - PCC = ", round(pcc1,4), ", ", round(pcc2,4)),
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
  
  tsim1 = X1%*%log(p[i1])
  tsim2 = X2%*%log(p[i2])
  
  # compute likelihood
  SD = sigma
  likelihood = sapply(1:numRep, function(r){
    rr1 = which(!is.na(t1[,r]))
    l1 = sum(dnorm(x=tsim1[rr1,],mean=t1[rr1,r],sd=SD,log=TRUE))
    
    rr2 = which(!is.na(t2[,r]))
    l2 = sum(dnorm(x=tsim2[rr2,],mean=t2[rr2,r],sd=SD,log=TRUE))
    return(l1+l2)
  }) %>%
    sum(na.rm = T)
  
  if(is.na(likelihood) | !is.finite(likelihood)){
    likelihood = -10**11
  }
  
  return(likelihood)
}

# ----- initiate MCMC -----
# define prior and likelihood
prior <- createUniformPrior(mini, maxi)
bayesianSetup <- createBayesianSetup(likelihood = likelihoodFun, prior = prior)

folderN = "results/Bayesian_ProteaSMM/PSPpoly_0614/"
suppressWarnings(dir.create(folderN))

# initialize and run sampler
settings <- list(iterations = Niter,
                 consoleUpdates = 5000,
                 nrChains = 1,
                 Z = NULL,  # starting population
                 startValue=5,  # number of chains
                 pSnooker = 1e-06,  # probability of Snooker update
                 burnin = 0,  # number of iterations as burn in (not recorded)
                 thin = 5,  # thinning parameter
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

save(out, file=paste(folderN,"/out.RData",sep=""))
