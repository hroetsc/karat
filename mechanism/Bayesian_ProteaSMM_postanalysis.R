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
source("src/invitroSPI_utils.R")

theme_set(theme_classic())
AAchar_here = c("P","G","C","M","A","V","L","F","Y","W","H","R","K","D","E","N","Q","S","T","X")

### INPUT ###

inpFolder = "data/ProteaSMM/PCP_SR1feat_P1/"
load(paste0(inpFolder,"DATA.RData"))
X = DATA$X

# t = DATA$t
t = rowMeans(DATA$t, na.rm = T) %>% as.matrix()
t = log(t/100)

k = grep("MM", DATA$substrateIDs)
X = X[k,]
t = t[k,] %>% as.matrix()


load("Bayesian_ProteaSMM/server/Bayesian_ProteaSMM/PCP_SR1feat_P1_logModel_server6/posterior.RData")
N = dim(posterior)[1]
burnIn = round(0.9*N)
chain = posterior[-(1:burnIn),]
colnames(chain) = c(colnames(X),"sigma")

### MAIN PART ###
suppressWarnings(dir.create("results/Bayesian_ProteaSMM/PLOTS/"))
suppressWarnings(dir.create("results/Bayesian_ProteaSMM/PLOTS/PCP/"))

# ----- correlation between true and predicted -----
NN = min(c(dim(chain)[1],100))

sampledIndex = sample(c(1:dim(chain)[1]),size=NN,replace=FALSE)
param = chain[sampledIndex,]
p = param[,-dim(param)[2]]

tsim = matrix(NA,dim(param)[1],nrow(t))
for(i in 1:dim(param)[1]){
  tsim[i,] = X%*%log(p[i,])
}

tsim_mean = apply(tsim,2,mean)
pcc = cor(x = t, y = tsim_mean, method = "pearson")

plot(x=t, y=tsim_mean, pch = 16,
     xlim = c(-14,0), ylim = c(-14,0),
     xlab = "log SCS-P1", ylab = "simulated log SCS-P1",
     main = paste0("PCC = ", round(pcc,4)))
abline(a=0,b=1,col = "red")

# ----- boxplots -----
PCPparam = tidyr::gather(chain %>% as.data.frame()) %>%
  mutate(reactant = "PCP")

PCPparam$position = str_split_fixed(PCPparam$key,";",Inf)[,1]
PCPparam$aa = str_split_fixed(PCPparam$key,";",Inf)[,2]

PCPparam = PCPparam[-which(PCPparam$position == "sigma"), ]

# reorder amino acids
PCPparam$aa = factor(PCPparam$aa, levels = AAchar_here)


# plot per position
pos = c("P4", "P3", "P2", "P1","P-4","P-3","P-2","P-1")
allP = list()

for (p in 1:length(pos)) {
  cntR = PCPparam[PCPparam$position == as.character(pos[p]), ]
  cntP = ggplot(cntR, aes(x = aa, y = value, fill = reactant)) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(0.5), aes(alpha = .8)) +
    geom_hline(yintercept = c(mean(cntR$value[cntR$reactant == "PCP"])),
               col = c(plottingCols[["PCP"]]), lty = "dashed", lwd = 0.8) +
    # ylim(c(-6,6)) +
    scale_fill_manual(values = c(plottingCols[["PCP"]])) +
    xlab("amino acid") +
    ylab("posterior distribution") +
    ggtitle(pos[p])
  
  allP[[p]] = cntP
  
}

ggsave(filename = "results/Bayesian_ProteaSMM/PLOTS/PCP/_paramDistributions.pdf", 
       plot = gridExtra::marrangeGrob(allP, nrow=4, ncol=2, byrow = F), 
       width = 15, height = 18, dpi = "retina")

# plot per amino acid
PCPparam$position = factor(PCPparam$position, levels = c("P4", "P3", "P2", "P1", "P-1", "P-2", "P-3", "P-4"))
allPaa = list()
for (a in 1:length(AAchar_here)) {
  
  cntR = PCPparam[PCPparam$aa == as.character(AAchar_here[a]), ]
  cntP = ggplot(cntR, aes(x = position, y = value, fill = reactant)) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(0.5), aes(alpha = .8)) +
    geom_hline(yintercept = c(mean(cntR$value[cntR$reactant == "PCP"])),
               col = c(plottingCols[["PCP"]]), lty = "dashed", lwd = 0.8) +
    # ylim(c(-6,6)) +
    scale_fill_manual(values = c(plottingCols[["PCP"]])) +
    xlab("position") +
    ylab("regression weight") +
    ggtitle(AAchar_here[a])
  
  allPaa[[a]] = cntP
}

ggsave(filename = "results/Bayesian_ProteaSMM/PLOTS/PCP/_paramDistributions_aa-wise.pdf", 
       plot = gridExtra::marrangeGrob(allPaa, nrow=7, ncol=3, byrow = T), 
       width = 22, height = 35, dpi = "retina")


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

