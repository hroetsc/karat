### karat projetc - PCPS mechanism ###
# description:  analyse posterior distributions of Bayesian ProteaSMM
# input:        BayesianTools: posterior.RData, PCP dataset
# output:       posterior distributions across CV folds
# author:       HPR

library(dplyr)
library(stringr)
library(ggplot2)
library(BayesianTools)
library(numDeriv)
library(corrplot)
library(RColorBrewer)
library(transport)
library(tidymodels)
library(DescTools)
library(eulerr)
library(entropy)
library(dbscan)
source("src/invitroSPI_utils.R")
source("src/plotting_utils.R")
source("src/_ROC-curve.R")

theme_set(theme_classic())
AAchar_here = c("P","G","C","M","A","V","I","L","F","Y","W","H","R","K","D","E","N","Q","S","T","X")
AAchar_here_sorted = sort(AAchar_here)

rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
rcol <- rf(100)


### INPUT ###
# fs = list.files("Bayesian_ProteaSMM/server/singleSubs_PCP_0621/", pattern = "posterior.RData", recursive = T, full.names = T)
# load("data/ProteaSMM/PCP_SR1extfeat_P1/DATA.RData")

fs = list.files("Bayesian_ProteaSMM/server/singleSubs_SR1_0623/", pattern = "posterior.RData", recursive = T, full.names = T)
load("data/ProteaSMM/PSP_SR1extfeat_P1/DATA.RData")


### MAIN PART ###
suppressWarnings(dir.create("results/Bayesian_ProteaSMM/PLOTS/LOV/"))

# ----- load posterior and data -----
pseudo = 1e-05
X = DATA$X
t = log(DATA$t/100+pseudo)
t[!is.finite(t)] = NA

subIDs = DATA$substrateIDs
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

save(CHAINS, file = "results/Bayesian_ProteaSMM/PLOTS/LOV/PSP_chains_0623.RData")

# ----- plot posterior distributions -----
allParamLim = list()
for (j in 2:ncol(CHAINSdf)) {
  
  cnt = CHAINSdf[,c(1,j)]
  names(cnt) = c("file","param")
  cntP = ggplot(cnt, aes(x = param, col = file)) +
    geom_density() + 
    xlab("posterior value") +
    ggtitle(colnames(CHAINSdf)[j]) +
    scale_color_viridis_d() +
    xlim(c(0,2)) +
    theme(legend.position = "none")
  
  allParamLim[[j-1]] = cntP
}


ggsave(filename = "results/Bayesian_ProteaSMM/PLOTS/LOV/0623_PSPposteriors_priorRange.pdf", 
       plot = gridExtra::marrangeGrob(allParamLim, nrow=5, ncol=5, byrow = T), 
       width = 20, height = 20, dpi = "retina")

# ----- prediction on left-out data set -----
zscale = function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}

left_out = str_extract_all(fs, pattern = "[:alnum:]+(?=/posterior.RData)", simplify = T) %>% as.vector()

pdf("results/Bayesian_ProteaSMM/PLOTS/LOV/0623_PSPperformance.pdf", height = 5*4, width = 3*4)
par(mfrow = c(5,3))

PCC = rep(NA, length(CHAINS))
allROCs = list()

for (i in 1:length(CHAINS)) {
  
  print(left_out[i])
  keep = which(subIDs == left_out[i])
  
  # get prediction
  cntPosterior = CHAINS[[i]]
  cntPosterior = cntPosterior[,colnames(X)]
  tsims = apply(cntPosterior,1,function(p){
    tsim = X[keep,]%*%log(p)
    return(tsim)
  })
  
  # tsims_mean = apply((exp(tsims)-pseudo)*100,1,mean,na.rm = T)
  # tsims_sd = apply((exp(tsims)-pseudo)*100,1,sd,na.rm = T)
  # 
  # ttrue_mean = apply((exp(t[keep, ])-pseudo)*100,1,mean,na.rm = T)
  # ttrue_sd = apply((exp(t[keep, ])-pseudo)*100,1,sd,na.rm = T)
  
  tsims_mean = apply(tsims,1,mean,na.rm = T)
  tsims_sd = apply(tsims,1,sd,na.rm = T)
  
  ttrue_mean = apply(t[keep, ],1,mean,na.rm = T)
  ttrue_sd = apply(t[keep, ],1,sd,na.rm = T)
  
  
  # correlation coefficient
  pcc = cor(ttrue_mean, tsims_mean)
  
  plot(x = ttrue_mean, y = tsims_mean,
       pch = 16, xlab = "true", ylab = "predicted",
       main = paste0(left_out[i], ", PCC = ", round(pcc,4)))
  
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
  allROCs[[i]] = getROCcurve(zscale(ttrue_mean), zscale(tsims_mean), substrate = left_out[i], threshPerc = 0.5)
}

barplot(PCC)

dev.off()


pdf("results/Bayesian_ProteaSMM/PLOTS/LOV/0623_PSPbinaryperformance.pdf", height = 4, width = 8)
lapply(allROCs, function(cnt){
  gridExtra::grid.arrange(cnt[[1]], cnt[[2]], ncol = 2)
})
dev.off()

# ----- prediction of joint posterior on left out substrate -----
jointPosterior = as.matrix(CHAINSdf[,paramNames])

keep = which(subIDs == "MM582")
tsims = apply(jointPosterior,1,function(p){
  tsim = X[keep,]%*%log(p)
  return(tsim)
})

tsims_mean = apply(tsims,1,mean,na.rm = T)
tsims_sd = apply(tsims,1,sd,na.rm = T)

ttrue_mean = apply(t[keep, ],1,mean,na.rm = T)
ttrue_sd = apply(t[keep, ],1,sd,na.rm = T)


pdf("results/Bayesian_ProteaSMM/PLOTS/LOV/0623_PSP_MM582.pdf", height = 5, width = 8)
# correlation coefficient
pcc = cor(ttrue_mean, tsims_mean)

plot(x = ttrue_mean, y = tsims_mean,
     pch = 16, xlab = "true", ylab = "predicted",
     main = paste0("MM582, PCC = ", round(pcc,4)))

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

roc = getROCcurve(zscale(ttrue_mean), zscale(tsims_mean), substrate = "MM582", threshPerc = 0.5)
gridExtra::grid.arrange(roc[[1]], roc[[2]], ncol = 2)
dev.off()


# ----- PCA on joint posterior -----

# corM = cor(jointPosterior)
covM = cov(jointPosterior)

PCs = eigen(covM)

loadings = sapply(1:length(paramNames), function(j){
  PCs$vectors[,j] * sqrt(PCs$values)[j]
})
rownames(loadings) = paramNames

barplot(PCs$values, main = "eigenvalues", xlab = "principal components")

# ----- plot loadings -----
# colour by amino acid
aminoAcids = str_split_fixed(paramNames, pattern = ";", Inf)[,2] %>% as.data.frame()
AAcolors = data.frame(. = AAchar_here,
                      cols = c(rep("orange",2), rep("pink",2), rep("lightblue", 4), rep("gray", 3), rep("darkblue", 3), rep("firebrick", 4), rep("purple", 2), "black"))
AAcolors = left_join(aminoAcids, AAcolors)

# colour by position
positions = str_split_fixed(paramNames, pattern = ";", Inf)[,1] %>% as.data.frame()
Poscolors = data.frame(. = c("P6","P5","P4", "P3", "P2", "P1",
                             "P-1","P-2","P-3","P-4","P-5","P-6"),
                      cols = c(c("darkgoldenrod","blueviolet","lightgreen","orangered","lightblue","gray"),
                               rev(c("darkgoldenrod","blueviolet","lightgreen","orangered","lightblue","gray"))))
Poscolors = left_join(positions, Poscolors)


pdf("results/Bayesian_ProteaSMM/PLOTS/LOV/0623_PSP-PCAslopiness.pdf", height = 8, width = 18)

# ----- sloppy parameters
par(mfrow = c(2,1))
# colour by amino acid
k = order(loadings[,1], decreasing = T)
cnt = loadings[k,1][c(1:40, 210:250)]
plot(x = c(1:81), y = cnt, type = "h",
     col = AAcolors$cols[k][c(1:40, 210:250)], ylab = "PCA loadings", xlab = "parameters", main = "sloppy parameters",
     bty = "n", lwd = 5, frame.plot = F, axes = F)
text(paramNames[k][c(1:40, 210:250)],
     x = c(1:81), y = cnt, cex = .5, srt=90)
axis(2)

# colour by position
plot(x = c(1:81), y = cnt, type = "h",
     col = Poscolors$cols[k][c(1:40, 210:250)], ylab = "PCA loadings", xlab = "parameters", main = "sloppy parameters",
     bty = "n", lwd = 5, frame.plot = F, axes = F)
text(paramNames[k][c(1:40, 210:250)],
     x = c(1:81), y = cnt, cex = .5, srt = 90)
axis(2)

# ----- stiff parameters
par(mfrow = c(2,1))
# colour by amino acid
k = order(loadings[,250], decreasing = T)
cnt = loadings[k,250][c(1:40, 210:250)]
plot(x = c(1:81), y = cnt, type = "h",
     col = AAcolors$cols[k][c(1:40, 210:250)], ylab = "PCA loadings", xlab = "parameters", main = "stiff parameters",
     bty = "n", lwd = 5, frame.plot = F, axes = F)
text(paramNames[k][c(1:40, 210:250)],
     x = c(1:81), y = cnt, cex = .5, srt=90)
axis(2)

# colour by position
plot(x = c(1:81), y = cnt, type = "h",
     col = Poscolors$cols[k][c(1:40, 210:250)], ylab = "PCA loadings", xlab = "parameters", main = "stiff parameters",
     bty = "n", lwd = 5, frame.plot = F, axes = F)
text(paramNames[k][c(1:40, 210:250)],
     x = c(1:81), y = cnt, cex = .5, srt = 90)
axis(2)

dev.off()


# ----- contribution of parameters/positions to certain PCs -----

pdf("results/Bayesian_ProteaSMM/PLOTS/LOV/0623_PSP-PCAloadings.pdf", height = 12, width = 24)
positions = c("P6","P5","P4", "P3", "P2", "P1","P-1","P-2","P-3","P-4","P-5","P-6")
pcs = paste0("PC",seq(1,ncol(PCs$vectors)))

eigenvectors = loadings %>% as.data.frame()
colnames(eigenvectors) = pcs
eigenvectors$parameters = paramNames
eigenvectors$positions = str_split_fixed(paramNames,";",Inf)[,1]
eigenvectors$aas = str_split_fixed(paramNames,";",Inf)[,2]

eigenvectors = eigenvectors %>% tidyr::gather(key, value, -parameters, -positions, -aas)

# --- contribution of positions 
posGroup = eigenvectors %>% group_by(key, positions) %>%
  summarise(contribution = sum(value))
posGroup = posGroup %>% tidyr::spread(key, contribution) %>% as.data.frame()
posGroup = posGroup[match(posGroup$positions, positions), pcs] %>% as.matrix()

posGroup = apply(posGroup,2,function(x){
  (x - min(x)) / (max(x) - min(x))
})
image(posGroup %>% t(), axes = F, col = rcol,
      xlab = "principal component",
      ylab = "position",
      sub = "column-wise z-scaled PCA loadings - 0: blue, 1: red")
axis(2, at = seq(0,1,1/(length(positions)-1)), labels = positions)
axis(1, at = seq(0,1,1/(length(pcs)/10)), labels = c(pcs[seq(1,250,10)], "PC250"))

# # cluster based on first and last PC contributions
# cl = hdbscan(posGroup[,paste0("PC",seq(1,12))], minPts = 2)
# plot(cl$hc,
#      main = "sloppiness",
#      xlab = "HDBSCAN* Hierarchy",
#      labels = positions)
# cl = hdbscan(posGroup[,paste0("PC",seq(239,250))], minPts = 2)
# plot(cl$hc,
#      main = "stiffness",
#      xlab = "HDBSCAN* Hierarchy",
#      labels = positions)


# --- contribution of amino acids 
aaGroup = eigenvectors %>% group_by(key, aas) %>%
  summarise(contribution = sum(value))
aaGroup = aaGroup %>% tidyr::spread(key, contribution) %>% as.data.frame()
aaGroup = aaGroup[match(aaGroup$aas, AAchar_here), pcs] %>% as.matrix()

aaGroup = apply(aaGroup,2,function(x){
  (x - min(x)) / (max(x) - min(x))
})
image(aaGroup %>% t(), axes = F, col = rcol,
      xlab = "principal component",
      ylab = "amino acid",
      sub = "column-wise z-scaled PCA loadings - 0: blue, 1: red")
axis(2, at = seq(0,1,1/(length(AAchar_here)-1)), labels = AAchar_here)
axis(1, at = seq(0,1,1/(length(pcs)/10)), labels = c(pcs[seq(1,250,10)], "PC250"))

# # cluster based on first and last PC contributions
# cl = hdbscan(aaGroup[,paste0("PC",seq(1,12))], minPts = 2)
# plot(cl$hc,
#      main = "sloppiness",
#      xlab = "HDBSCAN* Hierarchy",
#      labels = AAchar_here)
# cl = hdbscan(aaGroup[,paste0("PC",seq(239,250))], minPts = 2)
# plot(cl$hc,
#      main = "stiffness",
#      xlab = "HDBSCAN* Hierarchy",
#      labels = AAchar_here)
dev.off()

# ----- distributions of stiff and sloppy parameters -----
numParam = length(paramNames)

firstLoadings = as.vector(abs(loadings[,1:12]))
names(firstLoadings) = rep(paramNames, 12)
sloppyParams = which(firstLoadings > quantile(firstLoadings, 0.9)) %>% names() %>% unique()

lastLoadings = as.vector(abs(loadings[,(numParam-12+1):numParam]))
names(lastLoadings) = rep(paramNames, 12)
stiffParams = which(lastLoadings > quantile(lastLoadings, 0.9)) %>% names() %>% unique()


# overlap between stiff and sloppy parameters
euler(list(stiff = stiffParams,
           sloppy = sloppyParams), shape = "ellipse") %>% plot(quantities = T)
stiffParams[stiffParams %in% sloppyParams]


# get parameter categories
jointPosteriorDF = jointPosterior %>% as.data.frame() %>% tidyr::gather()
head(jointPosteriorDF)

jointPosteriorDF$position = str_split_fixed(jointPosteriorDF$key,";",Inf)[,1]
jointPosteriorDF$aa = str_split_fixed(jointPosteriorDF$key,";",Inf)[,2]
jointPosteriorDF$aa = factor(jointPosteriorDF$aa, levels = AAchar_here)

jointPosteriorDF = jointPosteriorDF %>%
  mutate(class = ifelse(key %in% stiffParams, "stiff", "none"),
         class = ifelse(key %in% sloppyParams, "sloppy", class))


# ----- calculate entropy -----

# discretise marginal posteriors to get same bins everywhere
bins = seq(0,2,length.out = 1e03)

entropies = rep(NA, ncol(jointPosterior))
for (j in 1:ncol(jointPosterior)) {
  print(colnames(jointPosterior)[j])
  cntPosterior = jointPosterior[,j] %>% sort()
  y = sapply(bins, function(x){
    length(which(cntPosterior <= x))
  })
  entropies[j] = entropy(diff(y))
}

names(entropies) = colnames(jointPosterior)
barplot(entropies)

# get entropy of the prior distribution
prior = runif(n = 2e06, min = 0, max = 2)
y = sapply(bins, function(x){
  length(which(prior <= x))
})
priorEntropy = entropy(diff(y))

# difference between prior and posterior entropy
DeltaH = priorEntropy - entropies
barplot(DeltaH)

informativeParams = names(DeltaH)[DeltaH > quantile(DeltaH, 0.7)]
uninformativeParams = names(DeltaH)[DeltaH < quantile(DeltaH, 0.3)]


# euler diagram
png("results/Bayesian_ProteaSMM/PLOTS/LOV/0623_PSP_meaningfulParams.png", height = 4, width = 5, units = "in", res = 300)
euler(list(stiff = stiffParams, sloppy = sloppyParams,
           informative = informativeParams, uninformative = uninformativeParams),
      shape = "ellipse") %>%
  plot(quantities = T) %>%
  print()
dev.off()


# plot
jointPosteriorDF = jointPosteriorDF %>%
  mutate(infoclass = ifelse(key %in% informativeParams, "informative", "none"),
         infoclass = ifelse(key %in% uninformativeParams, "uninformative", infoclass))

save(DeltaH, file = "results/Bayesian_ProteaSMM/PLOTS/LOV/0622_PCP_DeltaH.RData")

# ----- plot parameter distributions -----
jointPosteriorDF$class = factor(jointPosteriorDF$class, levels = c("stiff", "sloppy", "none"))
jointPosteriorDF$infoclass = factor(jointPosteriorDF$infoclass, levels = c("informative", "uninformative", "none"))

fillVals = c(stiff = "palevioletred", sloppy = "paleturquoise", none = "gray")
colVals = c(informative = "violetred4", uninformative = "turquoise4", none = "black")

# plot per position
pos = c("P6","P5","P4", "P3", "P2", "P1", "P-6", "P-5", "P-4", "P-3","P-2", "P-1")
allP = list()
counter = 1
for (p in 1:length(pos)) {
  cntR = jointPosteriorDF[jointPosteriorDF$position == as.character(pos[p]), ]
  
  if (nrow(cntR) > 0) {
    cntP = ggplot(cntR, aes(x = aa, y = value, fill = class, color = infoclass)) +
      geom_boxplot(outlier.shape = NA, coef = 0, position = position_dodge(0.5), aes(alpha = .9), lwd = .8) +
      scale_fill_manual("PCA analysis", values = fillVals[names(fillVals) %in% levels(cntR$class)]) +
      scale_color_manual("entropy analysis", values = colVals[names(colVals) %in% levels(cntR$infoclass)]) +
      xlab("amino acid") +
      ylab("posterior parameter distribution") +
      ggtitle(pos[p])
    
    allP[[counter]] = cntP
    counter = counter + 1
  }
  
  
}

ggsave(filename = "results/Bayesian_ProteaSMM/PLOTS/LOV/0623_PSP_PARAMETERS.pdf", 
       plot = gridExtra::marrangeGrob(allP, nrow=6, ncol=2, byrow = F), 
       width = 15, height = 27, dpi = "retina")



allPaa = list()
counter = 1
jointPosteriorDF$position = factor(jointPosteriorDF$position, levels = c("P6","P5","P4", "P3", "P2", "P1",
                                                                         "P-1","P-2","P-3","P-4","P-5","P-6"))
for (a in 1:length(AAchar_here)) {
  
  cntR = jointPosteriorDF[jointPosteriorDF$aa == as.character(AAchar_here[a]), ]
  
  if (nrow(cntR) > 0) {
    cntP = ggplot(cntR, aes(x = position, y = value, fill = class, color = infoclass)) +
      geom_boxplot(outlier.shape = NA, coef = 0, position = position_dodge(0.5), aes(alpha = .9), lwd = .8) +
      scale_fill_manual("PCA analysis", values = fillVals[names(fillVals) %in% levels(cntR$class)]) +
      scale_color_manual("entropy analysis", values = colVals[names(colVals) %in% levels(cntR$infoclass)]) +
      xlab("position") +
      ylab("posterior parameter distribution") +
      ggtitle(AAchar_here[a])
    
    allPaa[[counter]] = cntP
    counter = counter + 1
  }
  
}

ggsave(filename = "results/Bayesian_ProteaSMM/PLOTS/LOV/0623_PSP_PARAMETERS_aawise.pdf", 
       plot = gridExtra::marrangeGrob(allPaa, nrow=7, ncol=3, byrow = T), 
       width = 28, height = 28, dpi = "retina")



# plot marginal posterior distributions of the overlap
stiffAndInformative = intersect(stiffParams, informativeParams)
sloppyAndUnInformative = intersect(sloppyParams, uninformativeParams)

stiff_info = CHAINSdf[,c(".id", stiffAndInformative)]
sloppy_uninfo = CHAINSdf[,c(".id", sloppyAndUnInformative)]


StiffInfoParam = list()
for (j in 2:ncol(stiff_info)) {
  
  cnt = stiff_info[,c(1,j)]
  names(cnt) = c("file","param")
  cntP = ggplot(cnt, aes(x = param, col = file)) +
    geom_density() + 
    xlab("posterior value") +
    ggtitle(colnames(stiff_info)[j]) +
    scale_color_viridis_d() +
    xlim(c(0,2)) +
    theme(legend.position = "none")
  
  StiffInfoParam[[j-1]] = cntP
}


ggsave(filename = "results/Bayesian_ProteaSMM/PLOTS/LOV/0623_STIFF+INFORMATIVE_PSPposteriors.pdf", 
       plot = gridExtra::marrangeGrob(StiffInfoParam, nrow=5, ncol=5, byrow = T), 
       width = 20, height = 20, dpi = "retina")



SloppyUninfoParam = list()
for (j in 2:ncol(sloppy_uninfo)) {
  
  cnt = sloppy_uninfo[,c(1,j)]
  names(cnt) = c("file","param")
  cntP = ggplot(cnt, aes(x = param, col = file)) +
    geom_density() + 
    xlab("posterior value") +
    ggtitle(colnames(sloppy_uninfo)[j]) +
    scale_color_viridis_d() +
    xlim(c(0,2)) +
    theme(legend.position = "none")
  
  SloppyUninfoParam[[j-1]] = cntP
}


ggsave(filename = "results/Bayesian_ProteaSMM/PLOTS/LOV/0623_SLOPPY+UNINFORMATIVE_PSPposteriors.pdf", 
       plot = gridExtra::marrangeGrob(SloppyUninfoParam, nrow=5, ncol=5, byrow = T), 
       width = 20, height = 20, dpi = "retina")


### OUTPUT ###
params = stiffAndInformative
save(params, file = "data/ProteaSMM/PSP_SR1extfeat_P1/stiff_informative_params.RData")
params = sloppyAndUnInformative
save(params, file = "data/ProteaSMM/PSP_SR1extfeat_P1/sloppy_uninformative_params.RData")

save(jointPosterior, file = "results/Bayesian_ProteaSMM/PLOTS/LOV/0623_PSPposteriors.RData")


