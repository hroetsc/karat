### karat projetc - PCPS mechanism ###
# description:  get predictive performance of model trained only on informative parameters
# input:        posterior.RData, true cleavage/splicing strengths
# output:       performance of informative SCS-P1 predictor
# author:       HPR

library(dplyr)
library(stringr)
source("src/plotting_utils.R")

pseudo = 1e-05

### INPUT ###
load("data/ProteaSMM/PSP_SR1extfeat_P1/DATA.RData")
load("Bayesian_ProteaSMM/server/SR1_0627/stiffAndInformative/posterior.RData")
load("data/ProteaSMM/PSP_SR1extfeat_P1/stiff_informative_params.RData")

### MAIN PART ###
X = DATA$X
t = log(DATA$t/100+pseudo)

X = X[, colnames(X) %in% params]

# ----- preprocessing -----
N = dim(posterior)[1]
burnIn = round(0.3*N)
chain = posterior[-(1:burnIn),]
colnames(chain) = c(colnames(X),"sigma")

parameters = chain[,colnames(X)]

# ----- prediction on left-out substrate -----
keep = which(DATA$substrateIDs == "MM582")
tsims = apply(parameters,1,function(p){
  tsim = X[keep,]%*%log(p)
  return(tsim)
})

tsims_mean = apply(tsims,1,mean,na.rm = T)
tsims_sd = apply(tsims,1,sd,na.rm = T)

ttrue_mean = apply(t[keep, ],1,mean,na.rm = T)
ttrue_sd = apply(t[keep, ],1,sd,na.rm = T)

pdf("results/Bayesian_ProteaSMM/PLOTS/LOV/StiffAndInformativePSPs_MM582.pdf", height = 6, width = 6)
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

dev.off()

save(parameters, file = "results/Bayesian_ProteaSMM/PLOTS/LOV/0623_PSPposteriors_stiff+informative.RData")


# ----- demonstrate fits -----
load("data/ProteaSMM/PSP_SR1extfeat_P1/DATA.RData")
load("data/ProteaSMM/PSP_SR1extfeat_P1/stiff_informative_params.RData")
load("Bayesian_ProteaSMM/server/SR1_0627/stiffAndInformative/posterior.RData")

# load("data/ProteaSMM/PSP_SR1extfeat_P1/DATA.RData")
# load("data/ProteaSMM/PSP_SR1extfeat_P1/sloppy_uninformative_params.RData")
# load("Bayesian_ProteaSMM/server/SR1_0627/stiffAndInformative/posterior.RData")


X = DATA$X
t = log(DATA$t/100+pseudo)
X = X[, colnames(X) %in% params]

N = dim(posterior)[1]
burnIn = round(0.3*N)
chain = posterior[-(1:burnIn),]
colnames(chain) = c(colnames(X),"sigma")

parameters = chain[,colnames(X)]


keep = which(DATA$substrateIDs != "MM582")
tsims = apply(parameters,1,function(p){
  X[keep,]%*%log(p)
})

tsims_mean = apply(tsims,1,mean,na.rm = T)
tsims_sd = apply(tsims,1,sd,na.rm = T)

ttrue_mean = apply(t[keep, ],1,mean,na.rm = T)
ttrue_sd = apply(t[keep, ],1,sd,na.rm = T)

# png("results/Bayesian_ProteaSMM/PLOTS/LOV/FIT_sloppy+uninformative_PCP.png", height = 5, width = 5, units = "in", res = 300)
png("results/Bayesian_ProteaSMM/PLOTS/LOV/FIT_stiff+informative_SR1.png", height = 5, width = 5, units = "in", res = 300)
lim = c(min(ttrue_mean, tsims_mean), max(ttrue_mean, tsims_mean))

# correlation coefficient
pcc = cor(ttrue_mean, tsims_mean)

plot(x = ttrue_mean, y = tsims_mean,
     pch = 16, xlab = "true", ylab = "predicted",
     col = add.alpha("black", .7), cex = .7,
     xlim = lim, ylim = lim,
     main = paste0("PCC = ", round(pcc,4)))

abline(a = 0, b = 1, col = "red")
# abline(v = log(0.1+pseudo), lty = "dashed", col = "blue")
# abline(h = log(0.1+pseudo), lty = "dashed", col = "blue")

arrows(x0 = ttrue_mean+ttrue_sd,
       x1 = ttrue_mean-ttrue_sd,
       y0 = tsims_mean, y1 = tsims_mean,
       col = add.alpha("black", .7),
       code = 3, angle = 90, length = 0.03, lwd = .3) %>% suppressWarnings()

arrows(x0 = ttrue_mean, x1 = ttrue_mean,
       y0 = tsims_mean+tsims_sd, y1 = tsims_mean-tsims_sd,
       col = add.alpha("black", .7),
       code = 3, angle = 90, length = 0.03, lwd = .3) %>% suppressWarnings()

dev.off()



### OUTPUT ###
