### karat projetc - PCPS mechanism ###
# description:  input for NetChop, validation of NetChop predictions for PCP
# input:        quantitative DBs, cleavage and splicing strengths, posteriors
# output:       substrate sequences as input for NetChop, NetChop performance
# author:       HPR


library(dplyr)
library(stringr)
library(seqinr)
source("src/invitroSPI_utils.R")
source("src/_ROC-curve.R")

pseudo = 1e-05

### INPUT ###
load("../../proteinsPCPS/new/data/aSPIre.RData")
proteins = Kinetics
load("data/aSPIre.RData")
polypeps = Kinetics

load("data/ProteaSMM/PCP_SR1extfeat_P1/DATA_proteinsOnly.RData")
# load("results/Bayesian_ProteaSMM/PLOTS/LOV/0622_PCPposteriors.RData")
load("results/Bayesian_ProteaSMM/PLOTS/LOV/0622_PCPposteriors_stiff+informative.RData")
load("data/ProteaSMM/PCP_SR1extfeat_P1/stiff_informative_params.RData")


### MAIN PART ###
suppressWarnings(dir.create("results/Bayesian_ProteaSMM/PREDICTION/NetChop-vs-BayesianProteaSMM/", recursive = T))

# ----- get substrates -----
nm = intersect(names(proteins), names(polypeps))
Quant = rbind(proteins[,nm], polypeps[,nm])

SeqTBL = Quant %>%
  group_by(substrateID, substrateSeq) %>%
  summarise() %>%
  unique()

seqList = lapply(SeqTBL$substrateSeq, function(x){x})
names(seqList) = SeqTBL$substrateID
write.fasta(sequences = seqList, names = SeqTBL$substrateID, file.out = "data/NetChop/sequences.fasta")


# ----- parser for NetChop output -----
file = "data/NetChop/allPredictions.txt"
readInNetChop = function(file) {
  
  tbl = read.delim(file, header = F, sep = "\t")
  tbl = tbl[-which(tbl$V1 == "--------------------------------------"),]
  tbl = tbl[-grep("Number of cleavage sites", x = tbl)]
  tbl = tbl[-grep("NetChop 3.0", x = tbl)]
  tbl = tbl[-grep("pos", x = tbl)]
  
  
  tbl = str_split_fixed(tbl, pattern = "[:space:]+", Inf) %>%
    as.data.frame()
  tbl = tbl[,-1]
  
  names(tbl) = c("pos", "AA", "C", "score", "Ident")
  tbl$score = as.numeric(tbl$score)
  
  return(tbl)
}

netchop = readInNetChop(file)


# ----- validate NetChop prediction: correlations -----
zscale = function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}

subIDs = DATA$substrateIDs %>% unique()
X = DATA$X
X = X[, colnames(X) %in% params]

jointPosterior = parameters[sample(nrow(parameters), 10e03), ]

pdf("results/Bayesian_ProteaSMM/PREDICTION/NetChop-vs-BayesianProteaSMM/_proteins_correlations.pdf", height = 5, width = 10)
par(mfrow = c(1,2))

NetChopROCs = list()
ProteaSMMROCs = list()
for(i in 1:length(subIDs)) {
  print(subIDs[i])
  keep = which(DATA$substrateIDs == subIDs[i])
  
  cntNetChop = netchop$score[netchop$Ident == subIDs[i]]
  cntTrueNetChop = DATA$t[keep, ]
  
  cntPred = apply(jointPosterior,1,function(p){
    tsim = X[keep,]%*%log(p)
    return(tsim)
  })
  cntTrue = log(cntTrueNetChop/100+pseudo)
  
  # NetChop vs. true
  tmean = apply(cntTrueNetChop,1,mean,na.rm = T)
  tsd = apply(cntTrueNetChop,1,sd,na.rm = T)
  pcc = cor(tmean, cntNetChop)
  plot(x = tmean, y = cntNetChop, pch = 16,
       xlab = "SCS-P1 (%)", ylab = "NetChop score",
       main = paste0(subIDs[i], ", PCC = ", round(pcc, 4)),
       sub = "NetChop")
  arrows(x0 = tmean+tsd,
         x1 = tmean-tsd,
         y0 = cntNetChop, y1 = cntNetChop,
         code = 3, angle = 90, length = 0.03, lwd = .5) %>% suppressWarnings()
  
  NetChopROCs[[i]] = getROCcurve(ttrue_mean = zscale(tmean), tsims_mean = zscale(cntNetChop), substrate = subIDs[i], threshPerc = 0.5, retValsOnly = T)[1,c("pr_auc","roc_auc")]
  
  # Bayesian ProteaSMM vs. true
  tmean = apply(cntTrue,1,mean,na.rm = T)
  tsd = apply(cntTrue,1,sd,na.rm = T)
  tpredmean = apply(cntPred,1,mean,na.rm = T)
  tpredsd = apply(cntPred,1,sd,na.rm = T)
  
  pcc = cor(tmean, tpredmean)
  plot(x = tmean, y = tpredmean, pch = 16,
       xlab = "log-transformed SCS-P1 (%)", ylab = "Bayesian ProteaSMM prediction",
       main = paste0(subIDs[i], ", PCC = ", round(pcc, 4)),
       sub = "Bayesian ProteaSMM")
  arrows(x0 = tmean+tsd,
         x1 = tmean-tsd,
         y0 = tpredmean, y1 = tpredmean,
         code = 3, angle = 90, length = 0.03, lwd = .5) %>% suppressWarnings()
  arrows(x0 = tmean, x1 = tmean,
         y0 = tpredmean+tpredsd, y1 = tpredmean-tpredsd,
         code = 3, angle = 90, length = 0.03, lwd = .5) %>% suppressWarnings()
  
  ProteaSMMROCs[[i]] = getROCcurve(ttrue_mean = zscale(tmean), tsims_mean = zscale(tpredmean), substrate = subIDs[i], threshPerc = 0.5, retValsOnly = T)[1,c("pr_auc","roc_auc")]
  
}
dev.off()


# ----- validate NetChop prediction: binary performance -----

names(ProteaSMMROCs) = subIDs
names(NetChopROCs) = subIDs

ProteaSMMROCs = ProteaSMMROCs %>% plyr::ldply()
NetChopROCs = NetChopROCs %>% plyr::ldply()



### OUTPUT ###