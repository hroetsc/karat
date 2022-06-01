### karat projetc - PCPS mechanism ###
# description:  specificity vs lengths - propensity scale approach
# input:        invitroSPI: Roetschke SciData, EGFR ery, WT substrates --> spliced peptides only
# output:       most dissimilar amino acid descriptors between P1 and P1'
# author:       HPR

library(dplyr)
library(stringr)
library(seqinr)
library(protr)
library(Peptides)
library(transport)
library(ggplot2)

theme_set(theme_classic())

### INPUT ###
load("results/SRspecificity/DBproc.RData")
load("results/SRspecificity/randomDBproc.RData")


### MAIN PART ###
data("AADescAll")
force(AADescAll)
AADescAll$AA = rownames(AADescAll)

# get the most dissimilar amino acid descriptors between SR1 and SR2
getaapropensities = function(residues, prop, norm = T) {
  
  x1 = data.frame(AA = residues) %>%
    left_join(prop) %>% 
    select(-AA) %>% suppressMessages()
  
  if (norm) {
    x1 = apply(x1,2,function(c){
      return((c-min(c))/(max(c)-min(c)))
    })
  }
  
  return(x1)
}

P1prop = getaapropensities(residues = DBproc$P1, prop = AADescAll)
P1_prop = getaapropensities(residues = DBproc$P1_, prop = AADescAll)

AADescAll = AADescAll %>% select(-AA)
# get pairwise similarities between amino acid indices
Sims = list()
for (a in 1:(ncol(AADescAll))) {
  Sims[[a]] = wasserstein1d(P1prop[,a], P1_prop[,a])
}

names(Sims) = colnames(AADescAll)
Sims = unlist(Sims)
summary(Sims)

# most dissimilar amino acid indices between SR1 and SR2
plot(density(Sims))
extracted = Sims[Sims >= quantile(Sims, 0.98)]
extracted

descriptions = sapply(names(extracted), function(x){
  AADescAll[names(AADescAll) %in% x][[1]][["D"]]
})

extractedInd = data.frame(AAindex = names(extracted),
                          desc = descriptions,
                          sr1sr2_sim = extracted)


xtr8ed = AADescAll[, colnames(AADescAll) %in% names(extracted)]
save(xtr8ed, file = "results/SRspecificity/xtr8ed_AADecsAll.RData")
