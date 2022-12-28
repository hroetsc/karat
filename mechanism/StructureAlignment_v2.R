### karat project - PCPS mechanism ###
# description:  align proteasome subunits of yeast and human
# input:        b1, b2 and b5 subunits of yeast and human proteasomes
# output:       alignment at crucial positions: 1,17,33,47,130,167,170
# author:       HPR


library(seqinr)
library(dplyr)
library(stringr)
library(Biostrings)
library(msa)


### INPUT ###
allBeta = readAAStringSet("data/20S_sequences/yeast+human_catalytic.fasta",format = "fasta", use.names = T)


### MAIN PART ###
pnames = names(allBeta)

MSA = msaClustalOmega(allBeta)

print(MSA)
conMat = consensusMatrix(MSA)
consensus = msaConsensusSequence(MSA)


### OUTPUT ###
