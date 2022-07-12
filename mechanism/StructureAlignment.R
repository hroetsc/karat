### karat projetc - PCPS mechanism ###
# description:  align proteasome subunits of different organisms
# input:        b1, b2 and b5 subunits of different proteasomes
# output:       alignment at crucial positions: 1,17,33,47,130,167,170
# author:       HPR

library(dplyr)
library(stringr)
library(seqinr)
library(Biostrings)
library(msa)
library(ggmsa)

### INPUT ###
allBeta = readAAStringSet("data/structure/Uniprot_allProteasomes.fasta",format = "fasta", use.names = T)

### MAIN PART ###
# ----- disentangle beta subunits from each other -----
pnames = names(allBeta)
allBetaF = allBeta[grepl("Proteasome", pnames)]
pnames = names(allBetaF)

# beta 5 only
b5names = pnames[grepl("5 OS=", pnames)]
allBeta5 = allBetaF[b5names]

MSA = msaClustalOmega(allBeta5)

print(MSA)
conMat = consensusMatrix(MSA)

protseq = msaConvert(MSA, type = "bio3d::fasta")

consensus = msaConsensusSequence(MSA)
consensus = substr(consensus,89,nchar(consensus))

consensus[c(17,33,47,130,167,170)]


ggmsa(allBeta5, seq_name = T, start = 75, end = 100) + geom_seqlogo() + geom_msaBar()


### OUTPUT ###


