### karat projetc - PCPS mechanism ###
# description:  compare number of detected peptides with SCS/PSP-P1
# input:        aSPIre: Roetschke SciData, EGFR ery, WT substrates (quantitative DB)
# output:       correlation between peptide counts and SCS/PSP-P1
# author:       HPR

library(dplyr)
library(stringr)
library(pheatmap)
source("src/invitroSPI_utils.R")
source("src/_extract-aa.R")
source("src/SCS+PSP-P1.R")

suppressWarnings(dir.create("data/ProteaSMM/"))

AAchar_here = c("P","G","C","M","A","V","I","L","F","Y","W","H","R","K","D","E","N","Q","S","T","X")
AAchar_here_sorted = sort(AAchar_here)

### INPUT ###
load("data/aSPIre.RData")

### MAIN PART ###
# ----- preprocessing -----
Quant = Kinetics %>%
  # ILredundancy() %>%
  disentangleMultimappers.Type() %>%
  tidyr::separate_rows(digestTimes, intensities, sep=";") %>%
  rename(digestTime = digestTimes,
         intensity = intensities) %>%
  mutate(intensity = as.numeric(intensity),
         digestTime = as.numeric(digestTime)) %>%
  filter(digestTime == 4) %>%
  resolve_multimapper() %>%
  tidyr::separate_rows(positions, sep = ";")

# ----- get cleavage/splicing strength -----
# get SCS and PSP for P1 or P1' for each bio rep
outP1 = SCSandPSP_allSubs(DB = Quant, target = "P1", meanOverBio = F, Zscale = F)

yPredDF = plyr::ldply(outP1) %>%
  as.data.frame()
yPredDF$target = "P1"
names(yPredDF)[1] = "substrateID"

scs = str_split_fixed(yPredDF$scs_mean, ";", Inf) %>% as.matrix()
scs = apply(scs,2,as.numeric)

psp = str_split_fixed(yPredDF$psp_mean, ";", Inf) %>% as.matrix()
psp = apply(psp,2,as.numeric)


# ----- counts vs strength -----
png("results/SCS+PSP/counts-vs-strength.png", height = 4, width = 8, units = "in", res = 600)
par(mfrow = c(1,2))

plot(x = yPredDF$scs_n,
     y = rowMeans(scs),
     pch = 16, cex = .7, xlab = "number of mapping peptides", ylab = "SCS-P1 (%)")

plot(x = yPredDF$psp_n,
     y = rowMeans(psp),
     pch = 16, cex = .7, xlab = "number of mapping peptides", ylab = "PSP-P1 (%)")

dev.off()

### OUTPUT ###
