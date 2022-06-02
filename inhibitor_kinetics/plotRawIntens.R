### INHIBITOR KINETICS ###
# description:  plot raw intensities of selected peptides for 0,2,4 hrs
# input:        selectedPeps, qiSPI output for 0+4 and 0+2 hrs
# output:       plots for selected peptides
# author:       HR


library(dplyr)
library(stringr)
source("../../brainstorming/src/invitroSPI_utils.R")
source("src/aSPIre_plotting.R")


### INPUT ###
# selected peptides
selectedPeps = read.csv("results/b5/selectedPeps_intensities.csv", stringsAsFactors = F)
finalKinetics = read.csv("aSPIre_manual/results/TSN5inhibitor/finalKinetics.csv", stringsAsFactors = F)


### MAIN PART ###
# ----- preprocessing -----
DATA = finalKinetics
DATA$spliceType[is.na(DATA$spliceType)] = "PCP"
DATA = na.omit(DATA)
DATA = DATA %>%
  disentangleMultimappers.Type()

# add full replicate info
DATA$condition = paste0(str_replace_all(DATA$biological_replicate, "_(?=[:digit:]$)", "_bio"),
                        "_tech", DATA$technical_replicate)
DATA = DATA[order(DATA$pepSeq), ]

Qtable = DATA %>%
  tidyr::separate_rows(digestTimes, intensities, sep = ";") %>%
  mutate(digestTime = as.numeric(digestTimes),
         intensity = as.numeric(intensities))


# ----- plotting -----
pepSeqs = selectedPeps$pepSeq

plotKinetics(Qtable, outfile = "results/rawIntensities.pdf", meanTech = T, earlyOnly = T, sortByInt = T)
plotKinetics(Qtable[Qtable$pepSeq %in% selectedPeps$pepSeq, ], outfile = "results/b5/rawIntensities_selectedPeps.pdf", meanTech = T, earlyOnly = T, sortByInt = T)


# pdf("results/b5/rawIntensities-for-candidates.pdf", height = 16, width = 16)
# par(mfrow = c(4,4))

# png("~/Documents/Studium/Fachvertiefung+BA/Praktikumsbericht/plots/rawIntensities-for-candidates.png",
#     height = 26*2, width = 20*2, units = "cm", res = 600)




