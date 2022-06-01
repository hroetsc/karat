### karat projetc - PCPS mechanism ###
# description:  generate random quantitative random databases
# input:        qualitative random DB, identified sequences + MS1 intensities
# output:       quantitative random database
# author:       HPR

library(dplyr)
library(stringr)
source("src/invitroSPI_utils.R")

### INPUT ###
load("data/aSPIre.RData")
load("data/randomDB_aSPIre.RData")


### MAIN PART ###
# ----- preprocessing -----
Kinetics = Kinetics %>%
  tidyr::separate_rows(digestTimes, intensities) %>%
  mutate(digestTime = as.numeric(digestTimes),
         intensity = as.numeric(intensities)) %>%
  select(-digestTimes, -intensities) %>%
  filter(digestTime == 4)


# ----- sample MS1 intensities -----
rndDB$spliceType[rndDB$productType == "PCP"] = "PCP"
randomQuant = rndDB %>%
  disentangleMultimappers.Type()
randomQuant$intensity = NA

subIDs = Kinetics$substrateID %>% unique()
Types = Kinetics$spliceType %>% unique()

for (i in 1:length(subIDs)) {
  for (j in 1:length(Types)) {
    
    y = Kinetics$intensity[Kinetics$substrateID == subIDs[i] & Kinetics$spliceType == Types[j]]
    if (is.null(y) | length(y) == 0) {
      y = 0
    }
    k = which(randomQuant$substrateID == subIDs[i] & randomQuant$spliceType == Types[j])
    randomQuant$intensity[k] = sample(y, size = length(k), replace = T)
    
  }
}


### OUTPUT ###
save(rndDB, file = "data/randomDB_aSPIre.RData")
save(randomQuant, file = "data/randomDB_Quant_aSPIre.RData")

