### karat projetc - PCPS mechanism ###
# description:  fetch invitroSPI and aSPIre results for all substrates
#               basic statistics
# input:        aSPIre: EGFR ery, WT substrates
#               invitroSPI: Specht et al, Paes et al, gp100, EGFR, WT
# output:       joined qualitative and quantitative data sets
# author:       HPR


library(dplyr)
library(stringr)
source("src/invitroSPI_utils.R")


### INPUT ###
finalK = list.files("aSPIre_Mut/results/", pattern = "finalKinetics.csv",
                    recursive = T, full.names = T)



### MAIN PART ###
# ----- parse quantitative data set -----
Kinetics = lapply(finalK, function(x){
  read.csv(x,stringsAsFactors = F)
}) %>%
  plyr::ldply() %>%
  as.data.frame()

Kinetics$spliceType[is.na(Kinetics$spliceType)] = "PCP"
Kinetics = na.omit(Kinetics)
Kinetics = Kinetics %>%
  disentangleMultimappers.Type()

# remove peptides that have 0 intensity @ 4 hours
RemPeps = Kinetics %>%
  tidyr::separate_rows(digestTimes, intensities, sep=";") %>%
  rename(digestTime = digestTimes,
         intensity = intensities) %>%
  mutate(intensity = as.numeric(intensity),
         digestTime = as.numeric(digestTime)) %>%
  filter(digestTime == 4 & intensity == 0)
RemPeps = RemPeps$pepSeq %>% unique()

Kinetics = Kinetics[-which(Kinetics$pepSeq %in% RemPeps), ]

KU = Kinetics %>%
  distinct(substrateID, pepSeq, .keep_all = T)
table(KU$substrateID,KU$spliceType)
table(KU$spliceType)
table(KU$productType)

KU$substrateID %>% unique() %>% length()



### OUTPUT ###
save(Kinetics, file = "data/MutPairs_aSPIre.RData")
write.csv(Kinetics, file = "data/MutPairs_aSPIre.csv", row.names = F)
