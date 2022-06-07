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
finalK = list.files("invitroSPI+aSPIre/results/", pattern = "finalKinetics.csv",
                    recursive = T, full.names = T)
ProteasomeDB_SciData = read.csv("../../invitroSPI/revision/results/ProteasomeDB/ProteasomeDB.csv", stringsAsFactors = F)
ProteasomeDB_BScthesisQuant = read.csv("invitroSPI/OUTPUT/BScthesisQuant/ProteasomeDB.csv", stringsAsFactors = F)


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

KU = Kinetics %>%
  distinct(substrateID, pepSeq, .keep_all = T)
table(KU$substrateID,KU$spliceType)
table(KU$spliceType)
table(KU$productType)

KU$substrateID %>% unique() %>% length()

# ----- parse qualitative data set -----
nm = intersect(names(ProteasomeDB_SciData), names(ProteasomeDB_BScthesisQuant))
ProteasomeDB = rbind(ProteasomeDB_SciData[,nm],
                     ProteasomeDB_BScthesisQuant[,nm]) %>%
  filterPepLength() %>%
  disentangleMultimappers.Type() %>%
  remSynthErrors() %>%
  filterEarlyTimepoints() %>%
  filter20Sstandard()
ProteasomeDB$spliceType[ProteasomeDB$spliceType %in% c(NA,"")] = "PCP"

PU = ProteasomeDB %>%
  distinct(substrateID, pepSeq, .keep_all = T)
table(PU$substrateID,PU$spliceType)
table(PU$spliceType)
table(PU$productType)

PU$substrateID %>% unique() %>% length()

# ----- parse proteins -----


### OUTPUT ###
save(Kinetics, file = "data/aSPIre.RData")
write.csv(Kinetics, file = "data/aSPIre.csv", row.names = F)

save(ProteasomeDB, file = "data/invitroSPI.RData")
write.csv(ProteasomeDB, file = "data/invitroSPI.csv", row.names = F)
