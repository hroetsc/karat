
library(dplyr)
library(tidyr)
library(stringr)

setwd("Documents/QSB/activeSite/brainstorming/")
pdb = read.csv("../../invitroSPI/revision/results/ProteasomeDB/ProteasomeDB.csv",
               stringsAsFactors = F)


# catalytic triade
acid_residue = c("E", "D")
base_residue = c("H")
nucleophile = c("S","C","T")
catalytic_triade = crossing(acid_residue, base_residue, nucleophile)
catalytic_triade = apply(catalytic_triade, 1,paste, collapse="")

subs = pdb$substrateSeq %>% unique()
subs_triade = sapply(catalytic_triade, str_detect, string=subs)

peps = pdb$pepSeq %>% unique()
peps_triade = sapply(catalytic_triade, str_detect, string=peps)

k = apply(peps_triade, 1, any) %>% which()
peps[k]

pdb[pdb$pepSeq %in% peps[k], ] %>%
  unique() %>%
  View()

# catalytic duo
catalytic_duo = c("EH", "DH")
subs_duo = sapply(catalytic_duo, str_detect, string=subs)
peps_duo = sapply(catalytic_duo, str_detect, string=peps)

k = apply(peps_duo, 1, any) %>% which()
peps[k]

x = pdb[pdb$pepSeq %in% peps[k], ] %>%
  distinct(substrateID, productType, spliceType, .keep_all = T)

length(which(x$productType == "PSP")) / nrow(x)
length(which(x$productType == "PCP")) / nrow(x)

pdb_u = pdb %>%
  distinct(substrateID, productType, spliceType, .keep_all = T)

length(which(pdb_u$productType == "PSP")) / nrow(pdb_u)
length(which(pdb_u$productType == "PCP")) / nrow(pdb_u)


kk = apply(subs_duo, 1, any) %>% which()
xx = pdb[pdb$substrateSeq %in% subs[kk], ] %>%
  distinct(substrateID, productType, spliceType, .keep_all = T)

length(which(xx$productType == "PSP")) / nrow(xx)
length(which(xx$productType == "PCP")) / nrow(xx)


fracs_pcps = rep(NA, length(subs))
fracs_psps = rep(NA, length(subs))
for (s in 1:length(subs)) {
  
  cnt = pdb_u[pdb_u$substrateSeq == subs[s], ]
  
  pcp = cnt[cnt$productType == "PCP", ]
  psp = cnt[cnt$productType == "PSP", ]
  
  fracs_pcps[s] = (apply(sapply(catalytic_duo, str_detect, string=pcp$pepSeq) %>%
          as.data.frame(), 1, any) %>%
    which() %>%
    length()) / nrow(pcp)
  
  fracs_psps[s] = (apply(sapply(catalytic_duo, str_detect, string=psp$pepSeq) %>%
                           as.data.frame(), 1, any) %>%
                     which() %>%
                     length()) / nrow(psp)
  
}

summary(fracs_psps)
summary(fracs_pcps)

