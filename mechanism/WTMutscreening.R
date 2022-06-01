### karat projetc - PCPS mechanism ###
# description:  check in WT/Mut dataset for identical PSPs with mutation @ P1
# input:        invitroSPI results of WT/Mut dataset
#               peptides with mutations @ P1 (or P1')
# author:       HPR

library(dplyr)
library(stringr)
source("src/_extract-aa.R")
source("src/invitroSPI_utils.R")


### INPUT ###
WTMut = read.csv("../DATA/WT-Mut/invitroSPI/OUTPUT/WT_Mut_0_4/ProteasomeDB.csv",
                 stringsAsFactors = F)
WTMut_late = read.csv("../DATA/WT-Mut/invitroSPI/OUTPUT/WT_Mut_0_20/ProteasomeDB.csv",
                      stringsAsFactors = F)

WTMutsubNames = read.csv("../DATA/WT-Mut/qiSPI/INPUT/sample_list.csv", stringsAsFactors = F) %>%
  select(protein_name, substrateID) %>%
  unique()


### MAIN PART ###
suppressWarnings(dir.create("results/WTMut/"))

# ----- preprocessing -----
WTMut = rbind(WTMut, WTMut_late)
WTMut = left_join(rbind(WTMut,WTMut_late), WTMutsubNames)

Qual = WTMut %>%
  ILredundancy() %>%
  uniquePeptides() %>%
  disentangleMultimappers.Type() %>%
  tidyr::separate_rows(positions, sep = ";")
  
  
Qual$spliceType[is.na(Qual$spliceType)] = "PCP"

# get splice-reactants
extractSRs = function(DB) {
  
  pos = str_split_fixed(DB$positions, "_", Inf)[,c(1:4)]
  pos = apply(pos,2,as.numeric)
  
  DB$sr1 = substr(DB$substrateSeq, pos[,1], pos[,2])
  DB$sr2 = substr(DB$substrateSeq, pos[,3], pos[,4])
  
  DB$posP1 = pos[,2]
  DB$posP1_ = pos[,3]
  DB$posP1_[DB$spliceType == "PCP"] = pos[DB$spliceType == "PCP",2]+1
  return(DB)
}

# get amino acids
DB = extractSRs(Qual)
DB = left_join(DB, extract_aminoacids(DB, onlyValidSeq = T))


# ----- extract mutation position -----
DB$protein_name %>% unique()
DB$origSubs = str_extract_all(DB$protein_name, "^[:alnum:]+(?=_)", simplify = T)
DB$origSubs %>% unique()

MutPos = DB %>%
  group_by(origSubs, protein_name, substrateSeq) %>%
  summarise() %>% ungroup() %>%
  # filter(origSubs != "KRAS") %>%  # !!!
  mutate(WT = ifelse(grepl("WT", protein_name), T, F)) %>%
  as.data.frame()

MutPos$residue = NA
orgs = unique(MutPos$origSubs)
for (o in orgs) {
  wt = which(MutPos$origSubs == o & MutPos$WT)
  cntMut = which(MutPos$origSubs == o & !MutPos$WT)
  
  x = sapply(cntMut, function(mut) {
    mapply(function(wt,mut) which(bitwXor(utf8ToInt(wt), utf8ToInt(mut)) > 0), MutPos$substrateSeq[wt], MutPos$substrateSeq[mut]) %>%
      as.numeric()
  })
  print(x)
  MutPos$residue[cntMut] = x
  MutPos$residue[rep(wt, each = length(unique(x)))] = unique(x)
}

# add the KRAS position
MutPos[nrow(MutPos)+1,] = c("KRAS",
                            "KRAS_WT",
                            "TEYKLVVVGAGDVGKSALTLQLLQNHFVDEYDPT",
                            TRUE,
                            11)

DBj = left_join(DB, MutPos)

X = DBj %>%
  filter(posP1 == residue | posP1_ == residue) %>%
  mutate(identical = ifelse(posP1 == residue, "P1", "P1_")) %>%
  group_by(origSubs, positions) %>%
  mutate(same = length(unique(WT))) %>%
  filter(same >= 2) %>%
  select(substrateID, protein_name, origSubs, spliceType, pepSeq, positions, P1, P1_, identical) %>%
  ungroup() %>%
  arrange(positions)

table(X$identical)
table(X$spliceType) /2

### OUTPUT ###
write.csv(X, "results/WTMut/WT-Mut-paris.csv", row.names = F)
