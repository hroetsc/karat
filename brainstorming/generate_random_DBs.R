### PCPS mechanism brainstorming ###
# description:  generate random databases
# input:        ProteasomeDB, _allPSP fasta files
# output:       random databases
# author:       HR

library(dplyr)
library(stringr)
library(Biostrings)

source("src/invitroSPI_utils.R")
source("src/loadFunctions.r")

rndSize = 10


### INPUT ###
DBpath = "~/Documents/QSB/invitroSPI/revision/"
ProteasomeDB = read.csv(paste0(DBpath, "results/ProteasomeDB/ProteasomeDB.csv"),
                        stringsAsFactors = F)

FASTApath = "/Volumes/DATA16040/USERS/Juliane/ProteasomeDB/PCPS_DB/MasterDB_all_2022/FASTA/"
PSPdb = list.files(FASTApath, pattern = ".fasta",
                   recursive = T, full.names = T)

length(PSPdb)

### MAIN PART ###
# ----- 1) DB preprocessing -----
ProteasomeDB = ProteasomeDB %>%
  ILredundancy() %>%
  filterPepLength() %>%
  filter20Sstandard() %>%
  # filterEarlyTimepoints() %>% # !!!!!
  uniquePeptides.perTP() %>%
  disentangleMultimappers.Type() %>%
  disentangleMultimappers.SRlen() %>%
  disentangleMultimappers.IVSeqLen() %>%
  disentangleMultimappers.AA() %>%
  DBcosmetics()


# ----- 2) sample random DBs ----
tps = ProteasomeDB$digestTime %>% unique()
tps = c(2,4)  # !!!!!

S = ProteasomeDB$substrateID %>% unique()
protTypes = c("PCP","cis","revCis","trans")

randomDB.allSubs = list()
for (s in 1:length(S)) {
  
  print("----------")
  print(S[s])
  print("----------")
  
  if (length(PSPdb) == 1) {
    allPSP = readAAStringSet(filepath = PSPdb, format = "fasta")
  } else {
    allPSP = readAAStringSet(filepath = PSPdb[str_detect(PSPdb, paste0(S[s],"_allPSP.fasta"))][1], format = "fasta")
  }
  
  
  pos = str_split(allPSP@ranges@NAMES, "_", simplify = T)
  pos = pos[, (2:5)]
  pos = apply(pos, 2, as.numeric)
  
  cntpcp = which(is.na(pos[,4]))
  
  intv = pos[,3]-pos[,2]
  cntcis = which(intv > 0 & !is.na(intv) & pos[,3] > pos[,2])
  cntrevcis = which(intv <= 0 & !is.na(intv) & pos[,4] <= pos[,1])
  cnttrans = which(intv <= 0 & !is.na(intv) & pos[,4] > pos[,1])
  # length(cntpcp) + length(cntcis) + length(cntrevcis) + length(cnttrans) == nrow(pos) 
  
  allProds = as.character(allPSP)
  
  randomDB.allTPs = list()
  for (t in 1:length(tps)) {
    
    print("----------------------------------")
    print("generate random DB for time point:")
    print(tps[t])
    print("----------------------------------")
    
    cntDB = ProteasomeDB[ProteasomeDB$substrateID == S[s] & ProteasomeDB$digestTime == tps[t], ]
    
    if (nrow(cntDB) > 0) {
      
      # get fraction of products
      cisFrac = length(which(cntDB$spliceType == "cis")) * rndSize
      revCisFrac = length(which(cntDB$spliceType == "revCis")) * rndSize
      transFrac = length(which(cntDB$spliceType == "trans")) * rndSize
      PCPFrac = length(which(cntDB$spliceType == "PCP")) * rndSize
      
      # sample number of products
      cis = sample(cntcis, size = cisFrac)
      revcis = sample(cntrevcis, size = revCisFrac)
      trans = sample(cnttrans, size = transFrac)
      pcp = sample(cntpcp, size = PCPFrac, replace = T)
      
      # sample randomDB
      cntRnd = data.frame(matrix(ncol = ncol(ProteasomeDB), nrow = length(cis)+length(revcis)+length(trans)+length(pcp)))
      colnames(cntRnd) = colnames(ProteasomeDB)
      
      cntRnd$substrateID = cntDB$substrateID[1]
      cntRnd$digestTime = tps[t]
      cntRnd$substrateSeq = cntDB$substrateSeq[1]
      
      cntRnd$productType = c(rep("PSP", length(cis)+length(revcis)+length(trans)),
                             rep("PCP", length(pcp)))
      cntRnd$spliceType = c(rep("cis", length(cis)),
                            rep("revCis", length(revcis)),
                            rep("trans", length(trans)),
                            rep("PCP", length(pcp)))
      
      cntRnd$positions = apply(pos[c(cis, revcis, trans, pcp),], 1, paste, collapse = "_")
      cntRnd$pepSeq = c(allProds[c(cis, revcis, trans, pcp)])
      
      randomDB.allTPs[[t]] = cntRnd
      
    }
  }
  
  randomDB.allSubs[[s]] = plyr::ldply(randomDB.allTPs)
}

d = plyr::ldply(randomDB.allSubs)


# ----- 3) get multi-mappers ----
d = mapping(d)


# ----- 4) sanity check -----



### OUTPUT ####
save(d, file = "data/randomDB_earlyTPs.RData")
