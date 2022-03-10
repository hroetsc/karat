### INHIBITOR KINETICS ###
# description:  select canditates based on SCS/PSP-P1
# input:        finalKinetics, selected peptides, manually selected residues
# output:       candidate precursors
# author:       HR

library(dplyr)
library(stringr)
source("../brainstorming/src/invitroSPI_utils.R")

### INPUT ###
finalKinetics = read.csv("qiSPI/OUTPUT/TSN5_0+4/finalKinetics.csv", stringsAsFactors = F)
selectedPeps = read.csv("results/b5/selectedPeps_intensities.csv", stringsAsFactors = F)
candidates = read.csv("results/b5/candidates.csv", stringsAsFactors = F)

### MAIN PART ###
# get sP1s
# find out which peptides actually have this residue as (s)P1
# for each candidate residue, find a PCP precursor

# ----- preprocessing -----
DB = finalKinetics %>%
  ILredundancy() %>%
  disentangleMultimappers.AA() %>%
  removeMultimappers.AA()

names(DB)
DB = DB[, c("pepSeq","productType","spliceType","positions","substrateID","substrateSeq")] %>%
  unique()

DB_selected = finalKinetics %>%
  filter(pepSeq %in% selectedPeps$pepSeq) %>%
  ILredundancy() %>%
  disentangleMultimappers.AA() %>%
  removeMultimappers.AA()
DB_selected = DB_selected[, c("pepSeq","productType","spliceType","positions","substrateID","substrateSeq")] %>%
  unique()


# ----- extract P1s -----
extractP1 = function(DB) {
  pos = str_split_fixed(DB$positions, coll("_"), Inf) %>% as.data.frame()
  
  DB$Nterm = NA
  DB$P1 = NA
  DB$P1 = as.numeric(pos$V2)
  DB$Nterm = as.numeric(pos$V1)
  
  DB$P1_aa = NA
  DB$P1_aa = substr(DB$substrateSeq, start = pos$V2, stop = pos$V2)
  
  return(DB)
}

DB = extractP1(DB)
DB_selected = extractP1(DB_selected)

# ----- extract peptides that have candidate residue as P1 in selected peptides -----
candidates$pos_PSPs = NA
candidates$pos_PCPs = NA

for (i in 1:nrow(candidates)) {
  
  j_psp = which(DB_selected$P1 == candidates$position[i] & DB_selected$productType == "PSP")
  j_pcp = which(DB_selected$P1 == candidates$position[i] & DB_selected$productType == "PCP")
  
  candidates$pos_PSPs[i] = paste(DB_selected$positions[j_psp], collapse = ";")
  candidates$pos_PCPs[i] = paste(DB_selected$positions[j_pcp], collapse = ";")
}



# ----- extract precursors -----
# in entire DB and in selected DB
# check which peptides have the same N terminus as a detected peptide

selectPrecursors = function(DB, candidates) {
  
  candidates$precursors_PSPs = NA
  candidates$precursors_PCPs = NA
  
  for (i in 1:nrow(candidates)) {
    
    # PSPs
    cdts = paste(candidates$pos_PSPs[i], candidates$pos_PCPs[i], collapse = ";") %>%
      strsplit(";") %>%
      unlist()
    cdts = str_remove_all(cdts, coll(" "))
    
    # select precursors that contain the candidate residue
    # their C-term should include any other candidate residue
    j_psp = which((candidates$position[i] > DB$Nterm) &
                    (candidates$position[i] < DB$P1) &
                    (!DB$P1 %in% candidates$position[-i]) &
                    (DB$productType == "PSP"))
    j_pcp = which((candidates$position[i] > DB$Nterm) &
                    (candidates$position[i] < DB$P1) &
                    (!DB$P1 %in% candidates$position[-i]) &
                    (DB$productType == "PCP"))
    
    # select precursors with the same N-term as a candidate peptide
    if (cdts != "") {
      cdts_N = str_split_fixed(cdts, coll("_"), Inf)[,1] %>%
        as.numeric()
      
      # PSP
      if (length(j_psp) > 0) {
        jj1 = which((str_split_fixed(DB$positions[j_psp], coll("_"), Inf)[,1] %>% as.numeric()) == cdts_N)
        if (length(jj1)) {
          j_psp = j_psp[jj1]
        }
      }
      
      # PCP
      if (length(j_pcp) > 0) {
        jj2 = which((str_split_fixed(DB$positions[j_pcp], coll("_"), Inf)[,1] %>% as.numeric()) == cdts_N)
        if (length(jj2)) {
          j_pcp = j_pcp[jj2]
        }
      }
      
    }
    
    candidates$precursors_PSPs[i] = paste(DB$positions[j_psp], collapse = ";")
    candidates$precursors_PCPs[i] = paste(DB$positions[j_pcp], collapse = ";")
  }
  
  return(candidates)
}

candidates2 = selectPrecursors(DB_selected, candidates)

### OUTPUT ###
write.csv(candidates2, "results/b5/candidates_enriched.csv", row.names = F)
