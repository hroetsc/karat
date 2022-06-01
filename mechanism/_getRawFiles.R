### karat projetc - PCPS mechanism ###
# description:  fetch all raw files for quantification in Skyline
# input:        EGFR and WT-Mut invitroSPI assignments
# output:       list of .raw files
# author:       HPR

library(dplyr)

### INPUT ###
WT_samplelist = read.csv("../DATA/WT-Mut/qiSPI/INPUT/sample_list.csv", stringsAsFactors = F)
EGFR_samplelist = read.csv("../DATA/EGFR/INPUT/sample_list.csv", stringsAsFactors = F)
EGFR_metainformation = read.csv("../DATA/EGFR/INPUT/metainformation_EGFR.csv", stringsAsFactors = F)

load("../DATA/EGFR/OUTPUT/EGFR/tmp/allPSMs.RData")
EGFR_allPSMs = allPSMs

load("../DATA/WT-Mut/invitroSPI/OUTPUT/WT_Mut_0_4/tmp/allPSMs.RData")
WT1_allPSMs = allPSMs


### MAIN PART ###

# ----- get standard protea from EGFR -----
# no mutation and only standard proteasome
# 0-4 hours measurements on QE

EGFR_relevant = EGFR_metainformation[EGFR_metainformation$protIsotype == "20 standard (ery)" & 
                                       EGFR_metainformation$substrateID %in% c("MM835","MM836") &
                                       EGFR_metainformation$instrument == "QEx1", ]
EGFR_samplelist = EGFR_samplelist[EGFR_samplelist$filename %in% EGFR_relevant$filename, ]

EGFR_samplelist$substrateName = NA
EGFR_samplelist$substrateName[EGFR_samplelist$substrateID == "MM835"] = "EGFR_1-45"
EGFR_samplelist$substrateName[EGFR_samplelist$substrateID == "MM836"] = "EGFR_pepVIII"

EGFR_allPSMs = EGFR_allPSMs[EGFR_allPSMs$filename %in% EGFR_samplelist$filename, ]

# ----- get WT from WT-Mut -----
# only WT measurements
# 0-4 hours measurements on QE


WT_relevant = WT_samplelist
# WT_relevant = WT_samplelist[WT_samplelist$project_name %in% c("rest","WT_Mut_0_4") &
#                               WT_samplelist$substrateID %in% c("MM582","MM336","MM539","MM537","MM499","MM136"), ]
WT_relevant$substrateName = NA

# match MM IDs to common names
WT_relevant$substrateName[WT_relevant$substrateID == "MM582"] = "IDH1"
WT_relevant$substrateName[WT_relevant$substrateID == "MM336"] = "NRAS"
WT_relevant$substrateName[WT_relevant$substrateID == "MM539"] = "MPL"
WT_relevant$substrateName[WT_relevant$substrateID == "MM537"] = "JAK2"
WT_relevant$substrateName[WT_relevant$substrateID == "MM499"] = "BRAF"
WT_relevant$substrateName[WT_relevant$substrateID == "MM136"] = "KRAS"


# ----- combine both lists -----

masterlist = rbind(EGFR_samplelist, WT_relevant)
masterlist$project_name = "BScthesisQuant"
masterlist$metainformation = NA

# ----- fetch search results -----

fs = sapply(masterlist$filename,function(x){
  list.files(path = "../DATA/", pattern = x, full.names = T, recursive = T)
}) %>%
  unlist()

masterlist$filename[!masterlist$filename %in% basename(fs) %>% unique()]

able2copy = file.copy(from = fs, to = "invitroSPI/INPUT/search_results/", overwrite = T)
all(able2copy)


### OUTPUT ###
write.csv(masterlist, "invitroSPI/INPUT/sample_list.csv", row.names = F)



