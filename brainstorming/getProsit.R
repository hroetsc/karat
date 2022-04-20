# transform invitroSPI output to get spectral angles

library(dplyr)
library(stringr)
source("../../proteinsPCPS/new/src/invitroSPI_utils.R")

sample_list = read.csv("invitroSPI/INPUT/sample_list.csv", stringsAsFactors = F)
DB = read.csv("invitroSPI/OUTPUT/TSN5/ProteasomeDB.csv", stringsAsFactors = F)

# preprocessing
DB = DB %>% remSynthErrors()
DB$filename = str_extract_all(DB$runID, "F0[:digit:]+.csv",simplify = T)
ivSPI_sa = left_join(DB, sample_list, by = "filename")

ivSPI_sa = ivSPI_sa %>%
  rename(substrateID = substrateID.x,
         substrateSeq = substrateSeq.x,
         digestTime = digestTime.x) %>%
  select(-substrateID.y,-substrateSeq.y,-digestTime.y)

# no inhibitor
ivSPI_sa$sampleName = str_extract_all(ivSPI_sa$MSfile, "(?<=J_Liepe_)[:alpha:]{2}[:digit:]",simplify = T)
ivSPI_noInh = ivSPI_sa[grepl("A",ivSPI_sa$sampleName),]

# get search results
sru = unique(ivSPI_noInh$filename)
SResults = lapply(sru, function(x){
  
  con = file(paste0("invitroSPI/INPUT/search_results/", x),"r")
  i1 = readLines(con)
  close(con)
  
  start = grep("prot_hit_num",i1)
  currentSearchFile = read.table(paste0("invitroSPI/INPUT/search_results/", x),
                                 skip = start-1,
                                 sep = ",", header = TRUE, fill = TRUE)
  
  kk = which(currentSearchFile$prot_hit_num == "")
  if (length(kk) > 0) {
    kk = min(kk)
    currentSearchFile = currentSearchFile[c(1:(kk-2)),]
  }
  
  return(currentSearchFile)
})

names(SResults) = sru
SResultsTbl = plyr::ldply(SResults)
SResultsTbl$scanNum = sapply(SResultsTbl$pep_scan_title, function(x){
  as.numeric(unlist(strsplit(x, " "))[3])
})

SResultsTbl$pepSeq = SResultsTbl$pep_seq
SResultsTbl$filename = SResultsTbl$.id
SResultsTbl$ionScore = SResultsTbl$pep_score

m_invitroSPI = left_join(ivSPI_noInh, SResultsTbl)
dir.create("data/spectralAngles/")
write.csv(m_invitroSPI, "data/spectralAngles/TSN5noinhibitor_for_spAngles.csv", row.names = F)

m_invitroSPI$MSfile %>% unique()

