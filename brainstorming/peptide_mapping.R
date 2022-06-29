### PCPS brainstorming ###
# description:  mapping of peptide sequences to the substrate - in vitro digestions
# input:        substrate sequence, peptides
# output:       PCPs and PSPs with coordinates and splice types
# author:       HPR, adapted from JL

library(dplyr)
library(stringr)

# for testing
pdb = read.csv("../../invitroSPI/revision/results/ProteasomeDB/ProteasomeDB.csv", stringsAsFactors = F)
subSeq = pdb$substrateSeq[1]
peps = pdb$pepSeq[pdb$substrateSeq == subSeq] %>% unique()

test = pdb[pdb$substrateSeq == subSeq, ]
implement = locate_peps(subSeq, peps = unique(test$pepSeq))

test2 = test %>%
  distinct(pepSeq, .keep_all = T) %>%
  tidyr::separate_rows(spliceType, positions, sep = ";") %>%
  select(pepSeq, positions, spliceType) %>%
  arrange(pepSeq, positions)
test2$spliceType[is.na(test2$spliceType)] = "PCP"

implement2 = implement %>%
  tidyr::separate_rows(spliceType, positions, sep = ";") %>%
  select(pepSeq, positions, spliceType) %>%
  arrange(pepSeq, positions)

k = which(test2$spliceType != implement2$spliceType)

test2[k,]
implement2[k,]

# ----- PCP location -----
locate_PCP <- function(peps, subSeq){
  
  # check which peptides are detected as PCP
  k = sapply(peps, function(x){
    grepl(pattern=x,x=subSeq)
  }) %>% which()
  
  # for those that are, do mapping
  if(length(k)>0){
    
    pcp = lapply(k, function(j){
      
      # get all positions of a PCP to account foe multi-mapping
      pos = str_locate_all(subSeq,peps[j]) %>%
        plyr::ldply() %>%
        as.data.frame()
      id = data.frame(pepSeq = peps[j], pos1 = pos$start, pos2 = pos$end)
      
      return(id)
    }) %>%
      plyr::ldply()
    
  } else {
    pcp <- NULL
  }
  return(pcp %>% select(-.id))
}

# ----- PSP location -----
locate_PSP <- function(PSPpeps, subSeq) {
  
  # sort by N-mers
  Nmers = nchar(PSPpeps) %>% unique() %>% sort()
  
  allPSPpos = list()
  for (k in 1:length(Nmers)) {
    N = Nmers[k]
    
    # get all possible splits of N into two splice-reactants
    q = suppressMessages(tidyr::crossing(c(1:N), c(1:N), .name_repair = "unique"))
    q = q[which(rowSums(q) == N), ] %>%
      as.matrix()
    
    # get all PSP candidates of length N and split them
    cntCand = PSPpeps[nchar(PSPpeps) == N]
    
    PSPpos = lapply(cntCand, function(s){
      
      # get all SRs
      P = strsplit(s,"") %>% unlist()
      cntSRs = sapply(seq(1,nrow(q)), function(i){
        srs = data.frame(pepSeq = s,
                         SR1 = paste(P[1:q[i,1]], collapse = ""),
                         SR2 = paste(P[(q[i,1]+1):N], collapse = ""))
      }) %>% 
        t() %>% 
        as.data.frame()
      
      return(cntSRs)
    }) %>%
      plyr::ldply()
    
    PSPpos = PSPpos %>%
      mutate(pepSeq = unlist(pepSeq),
             SR1 = unlist(SR1),
             SR2 = unlist(SR2))
    
    # map SRs as PCP
    sr1_loc = locate_PCP(peps = PSPpos$SR1, subSeq = subSeq) %>%
      rename(SR1 = pepSeq)
    
    sr2_loc = locate_PCP(peps = PSPpos$SR2, subSeq = subSeq) %>%
      rename(SR2 = pepSeq,
             pos3 = pos1, pos4 = pos2)
    
    POS = suppressMessages(left_join(PSPpos, sr1_loc)) %>%
      na.omit()
    POS = suppressMessages(left_join(POS, sr2_loc)) %>%
      na.omit() %>%
      unique()
    
    # get splice types
    POS$type = NA
    intv = POS$pos3-POS$pos2
    
    POS$type[intv > 0 & POS$pos3 > POS$pos2] = "cis"
    POS$type[intv <= 0 & POS$pos4 < POS$pos1] = "revCis"
    POS$type[intv <= 0 & POS$pos4 >= POS$pos1] = "trans"
    
    # collapse to single assignment per peptide
    POS$allpos = do.call(paste, c(POS[c("pos1","pos2","pos3","pos4")], sep = "_"))
    
    POS = POS %>%
      group_by(pepSeq) %>%
      summarise(spliceType = paste(type, collapse = ";"),
                positions = paste(allpos, collapse = ";"))
    
    allPSPpos[[k]] = POS
  }
  
  pspMAP = plyr::ldply(allPSPpos) %>% as.data.frame()
  return(pspMAP)
}


# ----- actual mapping -----
locate_peps = function(subSeq, peps) {
  
  # account for I/L redundancy
  subSeq = gsub("I","L",subSeq)
  peps = gsub("I","L",peps)
  
  # PCP mapping
  pcpMAP = locate_PCP(peps, subSeq)
  pcpMAP = pcpMAP %>%
    mutate(positions = do.call(paste, c(pcpMAP[c("pos1","pos2")], sep = "_")),
           spliceType = "PCP",
           productType = "PCP") %>%
    select(-pos1, -pos2)
  
  # remove all PSPs from pool of spliced peptides
  PSPpeps = peps[!peps %in% pcpMAP$pepSeq]
  pspMAP = locate_PSP(PSPpeps, subSeq) %>%
    mutate(productType = "PSP")
  
  MAP = rbind(pcpMAP, pspMAP)
  return(MAP)
}

