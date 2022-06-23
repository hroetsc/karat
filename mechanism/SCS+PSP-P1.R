### protein PCPS ###
# description:  provide functions for SCS-P1 and PSP-P1
# input:        qiSPI identified sequences, protein sequences
# output:       cleavage and splicing strength for each residue
# author:       HPR


library(dplyr)
library(stringr)
library(seqinr)

source("src/invitroSPI_utils.R")

# ---------- relevant info ----------
# positions
SR1pos = c("P6"=-5,"P5"=-4,"P4"=-3, "P3"=-2, "P2"=-1, "P1"=0,
           "P-1"=1, "P-2"=2, "P-3"=3, "P-4"=4, "P-5"=5, "P-6"=6)

SR2pos = c("P-6_"=-6,"P-5_"=-5,"P-4_"=-4, "P-3_"=-3, "P-2_"=-2, "P-1_"=-1,
           "P1_"=0, "P2_"=1, "P3_"=2, "P4_"=3, "P5_"=4, "P6_"=5)

PCPpos = c("P6"=-5,"P5"=-4,"P4"=-3, "P3"=-2, "P2"=-1, "P1"=0,
           "P-1"=1, "P-2"=2, "P-3"=3, "P-4"=4, "P-5"=5, "P-6"=6,
           "P1_"=1, "P2_"=2, "P3_"=3, "P4_"=4, "P5_"=5, "P6_"=6)

SRnames = c(names(SR1pos), names(SR2pos))
types = c("cis", "revCis", "trans", "PCP")


# ---------- extract amino acids ----------
# extract coordinates of different positions in the substrate
extract_coordinates = function(tbl = ""){
  
  tbl$spliceType[(tbl$spliceType == "") | (is.na(tbl$spliceType))] = "PCP"
  
  # table with position indices
  pos = str_split_fixed(tbl$positions, coll("_"), Inf) %>% as.data.frame()
  pos = apply(pos, 2, function(x){as.numeric(as.character(x))})
  
  pcp = which(tbl$spliceType == "PCP")
  psp = which(tbl$spliceType != "PCP")
  
  
  # PCPs
  pcpTBL = sapply(PCPpos, function(x){
    # substr(tbl$substrateSeq[pcp], start = pos[pcp,2]+x, stop = pos[pcp,2]+x)
    pos[pcp,2]+x
  }) %>%
    as.data.frame() %>%
    mutate(spliceType = tbl$spliceType[pcp],
           positions = tbl$positions[pcp],
           pepSeq = tbl$pepSeq[pcp],
           substrateID = tbl$substrateID[pcp],
           substrateSeq = tbl$substrateSeq[pcp])
  
  
  # PSPs
  pspSR1TBL = sapply(SR1pos, function(x){
    # substr(tbl$substrateSeq[psp], start = pos[psp,2]+x, stop = pos[psp,2]+x)
    pos[psp,2]+x
  }) %>%
    as.data.frame() %>%
    mutate(spliceType = tbl$spliceType[psp],
           positions = tbl$positions[psp],
           pepSeq = tbl$pepSeq[psp],
           substrateID = tbl$substrateID[psp],
           substrateSeq = tbl$substrateSeq[psp])
  
  
  pspSR2TBL = sapply(SR2pos, function(x){
    # substr(tbl$substrateSeq[psp], start = pos[psp,3]+x, stop = pos[psp,3]+x)
    pos[psp,3]+x
  }) %>%
    as.data.frame() %>%
    mutate(spliceType = tbl$spliceType[psp],
           positions = tbl$positions[psp],
           pepSeq = tbl$pepSeq[psp],
           substrateID = tbl$substrateID[psp],
           substrateSeq = tbl$substrateSeq[psp])
  
  # merge all tables
  pspTBL = cbind(pspSR1TBL[,names(SR1pos)], pspSR2TBL)
  
  pcpPlaceholder = matrix("", length(pcp), length(SR2pos)) %>%
    as.data.frame()
  names(pcpPlaceholder) = SRnames[! SRnames %in% names(PCPpos)]
  pcpTBL2 = cbind(cbind(pcpTBL[,names(PCPpos)], pcpPlaceholder)[, c(names(SR1pos), names(SR2pos))],
                  pcpTBL[,which(!names(pcpTBL) %in% names(PCPpos))])
  
  TBL = rbind(pspTBL, pcpTBL2) %>% as.data.frame()
  return(TBL)
}


# ----- resolve multi-mappers -----
# assign weight to position multi-mappers
resolve_multimapper = function(ProteasomeDB) {
  
  k = which(str_detect(ProteasomeDB$positions, coll(";")))
  if (length(k) > 0) {
    
    DB_mm = ProteasomeDB[k, ]
    DB = ProteasomeDB[-k, ]
    
    for (r in 1:nrow(DB_mm)) {
      cnt = DB_mm[r, ] %>%
        tidyr::separate_rows(positions, sep=";") %>%
        as.data.frame()
      
      cnt$intensity = cnt$intensity / nrow(cnt)
      
      DB = rbind(DB, cnt)
    }
    
    return(as.data.frame(DB))
    
  } else {
    return(ProteasomeDB)
  }
  
}


# ----- compute SCS-P1 and PSP-P1 -----

SCS_and_PSP = function(DBMaster,target,meanOverBio=T,Zscale=F,rawVals=F,SR2forSCS=T) {
  
  pcpidx = which(DBMaster$productType == "PCP")
  pspidx = which(DBMaster$productType == "PSP")
  
  res = str_split_fixed(DBMaster$positions, coll("_"), n = Inf) %>%
    as.data.frame()
  
  d = data.frame(residue = c(1:nchar(DBMaster$substrateSeq[1])),
                 scs_mean = 0,
                 scs_sd = 0,
                 scs_n = 0,
                 psp_mean = 0,
                 psp_sd = 0,
                 psp_n = 0)
  
  for (i in 1:nrow(d)) {
    
    kpcp = which(DBMaster[,target] == d$residue[i] & DBMaster$productType == "PCP")
    kpsp = which(DBMaster[,target] == d$residue[i] & DBMaster$productType == "PSP")
    # add the C-terminus of spliced peptides as cleaved residue
    if (target == "P1" & SR2forSCS) {
      kpcp = c(kpcp, which(res$V4[pspidx] == d$residue[i]))
    }
    
    cntscs = DBMaster[kpcp,] %>%
      dplyr::group_by(biological_replicate) %>%
      dplyr::summarise(int = sum(intensity))
    
    d$scs_mean[i] = paste(cntscs$int, collapse = "_")
    d$scs_n[i] = DBMaster$pepSeq[kpcp] %>%
      unique() %>%
      length()
    
    cntpsp = DBMaster[kpsp,] %>%
      dplyr::group_by(biological_replicate) %>%
      dplyr::summarise(int = sum(intensity))
    
    d$psp_mean[i] = paste(cntpsp$int, collapse = "_")
    d$psp_n[i] = DBMaster$pepSeq[kpsp] %>%
      unique() %>%
      length()
  }
  
  d[is.na(d)] = 0
  d[d == ""] = 0
  
  # remove substrate's C terminus from the SCS
  d$scs_mean[nchar(DBMaster$substrateSeq[1])] = 0
  d$scs_sd[nchar(DBMaster$substrateSeq[1])] = 0
  
  # normalise by sum for each biological replicate
  scs = apply(str_split(d$scs_mean, pattern = "_", simplify = T), 2, as.numeric) %>%
    as.data.frame()
  
  psp = apply(str_split(d$psp_mean, pattern = "_", simplify = T), 2, as.numeric) %>%
    as.data.frame()
  
  if (meanOverBio) {
    if (!Zscale) {
      scs = sweep(scs, 2, colSums(scs,na.rm = T), FUN = "/")
      psp = sweep(psp, 2, colSums(psp,na.rm = T), FUN = "/")
    } else if (rawVals) {
      
      scs = scs
      psp = psp
      
    } else {
      scs = apply(scs,2,function(x){
        return((x - min(x)) / (max(x)-min(x)))
      })
      psp = apply(psp,2,function(x){
        return((x - min(x)) / (max(x)-min(x)))
      })
    }
    
    d$scs_mean = rowMeans(scs,na.rm = T) * 100
    d$scs_sd = apply(scs,1,sd,na.rm=T) * 100
    
    d$psp_mean = rowMeans(psp,na.rm = T) * 100
    d$psp_sd = apply(psp,1,sd,na.rm=T) * 100
    
  } else {
    if (!Zscale) {
      scs = sweep(scs, 2, colSums(scs,na.rm = T), FUN = "/") * 100
      psp = sweep(psp, 2, colSums(psp,na.rm = T), FUN = "/") * 100
    } else if (rawVals) {
      
      scs = scs
      psp = psp
      
    } else {
      scs = apply(scs,2,function(x){
        return(((x - min(x,na.rm = T)) / (max(x,na.rm = T)-min(x,na.rm = T)))*100)
      })
      psp = apply(psp,2,function(x){
        return(((x - min(x,na.rm = T)) / (max(x,na.rm = T)-min(x,na.rm = T)))*100)
      })
    }
    
    scs[is.na(scs)] = 0
    d$scs_mean = apply(scs,1,paste,collapse = ";")
    d$scs_sd = 0
    
    psp[is.na(psp)] = 0
    d$psp_mean = apply(psp,1,paste,collapse = ";")
    d$psp_sd = 0
  }
  
  d[is.na(d)] = 0
  return(d)
}

# ----- compute STS -----
# substrate-specific transpeptidation strength
STS = function(DBMaster, meanOverBio=T) {
  
  tbl = DBMaster[DBMaster$productType == "PSP", ]
  subSeq = tbl$substrateSeq[1]
  
  res = str_split_fixed(DBMaster$positions, coll("_"), n = Inf) %>%
    as.data.frame()
  
  d = tidyr::crossing(c(1:nchar(subSeq)), c(1:nchar(subSeq)), .name_repair = "unique") %>%
    as.data.frame()
  names(d) = c("P1","P1_")
  
  d = d %>%
    mutate(sts_mean = 0,
           sts_sd = 0,
           sts_n = 0)
  
  for (i in 1:nrow(d)) {
    
    kpsp = which(DBMaster$P1 == d$P1[i] & DBMaster$P1_ == d$P1_[i])
    
    cntsts = DBMaster[kpsp, ] %>%
      dplyr::group_by(biological_replicate) %>%
      dplyr::summarise(int = sum(intensity))
    
    d$sts_mean[i] = paste(cntsts$int, collapse = "_")
    d$sts_n[i] = DBMaster$pepSeq[kpsp] %>%
      unique() %>%
      length()
  }
  
  # normalise by sum for each biological replicate
  sts = apply(str_split(d$sts_mean, pattern = "_", simplify = T), 2, as.numeric) %>%
    as.data.frame()
  sts = sweep(sts, 2, colSums(sts,na.rm = T), FUN = "/") * 100
  
  if (meanOverBio) {
    
    d$sts_mean = rowMeans(sts,na.rm = T) * 100
    d$sts_sd = apply(sts,1,sd,na.rm=T) * 100
    
  } else {
    
    sts[is.na(sts)] = 0
    d$sts_mean = apply(sts,1,paste,collapse = ";")
    d$sts_sd = 0
    
  }
  
  d[is.na(d) | d == ""] = 0
  return(d)
}

# ----- iterate substrate sequences -----

SCSandPSP_allSubs = function(DB, target,meanOverBio=T,Zscale=F,rawVals=F,SR2forSCS=T) {
  
  subs = DB$substrateID %>% unique()
  
  print("resolving multimappers")
  
  print("extracting coordinates")
  DBaa = extract_coordinates(DB)
  DBMaster = inner_join(DB, DBaa, by=c("pepSeq", "positions")) %>%
    unique() %>%
    rename(spliceType = spliceType.x,
           substrateID = substrateID.x,
           substrateSeq = substrateSeq.x)
  
  DBMaster$productType = toupper(DBMaster$productType)
  
  print("calculating SCS and PSP-P1")
  out = lapply(subs, function(x){
    SCS_and_PSP(DBMaster[DBMaster$substrateID == x, ], target,meanOverBio,Zscale,rawVals,SR2forSCS)
  })
  
  names(out) = subs
  return(out)
}


STS_allSubs = function(DB, meanOverBio = T) {
  subs = DB$substrateID %>% unique()
  
  print("resolving multimappers")
  
  print("extracting coordinates")
  DBaa = extract_coordinates(DB)
  DBMaster = inner_join(DB, DBaa, by=c("pepSeq", "positions")) %>%
    unique() %>%
    rename(spliceType = spliceType.x,
           substrateID = substrateID.x,
           substrateSeq = substrateSeq.x)
  
  DBMaster$productType = toupper(DBMaster$productType)
  
  print("calculating SCS and PSP-P1")
  out = lapply(subs, function(x){
    STS(DBMaster[DBMaster$substrateID == x, ],meanOverBio)
  })
  
  names(out) = subs
  return(out)
}

