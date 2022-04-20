### PCPS mechanism brainstorming ###
# description:  AA frequency at P2, P2, P3, P4
#               cis, revCis and trans
# input:        ProteasomeDB, random DB (SciData 2022)
# output:       bootstrap intervals for AA frequency
# author:       HR

library(dplyr)
library(stringr)
library(ggplot2)
library(grid)
library(gridExtra)
library(berryFunctions)

source("src/invitroSPI_utils.R")
theme_set(theme_classic())

### HYPERPARAMETRS ###
no_bootstraps = 200
fraction = .8
q = c(0.05, 0.95)
qymin = q[1]
qymax = q[2]


### INPUT ###
DBpath = "~/Documents/QSB/invitroSPI/revision/"
DB = read.csv(paste0(DBpath, "results/ProteasomeDB/ProteasomeDB.csv"),
              stringsAsFactors = F)

# load(paste0(DBpath,"data/wholeDB_random/wholeDB_random_AAMultimapper.RData"))
# randomDB = randomDB.whole_processed_AA
# randomDB = plyr::ldply(randomDB[1:3])
load("data/randomDB_earlyTPs.RData")
randomDB = d

# --------- preprocessing ----------
DB = DB %>%
  ILredundancy() %>%
  filterPepLength() %>%
  disentangleMultimappers.Type() %>%
  removeMultimappers.Type() %>%
  disentangleMultimappers.AA() %>%
  removeMultimappers.AA() %>%
  remSynthErrors() %>%
  filterEarlyTimepoints() %>%
  filter20Sstandard() %>%
  uniquePeptides () %>%
  DBcosmetics()

# !!!
randomDB = randomDB %>%
  ILredundancy() %>%
  filterPepLength() %>%
  disentangleMultimappers.Type() %>%
  removeMultimappers.Type() %>%
  disentangleMultimappers.AA() %>%
  removeMultimappers.AA() %>%
  filterEarlyTimepoints() %>%
  uniquePeptides () %>%
  DBcosmetics()


# ---------- relevant info ----------
# get normalised tables
# all tables should contain the same amino acids
AA = c("A","D","E","F","G","H","K","L","N","P","Q","R","S","T","V","W","Y", "M", "C")
AAchar = c("P","G","A","V","I","L","M","F","Y","W","H","R","K","D","E","N","Q","S","T","C")
# AA = AA[order(AA)]

# positions
SR1pos = c("P4"=-3, "P3"=-2, "P2"=-1, "P1"=0,
           "P-1"=1, "P-2"=2, "P-3"=3, "P-4"=4)

SR2pos = c("P-4_"=-4, "P-3_"=-3, "P-2_"=-2, "P-1_"=-1,
           "P1_"=0, "P2_"=1, "P3_"=2, "P4_"=3)

PCPpos = c("P4"=-3, "P3"=-2, "P2"=-1, "P1"=0,
           "P1_"=1, "P2_"=2, "P3_"=3, "P4_"=4)

SRnames = c(names(SR1pos), names(SR2pos))
types = c("cis", "revCis", "trans", "PCP")

# ---------- extract amino acids and normalise aa frequencies ----------
# extract amino acids at different positions
# P1 (or C-term)
# P1' (or N-term)

extract_aa = function(tbl = ""){
  tbl$spliceType[(tbl$spliceType == "") | (is.na(tbl$spliceType))] = "PCP"
  
  # table with position indices
  pos = str_split_fixed(tbl$positions, coll("_"), Inf) %>% as.data.frame()
  pos = apply(pos, 2, function(x){as.numeric(as.character(x))})
  
  pcp = which(tbl$spliceType == "PCP")
  psp = which(tbl$spliceType != "PCP")
  
  
  # PCPs
  pcpTBL = sapply(PCPpos, function(x){
    substr(tbl$substrateSeq[pcp], start = pos[pcp,2]+x, stop = pos[pcp,2]+x)
  }) %>%
    as.data.frame() %>%
    mutate(spliceType = tbl$spliceType[pcp],
           positions = tbl$positions[pcp],
           pepSeq = tbl$pepSeq[pcp],
           substrateID = tbl$substrateID[pcp],
           substrateSeq = tbl$substrateSeq[pcp])
  
  
  # PSPs
  pspSR1TBL = sapply(SR1pos, function(x){
    substr(tbl$substrateSeq[psp], start = pos[psp,2]+x, stop = pos[psp,2]+x)
  }) %>%
    as.data.frame() %>%
    mutate(spliceType = tbl$spliceType[psp],
           positions = tbl$positions[psp],
           pepSeq = tbl$pepSeq[psp],
           substrateID = tbl$substrateID[psp],
           substrateSeq = tbl$substrateSeq[psp])
  
  
  pspSR2TBL = sapply(SR2pos, function(x){
    substr(tbl$substrateSeq[psp], start = pos[psp,3]+x, stop = pos[psp,3]+x)
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


same_pos = function(tbl) {
  
  if (! all(AA %in% names(tbl))) {
    k = which(! AA %in% names(tbl))
    
    for (i in 1:length(k)) {
      tbl = append(tbl, 0, k[i]-1)
    }
    names(tbl) = AA
  }
  
  return(tbl)
}


# normalised counts
table_aa = function(type, truePos, random) {
  
  if (type %in% truePos$spliceType) {
    
    out = lapply(SRnames, function(x) {
      
      tb = table(truePos[,x])
      tb = tb[! names(tb) %in% c("", NA)]
      
      tbr = table(random[,x])
      tbr = tbr[! names(tbr) %in% c("", NA)]
      
      ind = order(names(tb))
      
      res = same_pos(tb[ind]) / same_pos(tbr[ind])
      res[!is.finite(res)] = 0
      res = res/sum(res)
      
      # res = res[order(match(names(res), AAchar))]
      
      return(res)
    })
    names(out) = SRnames
    return(out) 
    
  } else {
    
    out = lapply(SRnames, function(x){
      rep(NA, length(AA))
    })
    names(out) = SRnames
    
    return(out)
  }
  
}

# ---------- bootstrapping / plotting ----------
# run bootstrapping
bootstrapping = function(DB, randomDB) {
  
  print("extracting amino acids")
  truePos = extract_aa(DB)
  randomPos = extract_aa(randomDB)
  
  cidx = which(randomDB$spliceType == "cis")
  ridx = which(randomDB$spliceType == "revCis")
  tidx = which(randomDB$spliceType == "trans")
  pidx = which(randomDB$spliceType %in% c("PCP", NA))
  
  # store bootstrap results
  # list of arrays: every list entry corresponds to one position
  # dim1: product/splice type, dim2: AAs, dim3: bootstraps
  MASTER = lapply(SRnames, function(x){
    lst = array(dim = c(4, length(AA), no_bootstraps),
                  dimnames = list(c("cis", "revCis", "trans", "PCP"),
                                  AA,
                                  c(1:no_bootstraps)))
  })
  names(MASTER) = SRnames
  
  
  print("bootstrapping")
  
  pb = txtProgressBar(min = 0, max = no_bootstraps, style = 3)
  for (i in 1:no_bootstraps) {
    setTxtProgressBar(pb,i)
    
    # do boostrapping for all positions
    cntTruePos = truePos[sample(nrow(truePos), ceiling(nrow(truePos) * fraction)), ]
    
    cis = table_aa(type = "cis", cntTruePos[cntTruePos$spliceType=="cis",], randomPos[cidx,])
    revCis = table_aa(type = "revCis", cntTruePos[cntTruePos$spliceType=="revCis",], randomPos[ridx,])
    trans = table_aa(type = "trans", cntTruePos[cntTruePos$spliceType=="trans",], randomPos[tidx,])
    pcp = table_aa(type = "PCP", cntTruePos[cntTruePos$spliceType=="PCP",], randomPos[pidx,])
    
    # iterate positions
    for (j in SRnames) {
      MASTER[[j]][1:dim(MASTER[[1]])[1],1:dim(MASTER[[1]])[2],i] = rbind(cis[[j]],
                                                                         revCis[[j]],
                                                                         trans[[j]],
                                                                         pcp[[j]])
    }
    
  }
  
  # add peptide counts
  n = truePos %>% 
    group_by(spliceType) %>% 
    dplyr::count()
  
  MASTER[[length(MASTER)+1]] = n
  
  return(MASTER)
}


db = bootstrapping(DB, randomDB)

# ---------- plotting position-wise ----------
# compare databases with each other
# for each amino acid
plotAAfreqs_pos = function(position, db) {
  
  # format data
  n = db[[length(SRnames)+1]]
  
  # 4-dimensional array: product type, amino acids, bootstrap iterations, positions
  db2 = l2array(db[1:length(SRnames)])
  
  # iterate types and get frequency
  FREQS = numeric()
  k1 = which(dimnames(db2)[[4]] == position)
  
  for (t in types) {
    k2 = which(dimnames(db2)[[1]] == t)
    
    freqs = db2[k2,,,k1] %>% 
      t() %>% 
      as.data.frame() %>% 
      tidyr::gather() %>%
      mutate(type = t) %>%
      as.data.frame()
    
    FREQS = rbind(FREQS, freqs)
  }
  
  names(FREQS) = c("aa", "frequency","type")
  FREQS = FREQS[is.finite(FREQS$frequency), ]
  FREQS$type = factor(FREQS$type, levels = types)
  
  # get quantiles (bootstrap confidence interval)
  # plot: boxplot with confidence interval in box
  
  # quantiles = apply(db[[terminus]][k,,], MARGIN = 1, quantile, q)
  
  f = FREQS %>%
    dplyr::group_by(aa, type) %>%
    dplyr::summarise(lower = quantile(frequency, q[1]),
              upper = quantile(frequency, q[2]),
              middle = quantile(frequency, 0.5),
              IQR = diff(c(lower, upper)),
              ymin = max(quantile(frequency, qymin), lower - 1.5 * IQR),
              ymax = min(quantile(frequency, qymax), upper + 1.5 * IQR),
              outliers = list(frequency[which(frequency > upper + 1.5 * IQR | 
                                                frequency < lower - 1.5 * IQR)])) %>%
    dplyr::ungroup()
  
  f$aa = factor(f$aa, levels = AAchar[-which(AAchar == "I")])
  aaFreqs = f %>%
    ggplot(aes(x = aa, fill = type)) +
    geom_boxplot(aes(lower = lower,
                     upper = upper,
                     middle = middle,
                     ymin = ymin,
                     ymax = ymax),
                 stat = "identity",
                 alpha = .6,
                 # outlier.shape = NA
                 position = position_dodge(width = .5)) +
    scale_fill_manual("product type",
                      values = c(plottingCols[["cis"]], plottingCols[["revCis"]],
                                 plottingCols[["trans"]], plottingCols[["PCP"]])) +
    xlab("amino acid") +
    ylab("normalised frequency") +
    ylim(0,.3) +
    ggtitle(paste0("position: ", position),
            subtitle = paste0("n = ",sum(n$n)," peptides --- ",
                              no_bootstraps, " bootstrap iterations with ", fraction*100, "% of the data, ",
                              (q[2] - q[1])*100, "% confidence interval"))
  aaFreqs
  
  return(aaFreqs)
}


pdf("results/AAfreq_position-wise.pdf", height=4, width=10)
for (p in SRnames) {
  plotAAfreqs_pos(position = p, db) %>%
    print()
}
dev.off()


# ---------- plotting aa-wise ----------
# compare databases with each other
# for each position

plotAAfreqs_aa = function(aa, db) {
  
  # format data
  n = db[[length(SRnames)+1]]
  
  # 4-dimensional array: product type, amino acids, bootstrap iterations, positions
  db2 = l2array(db[1:length(SRnames)])
  
  # iterate types and get frequency
  FREQS = numeric()
  k1 = which(dimnames(db2)[[2]] == aa)
  
  for (t in types) {
    k2 = which(dimnames(db2)[[1]] == t)
    
    freqs = db2[k2,k1,,] %>% 
      as.data.frame() %>% 
      tidyr::gather() %>%
      mutate(type = t) %>%
      as.data.frame()
    
    FREQS = rbind(FREQS, freqs)
  }
  
  names(FREQS) = c("position", "frequency","type")
  FREQS = FREQS[is.finite(FREQS$frequency), ]
  FREQS$type = factor(FREQS$type, levels = types)
  
  # get quantiles (bootstrap confidence interval)
  # plot: boxplot with confidence interval in box
  
  # quantiles = apply(db[[terminus]][k,,], MARGIN = 1, quantile, q)
  
  f = FREQS %>%
    dplyr::group_by(position, type) %>%
    dplyr::summarise(lower = quantile(frequency, q[1]),
                     upper = quantile(frequency, q[2]),
                     middle = quantile(frequency, 0.5),
                     IQR = diff(c(lower, upper)),
                     ymin = max(quantile(frequency, qymin), lower - 1.5 * IQR),
                     ymax = min(quantile(frequency, qymax), upper + 1.5 * IQR),
                     outliers = list(frequency[which(frequency > upper + 1.5 * IQR | 
                                                       frequency < lower - 1.5 * IQR)])) %>%
    dplyr::ungroup()
  f$position = factor(f$position, levels = SRnames)
  
  aaFreqs = f %>%
    ggplot(aes(x = position, fill = type)) +
    geom_boxplot(aes(lower = lower,
                     upper = upper,
                     middle = middle,
                     ymin = ymin,
                     ymax = ymax),
                 stat = "identity",
                 alpha = .6,
                 # outlier.shape = NA
                 position = position_dodge(width = .5)) +
    geom_vline(xintercept = c(4.5, 8.5, 12.5), linetype="dashed", color="gray") +
    scale_fill_manual("product type",
                      values = c(plottingCols[["cis"]], plottingCols[["revCis"]],
                                 plottingCols[["trans"]], plottingCols[["PCP"]])) +
    xlab("position") +
    ylab("normalised frequency") +
    ylim(0,.3) +
    ggtitle(paste0("amino acid: ", aa),
            subtitle = paste0("n = ",sum(n$n)," peptides --- ",
                              no_bootstraps, " bootstrap iterations with ", fraction*100, "% of the data, ",
                              (q[2] - q[1])*100, "% confidence interval"))
  aaFreqs
  
  return(aaFreqs)
}


pdf("results/AAfreq_aa-wise.pdf", height=4, width=10)
for (a in AAchar[-which(AAchar == "I")]) {
  plotAAfreqs_aa(aa = a, db) %>%
    print()
}
dev.off()


# ---------- output ----------
save(db, file = "results/AAfreq.RData")


