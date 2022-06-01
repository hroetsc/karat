### protein PCPS ###
# description:  provide functions for random DB quality control
# input:        -
# output:       functions for random DB
# author:       HPR

{
  library(dplyr)
  library(data.table)
  library(tidyr)
  library(stringr)
  library(RColorBrewer)
  library(ggplot2)
  library(ggplotify)
  library(grDevices)
  library(gridExtra)
  library(lattice)
}

source("src/invitroSPI_utils.R")


### MAIN PART ###
# ----- extract amino acids -----
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


extract_aa = function(tbl = ""){
  
  tbl$spliceType[(tbl$spliceType == "") | (is.na(tbl$spliceType))] = "PCP"
  
  # table with position indices
  pos = str_split_fixed(tbl$positions, coll("_"), Inf) %>% as.data.frame()
  
  pcp = which(tbl$spliceType == "PCP")
  psp = which(tbl$spliceType != "PCP")
  
  
  # PCPs
  tbl[pcp, "C"] = str_sub(tbl$substrateSeq[pcp],
                          start = as.numeric(as.character(pos$V2[pcp])),
                          end = as.numeric(as.character(pos$V2[pcp])))
  tbl[pcp, "N"] = str_sub(tbl$substrateSeq[pcp],
                          start = as.numeric(as.character(pos$V1[pcp])),
                          end = as.numeric(as.character(pos$V1[pcp])))
  # cis, revCis and trans
  tbl[psp, "C"] = str_sub(tbl$substrateSeq[psp],
                          start = as.numeric(as.character(pos$V2[psp])),
                          end = as.numeric(as.character(pos$V2[psp])))
  tbl[psp, "N"] = str_sub(tbl$substrateSeq[psp],
                          start = as.numeric(as.character(pos$V3[psp])),
                          end = as.numeric(as.character(pos$V3[psp])))
  
  return(tbl)
}

# ----- matrices of joint probabilities -----
AA = c("A","D","E","F","G","H","K","L","N","P","Q","R","S","T","V","W","Y", "M", "C")
AA = AA[order(AA)]

# joint occurence of a given aa combination at sP1 and sP1'
jointProbs = function(tbl, type) {
  
  # only wanted splice type(s)
  tbl = tbl[tbl$spliceType %in% type, ]
  
  # 20x20 matrix
  m = matrix(ncol = length(AA), nrow = length(AA))
  rownames(m) = AA
  colnames(m) = AA
  
  for (i in 1:length(AA)) {
    for (j in 1:length(AA)) {
      
      m[i, j] = which(AA[i] == tbl$C & AA[j] == tbl$N) %>% length()
      
    }
  }
  
  return(m)
}

# normalise matrices by sum
normaliseMatrices = function(m, freq_matrix) {
  
  print("NORMALISE MATRICES BY BACKGROUND AA FREQUENCY AND BY SUM")
  
  m_norm = m / freq_matrix
  m_norm = m_norm / sum(m)
  
  return(m_norm)
}


# ----- derive amino acid combinations that are impossible -----
impossibleAAcombos = function(subSeqs) {
  
  print("REMOVE IMPOSSIBLE AA COMBINATIONS BASED ON SUBSTRATE STRUCTURE")
  
  m = matrix(ncol = length(AA), nrow = length(AA))
  rownames(m) = AA
  colnames(m) = AA
  
  aacombos = reshape2::melt(m) %>%
    as.data.frame()
  names(aacombos) = c("sP1", "sP1d", "cis_possible")
  aacombos$revCis_possible = NA
  
  for (r in 1:nrow(aacombos)) {
    
    aa1 = aacombos$sP1[r] %>% as.character()
    aa2 = aacombos$sP1d[r] %>% as.character()
    
    for (s in 1:length(subSeqs)) {
      
      cis = str_detect(pattern = paste0(aa1, "[:alpha:]+", aa2), string = subSeqs[s])
      revcis = str_detect(pattern = paste0(aa2, "[:alpha:]*", aa1), string = subSeqs[s])
      
      if (cis) { aacombos$cis_possible[r] = "yes" }
      if (revcis) { aacombos$revCis_possible[r] = "yes"}
      
      if (!is.na(aacombos$cis_possible[r]) & !is.na(aacombos$revCis_possible[r])) {
        break
      }
      
    }
    
  }
  
  aacombos_cis = dcast(aacombos[, c("sP1", "sP1d", "cis_possible")], sP1~sP1d) %>% as.matrix()
  rownames(aacombos_cis) = aacombos_cis[, 1]
  aacombos_cis = aacombos_cis[, c(2:ncol(aacombos_cis))]
  
  aacombos_revcis = dcast(aacombos[, c("sP1", "sP1d", "revCis_possible")], sP1~sP1d) %>% as.matrix()
  rownames(aacombos_revcis) = aacombos_revcis[, 1]
  aacombos_revcis = aacombos_revcis[, c(2:ncol(aacombos_revcis))]
  
  out = list(cis = aacombos_cis,
             revCis = aacombos_revcis)
  
  return(out)
}

# ----- normalise by amino acid frequency in substrate -----
AAfrequencyMatrix = function(subSeqs, iterations = 1e05) {
  
  print("SAMPLE AMINO ACID FREQUENCY FROM THE SUBSTRATES")
  print(paste0("number of sampling iterations: ", iterations))
  
  aacounts = paste(subSeqs, collapse = "") %>%
    strsplit("") %>% 
    unlist() %>%
    table()
  aafreq = aacounts / sum(aacounts)
  aafreq = same_pos(aafreq)
  
  m = matrix(ncol = length(AA), nrow = length(AA))
  m[,] = 0
  rownames(m) = AA
  colnames(m) = AA
  
  set.seed(2478)
  
  for (n in 1:iterations) {
    aa1 = AA[sample(length(AA), 1, prob = aafreq)]
    aa2 = AA[sample(length(AA), 1, prob = aafreq)]
    
    m[rownames(m) == aa1, colnames(m) == aa2] = m[rownames(m) == aa1, colnames(m) == aa2] + 1
    
  }
  
  return(m)
}

# ----- plotting stuff -----
spectralRamp = brewer.pal(11, "Spectral") %>%
  colorRampPalette()
spectral5000 = spectralRamp(5000)

set_matrix_NA = function(rndM, set_na) {
  
  if (nrow(set_na) > 0) {
    for (k in 1:nrow(set_na)) {
      rndM[set_na[k, 1], set_na[k, 2]] = NA
    }
  }
  
  return(rndM)
}


plotHeatmap = function(DB, max_color = 2e-04) {
  
  DB = DB %>%
    ILredundancy() %>%
    disentangleMultimappers.Type() %>%
    removeMultimappers.Type() %>%
    disentangleMultimappers.AA() %>%
    removeMultimappers.AA()
  
  tbl = extract_aa(tbl = DB)
  SubSeqs = DB$substrateSeq %>% unique()
  impossible.aacombos = impossibleAAcombos(subSeqs = SubSeqs)
  aafrequency.matrix = AAfrequencyMatrix(subSeqs = SubSeqs)
  
  rndM.cis = set_matrix_NA(rndM = jointProbs(tbl = tbl, type = "cis") %>%
                             normaliseMatrices(freq_matrix = aafrequency.matrix),
                           set_na = which(is.na(impossible.aacombos$cis), arr.ind = T))
  
  rndM.revCis = set_matrix_NA(rndM = jointProbs(tbl = tbl, type = "revCis") %>%
                             normaliseMatrices(freq_matrix = aafrequency.matrix),
                           set_na = which(is.na(impossible.aacombos$revCis), arr.ind = T))
  
  hm.cis = levelplot(rndM.cis %>% t(),
                 pretty = F,
                 col.regions = spectral5000,
                 # at = seq(0, max_color, length.out = 300),
                 main = "cis",
                 xlab = "sP1'",
                 ylab = "sP1") %>% as.grob()
  
  hm.revCis = levelplot(rndM.revCis %>% t(),
                     pretty = F,
                     col.regions = spectral5000,
                     # at = seq(0, max_color, length.out = 300),
                     main = "revCis",
                     xlab = "sP1'",
                     ylab = "sP1") %>% as.grob()
  
  random.hms = grid.arrange(hm.cis, hm.revCis, ncol = 2, nrow = 1,
                            top = DB$substrateID[1])
  return(random.hms)
}


