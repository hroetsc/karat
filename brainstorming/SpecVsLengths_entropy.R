### protein PCPS ###
# description:  specificity vs lengths - information theory approach
# input:        iaPSIre assigned and quantified peptides
# output:       specificity of the peptides vs SR/intervening sequence length
# author:       HPR

library(dplyr)
library(stringr)
library(seqinr)
library(ggplot2)
library(RColorBrewer)
source("../../proteinsPCPS/new/src/invitroSPI_utils.R")
source("../../proteinsPCPS/new/src/plotting_utils.R")

theme_set(theme_classic())
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(100)

suppressWarnings(dir.create("results/SpecVsLengths/"))

### INPUT ###
load("../../proteinsPCPS/new/data/aSPIre.RData")
load("../../proteinsPCPS/new/data/randomDB.RData")

polypeps = read.csv("../../invitroSPI/revision/results/ProteasomeDB/ProteasomeDB.csv", stringsAsFactors = F)
load("../../invitroSPI/revision/data/wholeDB_random/wholeDB_random_AAMultimapper.RData")

### MAIN PART ###
# ----- preprocessing -----
inSPIRE = Kinetics

proteins = inSPIRE %>%
  ILredundancy() %>%
  disentangleMultimappers.Type() %>%
  disentangleMultimappers.SRlen() %>%
  disentangleMultimappers.IVSeqLen() %>%
  removeMultimappers.Type() %>%
  removeMultimappers.SRlen() %>%
  filter(productType == "PSP") %>%
  uniquePeptides()

proteins_randomDB = rndDB %>%
  ILredundancy() %>%
  disentangleMultimappers.Type() %>%
  disentangleMultimappers.SRlen() %>%
  disentangleMultimappers.IVSeqLen() %>%
  removeMultimappers.Type() %>%
  removeMultimappers.SRlen() %>%
  filter(productType == "PSP") %>%
  uniquePeptides()

polypeps = polypeps %>%
  ILredundancy() %>%
  filterPepLength(cutoff = 5) %>%
  disentangleMultimappers.Type() %>%
  removeMultimappers.Type() %>%
  remSynthErrors() %>%
  filter20Sstandard() %>%
  DBcosmetics() %>%
  filter(productType == "PSP") %>%
  uniquePeptides()

polypeps_randomDB = rbind(randomDB.whole_processed_AA$cis,randomDB.whole_processed_AA$revCis, randomDB.whole_processed_AA$trans)


# ----- relevant information -----
AA = c("A","D","E","F","G","H","K","L","N","P","Q","R","S","T","V","W","Y", "M", "C", "X")
AAchar = c("P","G","A","V","L","M","F","Y","W","H","R","K","D","E","N","Q","S","T","C", "X")

SR1pos = c("P4"=-3, "P3"=-2, "P2"=-1, "P1"=0,
           "P-1"=1, "P-2"=2, "P-3"=3, "P-4"=4)

SR2pos = c("P-4_"=-4, "P-3_"=-3, "P-2_"=-2, "P-1_"=-1,
           "P1_"=0, "P2_"=1, "P3_"=2, "P4_"=3)

PCPpos = c("P4"=-3, "P3"=-2, "P2"=-1, "P1"=0,
           "P1_"=1, "P2_"=2, "P3_"=3, "P4_"=4)

SRnames = c(names(SR1pos), names(SR2pos))
types = c("cis", "revCis", "trans", "PCP")


# ----- get amino acids -----
extract_aminoacids = function(tbl = ""){
  
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


prots_DBaa = extract_aminoacids(proteins)
prots_rndDBaa = extract_aminoacids(proteins_randomDB)

polyp_DBaa = extract_aminoacids(polypeps)
polyp_rndDBaa = extract_aminoacids(polypeps_randomDB)


# ----- add lengths -----
getLen = function(DB) {
  
  DB$pepLen = nchar(DB$pepSeq)
  
  pos = str_split_fixed(DB$positions, pattern = "_", n = Inf)
  pspidx = which(DB$spliceType != "PCP")
  DB$SR1Len = NA
  DB$SR2Len = NA
  DB$SRshortLen = NA
  DB$SR1Len[pspidx] = as.numeric(pos[pspidx, 2]) - as.numeric(pos[pspidx, 1]) + 1
  DB$SR2Len[pspidx] = as.numeric(pos[pspidx, 4]) - as.numeric(pos[pspidx, 3]) + 1
  DB$SRshortLen[pspidx] = do.call(pmin, DB[pspidx,c("SR1Len", "SR2Len")])
  
  DB$IVlen = NA
  cistrans = which(DB$spliceType %in% c("cis","trans"))
  revcis = which(DB$spliceType == "revCis")
  DB$IVlen[c(cistrans,revcis)] = (abs(as.numeric(pos[c(cistrans,revcis), 3]) - as.numeric(pos[c(cistrans,revcis), 2])) - 1) %>%
    as.numeric()
  # DB$IVlen[revcis] = (abs(as.numeric(pos[revcis, 1]) - as.numeric(pos[revcis, 4])) - 1) %>%
  #   as.numeric()
  
  return(DB)
}


prots_DBaa = getLen(prots_DBaa)
prots_rndDBaa = getLen(prots_rndDBaa)

polyp_DBaa = getLen(polyp_DBaa)
polyp_rndDBaa = getLen(polyp_rndDBaa)


save(prots_DBaa, file = "results/SpecVsLengths/DATA_proteins.RData")
save(prots_rndDBaa, file = "results/SpecVsLengths/rndDATA_proteins.RData")
save(polyp_DBaa, file = "results/SpecVsLengths/DATA_polypeps.RData")
save(polyp_rndDBaa, file = "results/SpecVsLengths/rndDATA_polypeps.RData")


# ----- calculate entropy -----

getEntropy = function(DB, rndDB, m, plt = T) {
  
  getFreq = function(df) {
    entropy = apply(df[,SRnames],2,function(x){
      y = table(x)
      names(y)[names(y) == ""] = "X"
      z = same_pos(y)
      z = z[match(AAchar,names(z))]
      # print(names(z))
      return(z)
    }) %>%
      as.matrix()
    rownames(entropy) = AAchar
    
    Hs = apply(entropy,2,function(x){
      x/sum(x)
    })
    return(Hs)
  }
  
  htrue = getFreq(DB)
  hrnd = getFreq(rndDB)
  
  Hs = htrue[AAchar,]/hrnd[AAchar,]
  # Hs = Hs/sum(Hs)
  
  H = apply(Hs,2,function(pos){
    pos = (pos-min(pos, na.rm = T))/(max(pos, na.rm = T) - min(pos, na.rm = T))
    return(-1*sum(pos*log2(pos), na.rm = T))
  })
  
  if (plt) {
    image(Hs %>% t(), col = r,
          xlab = "position", ylab = "amino acid", main = m,
          axes = F)
    axis(1, at = seq(0, 1, 1/15), labels = names(H))
    axis(2, at = seq(0,1,1/19), labels = AAchar)
    
    plot(H, axes = F,
         xlab = "position", ylab = "H [bits]", main = m)
    axis(1, at = seq(1, 16, 1), labels = names(H))
    axis(2)
  }
  
  return(H)
}


# ----- iterate SR lengths -----
iterateLengths = function(DBaa, rndDBaa, type, suffix="") {
  
  
  # --- SR1 lengths
  SR1lengths = DBaa$SR1Len %>% unique() %>% sort()
  png(paste0("results/SpecVsLengths/SR1Len_",type,"_entropy",suffix,".png"), height = 12, width = 14, units = "cm", res = 300)
  par(xpd=T)
  H_SR1 = lapply(SR1lengths, function(x){
    getEntropy(DBaa[which(DBaa$SR1Len == x & DBaa$spliceType %in% type), ],
               rndDBaa[which(rndDBaa$SR1Len == x & rndDBaa$spliceType %in% type), ],
               paste0(type,"-",x), plt=F)
  })
  
  H_SR1m = plyr::ldply(H_SR1) %>% as.matrix()
  H_SR1m[H_SR1m == 0] = NA
  image(H_SR1m %>% t(), col = r,
        xlab = "position", ylab = "SR1 length", main = paste0("SMI: ",type, suffix),
        axes = F, bty = "L")
  legend("topleft", horiz = T, cex = .8,
         legend = c(min(H_SR1m, na.rm = T) %>% round(1), mean(unique(as.vector(H_SR1m)), na.rm = T) %>% round(1), max(H_SR1m, na.rm = T) %>% round(1)), fill = rf(3))
  axis(1, at = seq(0, 1, 1/15), labels = colnames(H_SR1m))
  axis(2, at = seq(0,1,1/(length(SR1lengths)-1)), labels = SR1lengths)
  dev.off()
  
  
  # --- SR2 lengths
  SR2lengths = DBaa$SR2Len %>% unique() %>% sort()
  png(paste0("results/SpecVsLengths/SR2Len_",type,"_entropy",suffix,".png"), height = 12, width = 14, units = "cm", res = 300)
  par(xpd=T)
  H_SR2 = lapply(SR2lengths, function(x){
    getEntropy(DBaa[which(DBaa$SR2Len == x & DBaa$spliceType %in% type), ],
               rndDBaa[which(rndDBaa$SR2Len == x & rndDBaa$spliceType %in% type), ],
               paste0(type,"-",x),plt=F)
  })
  
  H_SR2m = plyr::ldply(H_SR2) %>% as.matrix()
  H_SR2m[H_SR2m == 0] = NA
  image(H_SR2m %>% t(), col = r,
        xlab = "position", ylab = "SR2 length", main = paste0("SMI: ",type, suffix),
        axes = F)
  legend("topleft", horiz = T, cex = .8,
         legend = c(min(H_SR2m, na.rm = T) %>% round(1), mean(unique(as.vector(H_SR2m)), na.rm = T) %>% round(1), max(H_SR2m, na.rm = T) %>% round(1)), fill = rf(3))
  axis(1, at = seq(0, 1, 1/15), labels = colnames(H_SR2m))
  axis(2, at = seq(0,1,1/(length(SR2lengths)-1)), labels = SR2lengths)
  dev.off()
  
  
  # --- IV sequence lengths
  IVlengths = DBaa$IVlen %>% unique() %>% sort()
  IVlengths = IVlengths[IVlengths <= 22]
  png(paste0("results/SpecVsLengths/IVLen_",type,"_entropy",suffix,".png"), height = 12, width = 14, units = "cm", res = 300)
  par(xpd=T)
  H_IV = lapply(IVlengths, function(x){
    getEntropy(DBaa[which(DBaa$IVlen == x & DBaa$spliceType %in% type), ],
               rndDBaa[which(rndDBaa$IVlen == x & rndDBaa$spliceType %in% type), ],
               paste0(type,"-",x),plt=F)
  })
  
  H_IVm = plyr::ldply(H_IV) %>% as.matrix()
  H_IVm[H_IVm == 0] = NA
  image(H_IVm %>% t(), col = r,
        xlab = "position", ylab = "intervening sequence length", main = paste0("SMI: ",type, suffix),
        axes = F)
  legend("topleft", horiz = T, cex = .8,
         legend = c(min(H_IVm, na.rm = T) %>% round(1), mean(unique(as.vector(H_IVm)), na.rm = T) %>% round(1), max(H_IVm, na.rm = T) %>% round(1)), fill = rf(3))
  axis(1, at = seq(0, 1, 1/15), labels = colnames(H_IVm))
  axis(2, at = seq(0,1,1/(length(IVlengths)-1)), labels = IVlengths)
  dev.off()
  
  
}

iterateLengths(prots_DBaa, prots_rndDBaa, type = "cis", suffix = "_proteins")
iterateLengths(prots_DBaa, prots_rndDBaa, type = "revCis", suffix = "_proteins")
iterateLengths(prots_DBaa, prots_rndDBaa, type = "trans", suffix = "_proteins")


iterateLengths(polyp_DBaa, polyp_rndDBaa, type = "cis", suffix = "_polypeps")
iterateLengths(polyp_DBaa, polyp_rndDBaa, type = "revCis", suffix = "_polypeps")
iterateLengths(polyp_DBaa, polyp_rndDBaa, type = "trans", suffix = "_polypeps")

