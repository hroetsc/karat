### karat projetc - PCPS mechanism ###
# description:  show that SR1 is more specific than SR2
# input:        quantitative data set: EGFR, WT sequences of WT/Mut
#               qualitative data set: Roetschke et al. SciData, EGFR, WT sequences of WT/Mut
#               descriptions for Dragon data set
# output:       within-feature variance / diversity of SR1 vs. SR2
# author:       HPR

library(dplyr)
library(stringr)
library(seqinr)
library(protr)
library(Peptides)
library(transport)
library(ggplot2)
library(twosamples)
library(dgof)
library(data.table)
library(bettermc)
library(lsa)
library(readxl)
options(dplyr.summarise.inform = FALSE)
source("src/invitroSPI_utils.R")
source("src/plotting_utils.R")
source("src/_extract-aa.R")

theme_set(theme_classic())

### INPUT ###
load("data/aSPIre.RData")
load("data/invitroSPI.RData")
dragon = readxl::read_excel("results/dragon_molecular_descriptor_list.xlsx") %>%
  as.data.frame()
dragon$Name = str_replace_all(dragon$Name, "^[:digit:]{4} ", "")

### MAIN PART ###
suppressWarnings(dir.create("results/SR1vsSR2/"))

# ----- preprocessing -----
Qual = ProteasomeDB %>%
  ILredundancy() %>%
  disentangleMultimappers.Type() %>%
  disentangleMultimappers.SRlen() %>%
  disentangleMultimappers.IVSeqLen() %>%
  removeMultimappers.SRlen() %>%
  removeMultimappers.AA() %>%
  uniquePeptides()


Quant = Kinetics %>%
  tidyr::separate_rows(digestTimes, intensities, sep=";") %>%
  rename(digestTime = digestTimes,
         intensity = intensities) %>%
  mutate(intensity = as.numeric(intensity),
         digestTime = as.numeric(digestTime)) %>%
  filter(digestTime == 4) %>%
  ILredundancy() %>%
  disentangleMultimappers.Type() %>%
  disentangleMultimappers.SRlen() %>%
  disentangleMultimappers.IVSeqLen() %>%
  removeMultimappers.SRlen() %>%
  removeMultimappers.AA()

extractSRs = function(DB) {
  
  pos = str_split_fixed(DB$positions, "_", Inf)[,c(1:4)]
  pos = apply(pos,2,as.numeric)
  
  DB$sr1 = substr(DB$substrateSeq, pos[,1], pos[,2])
  DB$sr2 = substr(DB$substrateSeq, pos[,3], pos[,4])
  
  return(DB)
}

Qual = extractSRs(Qual)
Quant = extractSRs(Quant)

# ----- number of unique sequences -----
# ..... on qualitative data set

uniQSeqQual = Qual %>%
  filter(! spliceType %in% c("PCP","type_multi-mapper")) %>%
  group_by(substrateID, spliceType) %>%
  summarise(uniQ_SR1 = length(unique(sr1)),
            uniQ_SR2 = length(unique(sr2))) %>%
  tidyr::gather(sr, numSeq, -substrateID, -spliceType)

uniQSeq = ggplot(uniQSeqQual, aes(x = spliceType, y = numSeq, fill = sr)) +
  geom_split_violin() +
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar",  width = .2,
               position = position_dodge(width = .2)) +
  scale_fill_manual(labels = c("SR1","SR2"), values = c("gray", "lightblue")) +
  ggtitle("number of unique sequences")

ggsave(filename = "results/SR1vsSR2/NumUniqueSeq.png", plot = uniQSeq, dpi = "retina", height = 6, width = 6)

uniQSeqQual %>%
  group_by(spliceType, sr) %>%
  summarise(mean = mean(numSeq), median = median(numSeq), std = sd(numSeq)) %>%
  print.data.frame()

# ----- intensity distributions -----
# ..... on quantitative data set

SR1int = Quant %>%
  filter(! spliceType %in% c("PCP","type_multi-mapper")) %>%
  group_by(substrateID, spliceType, sr1) %>%
  summarise(std = sd(log10(intensity+1)),
            cv = sd(log10(intensity+1))/mean(log10(intensity+1))) %>%
  rename(sr = sr1) %>%
  mutate(type = "SR1")

SR2int = Quant %>%
  filter(! spliceType %in% c("PCP","type_multi-mapper")) %>%
  group_by(substrateID, spliceType, sr2) %>%
  summarise(std = sd(log10(intensity+1)),
            cv = sd(log10(intensity+1))/mean(log10(intensity+1))) %>%
  rename(sr = sr2) %>%
  mutate(type = "SR2")

SRint = rbind(SR1int, SR2int)

SDInt = ggplot(SRint, aes(x = spliceType, y = std, fill = type)) +
  geom_split_violin() +
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar",  width = .2,
               position = position_dodge(width = .2)) +
  scale_fill_manual(labels = c("SR1","SR2"), values = c("gray", "lightblue")) +
  ylab("standard deviation log10-transformed intensity") +
  ggtitle("intensity divergence per unique splice-reactant")
SDInt

ggsave(filename = "results/SR1vsSR2/IntensDivergence.png", plot = SDInt, dpi = "retina", height = 6, width = 6)

SRint %>%
  group_by(spliceType, type) %>%
  summarise(n = n(), mean = mean(std, na.rm = T), median = median(std, na.rm = T), std = sd(std, na.rm = T)) %>%
  print.data.frame()

ad_test(a = SR1int$std[SR1int$spliceType == "cis"], b = SR2int$std[SR2int$spliceType == "cis"])
ad_test(a = SR1int$std[SR1int$spliceType == "revCis"], b = SR2int$std[SR2int$spliceType == "revCis"])
ad_test(a = SR1int$std[SR1int$spliceType == "trans"], b = SR2int$std[SR2int$spliceType == "trans"])

# ----- feature mining -----
# ..... on qualitative data set
suppressWarnings(dir.create("results/SR1vsSR2/features/"))


################################################################################
# dbF: vector with feature values
# sr: vector with splice-reactants
getCosineSim = function(dbF, sr) {
  
  # get row indices that belong to each SR
  # (precompute to be faster)
  srU = unique(sr)
  idx = sapply(srU, function(k){
    which(sr == k)
  }) %>% t()
  
  # calculate pairwise feature similarities for each SR
  xx = sapply(seq(1,length(srU)), function(i){
    sapply(seq(i,length(srU)), function(j){
      
      if(i!=j) {
        res = cosine(dbF[idx[i,]],
                     dbF[idx[j,]])
      }
      
    })
  })
  
  xx = unlist(xx) %>% na.omit() %>% as.numeric()
  return(list(xx))
}


getWithinSRsimilarity = function(descriptor, features, cntDS) {
  # extract amino acids for both SRs
  # DB =  Qual[which(Qual$spliceType != "PCP"), ]
  DB = Qual
  DB = left_join(DB, extract_aminoacids(tbl = DB, onlyValidSeq = T))
  
  DBlarge = DB %>%
    select(substrateID, pepSeq, spliceType, positions, sr1, sr2, P4,P3,P2,P1,P1_,P2_,P3_,P4_) %>%
    tidyr::gather(residue, aa, -substrateID, -pepSeq, -spliceType, -positions, -sr1, -sr2)
  
  # get amino acid features
  DBfeat = left_join(DBlarge, descriptor)
  # replace missing amino acids by 0
  DBfeat[is.na(DBfeat)] = 0
  
  # group database into substrateID and spliceType
  types = c("cis","revCis","trans")
  
  SR1 = DBfeat %>%
    filter(spliceType %in% types) %>%
    group_by(substrateID, spliceType) %>%
    filter(residue %in% c("P4","P3","P2","P1")) %>%
    distinct(substrateID, spliceType, sr1, residue, .keep_all = T) %>%
    arrange(residue) %>%  # make sure that residues have the same order
    as.data.table()
  
  SR2 = DBfeat %>%
    filter(spliceType %in% types) %>%
    group_by(substrateID, spliceType) %>%
    filter(residue %in% c("P1_","P2_","P3_","P4_")) %>%
    distinct(substrateID, spliceType, sr2, residue, .keep_all = T) %>%
    arrange(residue) %>%  # make sure that residues have the same order
    as.data.table()
  
  PCP = DBfeat %>%
    filter(spliceType == "PCP") %>%
    group_by(substrateID, spliceType) %>%
    filter(residue %in% c("P4","P3","P2","P1")) %>%
    distinct(substrateID, spliceType, sr1, residue, .keep_all = T) %>%
    arrange(residue) %>%  # make sure that residues have the same order
    as.data.table()
  
  # for each feature: get similarity distributions
  # iterate features
  pdf(paste0("results/SR1vsSR2/features/",cntDS,"_median.pdf"), height = 40, width = 12)
  par(mfrow = c(10,length(types)))
  
  meanVars = lapply(features, function(f){
    print(f)
    
    sr1Var = SR1 %>%
      group_by(substrateID, spliceType) %>%
      rename(cntFeature = f) %>%
      summarise(featDistr = getCosineSim(dbF = cntFeature, sr = sr1))
    
    sr2Var = SR2 %>%
      group_by(substrateID, spliceType) %>%
      rename(cntFeature = f) %>%
      summarise(featDistr = getCosineSim(dbF = cntFeature, sr = sr2))
    
    pcpVar = PCP %>%
      group_by(substrateID) %>%
      rename(cntFeature = f) %>%
      summarise(featDistr = getCosineSim(dbF = cntFeature, sr = sr1))
    
    cntpcp = unlist(pcpVar$featDistr)
    mus = sapply(types, function(t){
      cntsr1 = unlist(sr1Var$featDistr[sr1Var$spliceType == t])
      cntsr2 = unlist(sr2Var$featDistr[sr1Var$spliceType == t])
      boxplot(cntsr1, cntsr2, cntpcp,
              col = c("gray", "lightblue", plottingCols["PCP"]), ylim = c(0,1),
              main = f, sub = paste0(t, " - SR1: ", mean(cntsr1) %>% round(4), ", SR2: ", mean(cntsr2) %>% round(4)))
      
      return(c(SR1 = median(cntsr1,na.rm = T), SR2 = median(cntsr2, na.rm = T), PCP = median(cntpcp, na.rm = T)))
    })
    
    return(mus)
  })
  
  
  dev.off()
  names(meanVars) = features
  save(meanVars, file = paste0("results/SR1vsSR2/features/",cntDS,"_median.RData"))
  
}

################################################################################


### load descriptor datasets
# data(package = "protr")

# topology descriptors
data("AATopo")
force(AATopo)
AATopo$aa = rownames(AATopo)
getWithinSRsimilarity(descriptor = AATopo,
                      features = names(AATopo)[names(AATopo) != "aa"],
                      cntDS = "AATopo")

# molecular properties
data("AAMolProp")
force(AAMolProp)
AAMolProp$aa = rownames(AAMolProp)
getWithinSRsimilarity(descriptor = AAMolProp,
                      features = names(AAMolProp)[names(AAMolProp) != "aa"],
                      cntDS = "AAMolProp")

# 3D molecular structure
data("AAWHIM")
force(AAWHIM)
AAWHIM$aa = rownames(AAWHIM)
getWithinSRsimilarity(descriptor = AAWHIM,
                      features = names(AAWHIM)[names(AAWHIM) != "aa"],
                      cntDS = "AAWHIM")

# AA index database
data("AAindex")
force(AAindex)
AAindex = AAindex[,colnames(AAindex) %in% c("AccNo", AA)] %>%
  t() %>%
  as.data.frame()
colnames(AAindex) = AAindex[1,]
AAindex = AAindex[-1,]
aa = rownames(AAindex)
AAindex = apply(AAindex,2,as.numeric) %>% as.data.frame()
AAindex$aa = aa

getWithinSRsimilarity(descriptor = AAindex,
                      features = names(AAindex)[names(AAindex) != "aa"],
                      cntDS = "AAindex")


# AA composition
data("AAConst")
force(AAConst)
AAConst$aa = rownames(AAConst)
getWithinSRsimilarity(descriptor = AAConst,
                      features = names(AAConst)[names(AAConst) != "aa"],
                      cntDS = "AAConst")

# AA geometric properties
data("AAGeom")
force(AAGeom)
AAGeom$aa = rownames(AAGeom)
getWithinSRsimilarity(descriptor = AAGeom,
                      features = names(AAGeom)[names(AAGeom) != "aa"],
                      cntDS = "AAGeom")

# AA walk
data("AAWalk")
force(AAWalk)
AAWalk$aa = rownames(AAWalk)
getWithinSRsimilarity(descriptor = AAWalk,
                      features = names(AAWalk)[names(AAWalk) != "aa"],
                      cntDS = "AAWalk")

data("AAEigIdx")
force(AAEigIdx)
AAEigIdx$aa = rownames(AAEigIdx)
getWithinSRsimilarity(descriptor = AAEigIdx,
                      features = names(AAEigIdx)[names(AAEigIdx) != "aa"],
                      cntDS = "AAEigIdx")

# ----- load results and plot -----

allDS = c("AAConst","AAEigIdx","AAGeom","AAindex","AAMolProp","AATopo","AAWalk","AAWHIM")

allFeat = list()
counter = 1
for (cntDS in allDS) {
  print(cntDS)
  
  load(paste0("results/SR1vsSR2/features/",cntDS,"_median.RData"))
  meanVarsDF = meanVars %>%
    plyr::ldply() %>%
    mutate(SR = rep(c("SR1","SR2","PCP"), length(meanVars))) %>%
    tidyr::gather(spliceType, meanvar, -.id, -SR)
  
  cnt = ggplot(meanVarsDF, aes(x = spliceType, y = meanvar, fill = factor(SR, levels = c("SR1","SR2","PCP")))) +
    geom_violin(draw_quantiles = 0.5, position = position_dodge(.5)) +
    # stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar",  width = .2,
    #              position = position_dodge(width = .2)) +
    scale_fill_manual("",labels = c("SR1","SR2","PCP"), values = c("gray", "lightblue", plottingCols[["PCP"]])) +
    ylim(c(0,1))+
    ylab("median cosine similarity within SR1/2 features") +
    ggtitle("within-SR feature similarity", subtitle = paste0(cntDS, ", n = ", length(meanVars)))
  
  ggsave(filename = paste0("results/SR1vsSR2/features/_distr",cntDS,"_median.png"), plot = cnt, dpi = "retina", height = 4, width = 4)
  
  meanVarsDF %>%
    group_by(spliceType, SR) %>%
    summarise(n = n(), mean = mean(meanvar), median = median(meanvar), std = sd(meanvar)) %>%
    print.data.frame()
  
  # p-values
  print("cis")
  ad_test(a = meanVarsDF$meanvar[meanVarsDF$spliceType == "cis" & meanVarsDF$SR == "SR1"],
          b = meanVarsDF$meanvar[meanVarsDF$spliceType == "cis" & meanVarsDF$SR == "SR2"]) %>% print()
  print("revCis")
  ad_test(a = meanVarsDF$meanvar[meanVarsDF$spliceType == "revCis" & meanVarsDF$SR == "SR1"],
          b = meanVarsDF$meanvar[meanVarsDF$spliceType == "revCis" & meanVarsDF$SR == "SR2"]) %>% print()
  print("trans")
  ad_test(a = meanVarsDF$meanvar[meanVarsDF$spliceType == "trans" & meanVarsDF$SR == "SR1"],
          b = meanVarsDF$meanvar[meanVarsDF$spliceType == "trans" & meanVarsDF$SR == "SR2"]) %>% print()
  
  allFeat[[counter]] = meanVarsDF
  counter = counter+1
}

allFeat = plyr::ldply(allFeat)

types = c("cis","revCis","trans")

png("results/SR1vsSR2/features/_SR1vsSR2_median.png", height = 5, width = 15, units = "in", res = 300)
par(mfrow = c(1,3))
for (t in types) {
  plot(x = allFeat$meanvar[allFeat$SR == "SR1" & allFeat$spliceType == t],
       y = allFeat$meanvar[allFeat$SR == "SR2" & allFeat$spliceType == t],
       col = plottingCols[t],
       pch = 16, xlim = c(0,1), ylim = c(0,1),
       main = t,
       xlab = "within-SR1 similarity", ylab = "within-SR2 similarity")
  abline(a=0,b=1)
}
dev.off()

## which revCis features are similar across SR1s but not across SR2s
revCisfeat = allFeat[allFeat$spliceType == "revCis", ] %>%
  tidyr::spread(SR, meanvar) %>%
  mutate(diff = SR1-SR2)

data("aaindex")
dragonNames = str_replace_all(dragon$Name, "[:punct:]", "_")
revCisfeat$descriptions = sapply(str_replace_all(revCisfeat$.id, "[:punct:]", "_"), function(x){
  if(x %in% names(aaindex)) {
    return(aaindex[names(aaindex) == x][[1]][["D"]] )
  } else if (x %in% dragonNames) {
    return(dragon$Description[dragonNames == x])
  } else {
    return(NA)
  }
})

### OUTPUT ###
write.csv(revCisfeat[order(revCisfeat$diff, decreasing = T),], "results/SR1vsSR2/features/revCis_featureSim_median.csv",
          row.names = F)

# ----- feature embedding -----

AAindex["I",] = 0
allProps = cbind(AAConst %>% arrange(aa) %>% select(-aa),
                 AAEigIdx %>% arrange(aa) %>% select(-aa),
                 AAGeom %>% arrange(aa) %>% select(-aa),
                 AAindex %>% arrange(aa) %>% select(-aa),
                 AAMolProp %>% arrange(aa) %>% select(-aa),
                 AATopo %>% arrange(aa) %>% select(-aa),
                 AAWalk %>% arrange(aa) %>% select(-aa),
                 AAWHIM %>% arrange(aa))

features = names(allProps)[-which(names(allProps) == "aa")]

DB =  Qual[which(Qual$spliceType != "PCP"), ]
DB = left_join(DB, extract_aminoacids(tbl = DB, onlyValidSeq = T))

DBlarge = DB %>%
  select(substrateID, pepSeq, spliceType, positions, sr1, sr2, P4,P3,P2,P1,P1_,P2_,P3_,P4_) %>%
  tidyr::gather(residue, aa, -substrateID, -pepSeq, -spliceType, -positions, -sr1, -sr2)

# get amino acid features
DBfeat = left_join(DBlarge, allProps)
# replace missing amino acids by 0
DBfeat[is.na(DBfeat)] = 0

SR1 = DBfeat %>%
  filter(spliceType %in% types) %>%
  group_by(substrateID, spliceType) %>%
  filter(residue %in% c("P4","P3","P2","P1")) %>%
  distinct(substrateID, spliceType, sr1, residue, .keep_all = T) %>%
  arrange(residue) %>%  # make sure that residues have the same order
  as.data.table()

SR2 = DBfeat %>%
  filter(spliceType %in% types) %>%
  group_by(substrateID, spliceType) %>%
  filter(residue %in% c("P1_","P2_","P3_","P4_")) %>%
  distinct(substrateID, spliceType, sr2, residue, .keep_all = T) %>%
  arrange(residue) %>%  # make sure that residues have the same order
  as.data.table()

# PCA on each feature
pdf("results/SR1vsSR2/_embeddingAll.pdf", height = 20, width = 12)
par(mfrow = c(5,3))

VarExpl = sapply(features, function(f){
  print(f)
  
  # get current features
  sr1 = SR1 %>%
    rename(cntFeature = f) %>%
    group_by(substrateID,spliceType,sr1) %>%
    summarise(P4 = cntFeature[residue == "P4"],
              P3 = cntFeature[residue == "P3"],
              P2 = cntFeature[residue == "P2"],
              P1 = cntFeature[residue == "P1"]) %>%
    mutate(co = add.alpha("gray", .5))
  
  sr2 = SR2 %>%
    rename(cntFeature = f) %>%
    group_by(substrateID,spliceType,sr2) %>%
    summarise(P4 = cntFeature[residue == "P4_"],  # !!!
              P3 = cntFeature[residue == "P3_"],
              P2 = cntFeature[residue == "P2_"],
              P1 = cntFeature[residue == "P1_"]) %>%
    mutate(co = add.alpha("lightblue", .5))
  both = rbind(sr1,sr2)
  
  # pca
  vare = sapply(types, function(t){
    k1 = which(sr1$spliceType == t)
    k2 = which(sr2$spliceType == t)
    pca1 = prcomp(sr1[k1, c("P4","P3","P2","P1")])
    pca2 = prcomp(sr2[k2, c("P4","P3","P2","P1")])
    pcs = rbind(pca1$x, pca2$x)
    
    plot(pcs[,c(1:2)],
         pch = 16, cex = .8,
         col = c(sr1$co[k1], sr2$co[k2]),
         main = paste0(f, " - ", t),
         sub = paste0("SR1: ", summary(pca1)$importance[3,2]*100, ", SR2: ", summary(pca2)$importance[3,2]*100))
    
    return(rbind(summary(pca1)$importance[3,],summary(pca2)$importance[3,]))
  })
  
  return(vare)
})

dev.off()


# distribution of cumulative variance explained by the first 2 PCs
Var = VarExpl %>%
  t() %>%
  as.data.frame() %>%
  select(V3,V4,V11,V12,V19,V20)
  # select(V1,V2,V9,V10,V17,V18)
names(Var) = c("cis_SR1","cis_SR2","revCis_SR1","revCis_SR2","trans_SR1","trans_SR2")  

Var = Var %>%
  tidyr::gather(type_SR, varexpl) %>%
  mutate(spliceType = str_extract_all(type_SR, "^[:alpha:]+(?=_)", simplify = T),
         SR = str_extract_all(type_SR, "(?<=_)[:alnum:]+", simplify = T))

VarE = ggplot(Var, aes(x = spliceType, y = varexpl, fill = SR)) +
  geom_split_violin() +
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar",  width = .2,
               position = position_dodge(width = .2)) +
  scale_fill_manual(labels = c("SR1","SR2"), values = c("gray", "lightblue")) +
  ylab("cumulative variance explained by PC1+2") +
  ggtitle("variance explained")
VarE
ggsave(filename = "results/SR1vsSR2/VarianceExpl.png", plot = VarE, dpi = "retina", height = 6, width = 6)
