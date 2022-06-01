### karat projetc - PCPS mechanism ###
# description:  determine differences between SR1 and PCP (precursors)
#               qualitative data set: Roetschke et al. SciData, EGFR, WT sequences of WT/Mut
# output:       similarity between SR1s and PCPs for different features
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
library(uwot)
source("src/invitroSPI_utils.R")
source("src/plotting_utils.R")
source("src/_extract-aa.R")

theme_set(theme_classic())

### INPUT ###
load("data/invitroSPI.RData")

dragon = readxl::read_excel("results/dragon_molecular_descriptor_list.xlsx") %>%
  as.data.frame()
dragon$Name = str_replace_all(dragon$Name, "^[:digit:]{4} ", "")

### MAIN PART ###
suppressWarnings(dir.create("results/SR1vsPCP/"))
suppressWarnings(dir.create("results/SR1vsPCP/features/"))


# ----- preprocessing -----
Qual = ProteasomeDB %>%
  ILredundancy() %>%
  disentangleMultimappers.Type() %>%
  disentangleMultimappers.SRlen() %>%
  disentangleMultimappers.IVSeqLen() %>%
  removeMultimappers.SRlen() %>%
  removeMultimappers.AA() %>%
  uniquePeptides()

extractSRs = function(DB) {
  
  pos = str_split_fixed(DB$positions, "_", Inf)[,c(1:4)]
  pos = apply(pos,2,as.numeric)
  
  DB$sr1 = substr(DB$substrateSeq, pos[,1], pos[,2])
  DB$sr2 = substr(DB$substrateSeq, pos[,3], pos[,4])
  
  return(DB)
}

Qual = extractSRs(Qual)

# ----- feature mining -----
# ..... on qualitative data set

###############################################################################
# dbF: vector with feature values
# sr: vector with splice-reactants
getCosineSim = function(PSPdbF, PCPdbF, sr, pcp) {
  
  if (length(sr) > 0 & length(pcp) > 0) {
    # get row indices that belong to each SR
    # (precompute to be faster)
    srU = unique(sr)
    sridx = lapply(srU, function(k){
      which(sr == k)
    })
    
    pcpU = unique(pcp)
    pcpidx = lapply(pcpU, function(k){
      which(pcp == k)
    })
    
    # calculate pairwise feature similarities for each SR
    xx = sapply(seq(1,length(srU)), function(i){
      sapply(seq(1,length(pcpU)), function(j){
        res = cosine(PSPdbF[sridx[i][[1]]],
                     PCPdbF[pcpidx[j][[1]]])
        
      })
    })
    
    xx = unlist(xx) %>% na.omit() %>% as.numeric()
    return(xx)
  } else {
    return(NA)
  }
  
}


getFeatureSimilarity = function(descriptor, features, cntDS) {
  # extract amino acids
  DB =  Qual
  DB = left_join(DB, extract_aminoacids(tbl = DB, onlyValidSeq = T))
  
  # get the P-1
  pos = str_split_fixed(DB$positions,"_",Inf)[,c(1:4)]
  DB$`P-1` = substr(DB$substrateSeq, as.numeric(pos[,2])+1, as.numeric(pos[,2])+1)
  
  DBlarge = DB %>%
    select(substrateID, pepSeq, spliceType, positions, sr1, sr2, P4,P3,P2,P1,`P-1`) %>%
    tidyr::gather(residue, aa, -substrateID, -pepSeq, -spliceType, -positions, -sr1, -sr2)
  
  # get amino acid features
  DBfeat = left_join(DBlarge, descriptor)
  # replace missing amino acids by 0
  DBfeat[is.na(DBfeat)] = 0
  
  # group database into substrateID and spliceType
  types = c("cis","revCis","trans")
  
  SR1 = DBfeat %>%
    filter(spliceType %in% types) %>%
    filter(residue %in% c("P4","P3","P2","P1","P-1")) %>%
    # filter(residue %in% c("P1","P-1")) %>%
    distinct(substrateID, spliceType, sr1, residue, .keep_all = T) %>%
    group_by(pepSeq, sr1) %>%
    mutate(valid = ifelse(length(residue) == 5,T,F)) %>%
    # mutate(valid = ifelse(length(residue) == 2,T,F)) %>%
    filter(valid) %>% ungroup() %>%
    select(-valid) %>%
    arrange(residue) %>%  # make sure that residues have the same order
    as.data.table()
  
  PCP = DBfeat %>%
    filter(spliceType == "PCP") %>%
    filter(residue %in% c("P4","P3","P2","P1","P-1")) %>%
    # filter(residue %in% c("P1","P-1")) %>%
    distinct(substrateID, spliceType, sr1, residue, .keep_all = T) %>%
    group_by(pepSeq, sr1) %>%
    mutate(valid = ifelse(length(residue) == 5,T,F)) %>%
    # mutate(valid = ifelse(length(residue) == 2,T,F)) %>%
    filter(valid) %>% ungroup() %>%
    select(-valid) %>%
    arrange(residue) %>%  # make sure that residues have the same order
    as.data.table()
  
  # split both data frames into substrate IDs and splice types
  subs = DB$substrateID %>% unique()
  
  ALL = lapply(types, function(t){
    cnt = lapply(subs, function(s) {
      return(list(PSP = SR1[which(SR1$substrateID == s & SR1$spliceType == t), ],
                  PCP = PCP[which(PCP$substrateID == s), ]))
    })
    names(cnt) = subs
    return(cnt)
  })
  names(ALL) = types
  
  # for each feature: get similarity distributions
  # iterate features
  pdf(paste0("results/SR1vsPCP/features/",cntDS,"_norm_median.pdf"), height = 40, width = 12)
  par(mfrow = c(10,length(types)))
  
  meanVars = lapply(features, function(f){
    print(f)
    
    meanVarTypes = lapply(ALL, function(cntType){
      
      meanVarSubs = lapply(cntType, function(cntSubs){
        
        # print(cntSubs$PCP$substrateID[1])
        xx = getCosineSim(PSPdbF = as.numeric(unlist(cntSubs$PSP[,..f])),
                          PCPdbF = as.numeric(unlist(cntSubs$PCP[,..f])),
                          sr = cntSubs$PSP$sr1,
                          pcp = cntSubs$PCP$sr1)
        
        return(unlist(xx))
      })
      
      return(unlist(meanVarSubs) %>% na.omit() %>% as.numeric())
    })
    
    mus = sapply(meanVarTypes, median, na.rm = T)
    
    boxplot(meanVarTypes, main = f,
            col = c(plottingCols["cis"], plottingCols["revCis"], plottingCols["trans"]),
            sub = paste0("cis: ", round(mus[1],4), ", revCis: ", round(mus[2],4), ", trans: ", round(mus[3],4)))
    
    return(mus)
    
  })
  
  dev.off()
  names(meanVars) = features
  save(meanVars, file = paste0("results/SR1vsPCP/features/",cntDS,"_norm_median.RData"))
  
}

################################################################################

### load descriptor datasets
# all Dragon descriptors
data("AADescAll")
force(AADescAll)
aa = rownames(AADescAll)
AADescAllnorm = apply(as.data.frame(AADescAll),2,function(x){
  (x - min(x))/(max(x) - min(x))
}) %>% as.data.frame()
AADescAllnorm$aa = aa
AADescAll$aa = aa
getFeatureSimilarity(descriptor = AADescAllnorm,
                     features = names(AADescAllnorm)[names(AADescAllnorm) != "aa"],
                     cntDS = "AADescAll")


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
AAindex[is.na(AAindex)] = 0
AAindexnorm = apply(as.data.frame(AAindex),2,function(x){
  (x - min(x))/(max(x) - min(x))
}) %>% as.data.frame()
AAindexnorm$aa = aa
AAindex$aa = aa

getFeatureSimilarity(descriptor = AAindexnorm,
                      features = names(AAindexnorm)[names(AAindexnorm) != "aa"],
                      cntDS = "AAindex")

# ---- load feature similarities -----

allDS = c("AADescAll","AAindex")

allFeat = list()
counter = 1
for (cntDS in allDS) {
  print(cntDS)
  
  load(paste0("results/SR1vsPCP/features/",cntDS,"_norm_median.RData"))
  meanVarsDF = meanVars %>%
    plyr::ldply() %>%
    tidyr::gather(spliceType, similarity, -.id)
  
  cnt = ggplot(meanVarsDF, aes(x = spliceType, y = similarity, fill = spliceType)) +
    geom_violin(draw_quantiles = c(0.25,0.5,0.75)) +
    scale_fill_manual(values = c(plottingCols[["cis"]], plottingCols[["revCis"]], plottingCols[["trans"]])) +
    ylab("mean cosine vector similarity") +
    ggtitle("PCP/PSP precursor feature similarity", subtitle = paste0(cntDS, ", n = ", length(meanVars)))
  
  ggsave(filename = paste0("results/SR1vsPCP/features/_distr",cntDS,"_norm_median.png"),
         plot = cnt, dpi = "retina", height = 6, width = 6)
  
  meanVarsDF %>%
    group_by(spliceType) %>%
    summarise(n = n(), mean = mean(similarity, na.rm = T),
              median = median(similarity, na.rm = T),
              std = sd(similarity, na.rm = T)) %>%
    print.data.frame()
  
  allFeat[[counter]] = meanVarsDF
  counter = counter+1
}

allFeatDF = plyr::ldply(allFeat) %>%
  tidyr::spread(spliceType, similarity)

data("aaindex")

dragonNames = str_replace_all(dragon$Name, "[:punct:]", "_")
allFeatDF$descriptions = sapply(str_replace_all(allFeatDF$.id, "[:punct:]", "_"), function(x){
  if(x %in% names(aaindex)) {
    return(aaindex[names(aaindex) == x][[1]][["D"]] )
  } else if (x %in% dragonNames) {
    return(dragon$Description[dragonNames == x])
  } else {
    return(NA)
  }
})

write.csv(allFeatDF, "results/SR1vsPCP/features/allTypes_featureSim_norm_median.csv", row.names = F)

# ----- feature embedding -----

allProps = cbind(AADescAll%>% filter (rownames(AADescAllnorm) != "I") %>% arrange(aa) %>% select(-aa),
                 AAindex %>% arrange(aa))
features = names(allProps)[-which(names(allProps) == "aa")]

DB =  Qual
DB = left_join(DB, extract_aminoacids(tbl = DB, onlyValidSeq = T))

# get the P-1
pos = str_split_fixed(DB$positions,"_",Inf)[,c(1:4)]
DB$`P-1` = substr(DB$substrateSeq, as.numeric(pos[,2])+1, as.numeric(pos[,2])+1)

DBlarge = DB %>%
  select(substrateID, pepSeq, spliceType, positions, sr1, sr2, P4,P3,P2,P1,`P-1`) %>%
  tidyr::gather(residue, aa, -substrateID, -pepSeq, -spliceType, -positions, -sr1, -sr2)

# get amino acid features
DBfeat = left_join(DBlarge, allProps)
# replace missing amino acids by 0
DBfeat[is.na(DBfeat)] = 0

# group database into substrateID and spliceType
types = c("cis","revCis","trans")

SR1 = DBfeat %>%
  filter(spliceType %in% types) %>%
  filter(residue %in% c("P4","P3","P2","P1","P-1")) %>%
  # filter(residue %in% c("P1","P-1")) %>%
  distinct(substrateID, spliceType, sr1, residue, .keep_all = T) %>%
  group_by(pepSeq, sr1) %>%
  mutate(valid = ifelse(length(residue) == 5,T,F)) %>%
  # mutate(valid = ifelse(length(residue) == 2,T,F)) %>%
  filter(valid) %>% ungroup() %>%
  select(-valid) %>%
  arrange(residue) %>%  # make sure that residues have the same order
  as.data.table()

PCP = DBfeat %>%
  filter(spliceType == "PCP") %>%
  filter(residue %in% c("P4","P3","P2","P1","P-1")) %>%
  # filter(residue %in% c("P1","P-1")) %>%
  distinct(substrateID, spliceType, sr1, residue, .keep_all = T) %>%
  group_by(pepSeq, sr1) %>%
  mutate(valid = ifelse(length(residue) == 5,T,F)) %>%
  # mutate(valid = ifelse(length(residue) == 2,T,F)) %>%
  filter(valid) %>% ungroup() %>%
  select(-valid) %>%
  arrange(residue) %>%  # make sure that residues have the same order
  as.data.table()


# PCA on each feature
{
pdf("results/SR1vsPCP/_embeddingAll.pdf", height = 20, width = 12)
par(mfrow = c(5,3))

sapply(features, function(f){
  print(f)
  
  # get current features
  sr1 = SR1 %>%
    rename(cntFeature = f) %>%
    group_by(substrateID,spliceType,sr1) %>%
    summarise(P4 = cntFeature[residue == "P4"],
              P3 = cntFeature[residue == "P3"],
              P2 = cntFeature[residue == "P2"],
              P1 = cntFeature[residue == "P1"],
              `P-1` = cntFeature[residue == "P-1"]) %>%
    mutate(co = ifelse(spliceType == "cis", add.alpha(plottingCols["cis"], .5), add.alpha(plottingCols["revCis"], .5)),
           co = ifelse(spliceType == "trans", add.alpha(plottingCols["trans"], .5), co))
  
  pcp = PCP %>%
    rename(cntFeature = f) %>%
    group_by(substrateID,spliceType,sr1) %>%
    summarise(P4 = cntFeature[residue == "P4"],
              P3 = cntFeature[residue == "P3"],
              P2 = cntFeature[residue == "P2"],
              P1 = cntFeature[residue == "P1"],
              `P-1` = cntFeature[residue == "P-1"]) %>%
    mutate(co = add.alpha(plottingCols["PCP"], .5))
  both = rbind(sr1,pcp)
  
  # pca
  pca1 = prcomp(both[, c("P4","P3","P2","P1","P-1")])
  # um = umap(both[, c("P4","P3","P2","P1","P-1")], verbose = T,
  #           init = "agspectral", metric = "cosine", n_neighbors = 10)
  
  plot(pca1$x[,c(1:2)], # um[,1], um[,2]
       pch = 16, cex = .8,
       col = both$co, main = f)
})

dev.off()
}

# ----- UMAP embedding of selected features -----
impFeatures = read.csv("results/SR1vsPCP/features/allTypes_featureSim_norm.csv", stringsAsFactors = F)
cisF = impFeatures$.id[impFeatures$cis <= quantile(impFeatures$cis, 0.05)]
revCisF = impFeatures$.id[impFeatures$revCis <= quantile(impFeatures$revCis, 0.05)]
transF = impFeatures$.id[impFeatures$trans <= quantile(impFeatures$trans, 0.05)]

feat = c(cisF, revCisF, transF) %>% unique()
feat


# UMAP
set.seed(42)
pdf("results/SR1vsPCP/_embeddingSelected_UMAP.pdf", height = 5, width = 5)
# par(mfrow = c(5,2))

sapply(feat, function(f){
  print(f)
  
  # get current features
  sr1 = SR1 %>%
    rename(cntFeature = f) %>%
    group_by(substrateID,spliceType,sr1) %>%
    summarise(P4 = cntFeature[residue == "P4"],
              P3 = cntFeature[residue == "P3"],
              P2 = cntFeature[residue == "P2"],
              P1 = cntFeature[residue == "P1"],
              `P-1` = cntFeature[residue == "P-1"]) %>%
    mutate(co = ifelse(spliceType == "cis", add.alpha(plottingCols["cis"], .5), add.alpha(plottingCols["revCis"], .5)),
           co = ifelse(spliceType == "trans", add.alpha(plottingCols["trans"], .5), co))
  
  pcp = PCP %>%
    rename(cntFeature = f) %>%
    group_by(substrateID,spliceType,sr1) %>%
    summarise(P4 = cntFeature[residue == "P4"],
              P3 = cntFeature[residue == "P3"],
              P2 = cntFeature[residue == "P2"],
              P1 = cntFeature[residue == "P1"],
              `P-1` = cntFeature[residue == "P-1"]) %>%
    mutate(co = add.alpha(plottingCols["PCP"], .5))
  both = rbind(sr1,pcp)
  
  # umap
  um = umap(both[, c("P4","P3","P2","P1","P-1")], verbose = T,
            init = "agspectral", metric = "cosine", n_neighbors = 50)
  
  plot(um[,1], um[,2],
       pch = 16, cex = .8,
       xlab = "UMAP 1", ylab = "UMAP 2",
       col = both$co, main = f, sub = impFeatures$descriptions[impFeatures$.id == f])
  
  # feature value distribution
  print(ggplot(both %>% tidyr::gather(residue, val, -substrateID, -spliceType, -sr1, -co),
         aes(x = residue, y = val, fill = factor(spliceType, levels = c("PCP","cis","revCis","trans")))) +
    geom_boxplot() +
    xlab("residue") +
    ylab("z-scaled feature value") +
    ggtitle(f, subtitle = impFeatures$descriptions[impFeatures$.id == f]) +
    scale_fill_manual("type",values = c(plottingCols[["PCP"]], plottingCols[["cis"]],
                                 plottingCols[["revCis"]], plottingCols[["trans"]])))
})

dev.off()

# ----- feature value distributions -----
featureDesc = read.csv("results/SR1vsPCP/features/allTypes_featureSim.csv", stringsAsFactors = F)

allProps = cbind(AADescAll %>% filter (rownames(AADescAll) != "I") %>% arrange(aa) %>% select(-aa),
                 AAindex %>% arrange(aa))
features = names(allProps)[-which(names(allProps) == "aa")]

DB =  Qual
DB = left_join(DB, extract_aminoacids(tbl = DB, onlyValidSeq = T))

# get the P-1
pos = str_split_fixed(DB$positions,"_",Inf)[,c(1:4)]
DB$`P-1` = substr(DB$substrateSeq, as.numeric(pos[,2])+1, as.numeric(pos[,2])+1)

DBlarge = DB %>%
  select(substrateID, pepSeq, spliceType, positions, sr1, sr2, P4,P3,P2,P1,`P-1`) %>%
  tidyr::gather(residue, aa, -substrateID, -pepSeq, -spliceType, -positions, -sr1, -sr2)

# get amino acid features
DBfeat = left_join(DBlarge, allProps)
# replace missing amino acids by 0
# DBfeat[is.na(DBfeat)] = 0

# group database into substrateID and spliceType
types = c("cis","revCis","trans")

SR1 = DBfeat %>%
  filter(spliceType %in% types) %>%
  filter(residue %in% c("P4","P3","P2","P1","P-1")) %>%
  distinct(substrateID, spliceType, sr1, residue, .keep_all = T) %>%
  group_by(pepSeq, sr1) %>%
  mutate(valid = ifelse(length(residue) == 5,T,F)) %>%
  filter(valid) %>% ungroup() %>%
  select(-valid) %>%
  arrange(residue) %>%  # make sure that residues have the same order
  as.data.table()

PCP = DBfeat %>%
  filter(spliceType == "PCP") %>%
  filter(residue %in% c("P4","P3","P2","P1","P-1")) %>%
  distinct(substrateID, spliceType, sr1, residue, .keep_all = T) %>%
  group_by(pepSeq, sr1) %>%
  mutate(valid = ifelse(length(residue) == 5,T,F)) %>%
  filter(valid) %>% ungroup() %>%
  select(-valid) %>%
  arrange(residue) %>%  # make sure that residues have the same order
  as.data.table()


p = lapply(features, function(f){
  print(f)
  
  # get current features
  sr1 = SR1 %>%
    rename(cntFeature = f) %>%
    group_by(substrateID,spliceType,sr1) %>%
    summarise(P4 = cntFeature[residue == "P4"],
              P3 = cntFeature[residue == "P3"],
              P2 = cntFeature[residue == "P2"],
              P1 = cntFeature[residue == "P1"],
              `P-1` = cntFeature[residue == "P-1"]) %>%
    mutate(co = ifelse(spliceType == "cis", add.alpha(plottingCols["cis"], .5), add.alpha(plottingCols["revCis"], .5)),
           co = ifelse(spliceType == "trans", add.alpha(plottingCols["trans"], .5), co))
  
  pcp = PCP %>%
    rename(cntFeature = f) %>%
    group_by(substrateID,spliceType,sr1) %>%
    summarise(P4 = cntFeature[residue == "P4"],
              P3 = cntFeature[residue == "P3"],
              P2 = cntFeature[residue == "P2"],
              P1 = cntFeature[residue == "P1"],
              `P-1` = cntFeature[residue == "P-1"]) %>%
    mutate(co = add.alpha(plottingCols["PCP"], .5))
  both = rbind(sr1,pcp) %>%
    tidyr::gather(residue, val, -substrateID, -spliceType, -sr1, -co)
  
  # feature value distribution
  pl = ggplot(both, aes(x = residue, y = val, fill = factor(spliceType, levels = c("PCP","cis","revCis","trans")))) +
          geom_boxplot() +
          xlab("residue") +
          ylab("feature value") +
          ggtitle(f, subtitle = featureDesc$descriptions[featureDesc$.id == f]) +
          scale_fill_manual("type",values = c(plottingCols[["PCP"]], plottingCols[["cis"]],
                                              plottingCols[["revCis"]], plottingCols[["trans"]]))
  
  return(pl)
})

ggsave(filename = "results/SR1vsPCP/_featureDistr.pdf", 
       plot = gridExtra::marrangeGrob(p, nrow=10, ncol=5), 
       width = 20, height = 40, dpi = "retina")


# select features that have a high difference between PCP and PSP
# ignore missing values
# SR1norm = apply(SR1[,features], 2, function(x){
#   (x - min(x, na.rm = T))/(max(x, na.rm = T) - min(x, na.rm = T))
# })
# PCPnorm = apply(PCP[,features], 2, function(x){
#   (x - min(x, na.rm = T))/(max(x, na.rm = T) - min(x, na.rm = T))
# })


# ----- relevant feature value distributions -----
manually = read.table("results/SR1vsPCP/features/distance_approach/features_manually")
featManually = manually$V1 %>% unique()

pMan = lapply(featManually, function(f){
  print(f)
  
  # get current features
  sr1 = SR1 %>%
    rename(cntFeature = f) %>%
    group_by(substrateID,spliceType,sr1) %>%
    summarise(P4 = cntFeature[residue == "P4"],
              P3 = cntFeature[residue == "P3"],
              P2 = cntFeature[residue == "P2"],
              P1 = cntFeature[residue == "P1"],
              `P-1` = cntFeature[residue == "P-1"]) %>%
    mutate(co = ifelse(spliceType == "cis", add.alpha(plottingCols["cis"], .5), add.alpha(plottingCols["revCis"], .5)),
           co = ifelse(spliceType == "trans", add.alpha(plottingCols["trans"], .5), co))
  
  pcp = PCP %>%
    rename(cntFeature = f) %>%
    group_by(substrateID,spliceType,sr1) %>%
    summarise(P4 = cntFeature[residue == "P4"],
              P3 = cntFeature[residue == "P3"],
              P2 = cntFeature[residue == "P2"],
              P1 = cntFeature[residue == "P1"],
              `P-1` = cntFeature[residue == "P-1"]) %>%
    mutate(co = add.alpha(plottingCols["PCP"], .5))
  both = rbind(sr1,pcp) %>%
    tidyr::gather(residue, val, -substrateID, -spliceType, -sr1, -co)
  
  # feature value distribution
  pl = ggplot(both, aes(x = residue, y = val, fill = factor(spliceType, levels = c("PCP","cis","revCis","trans")))) +
    geom_boxplot() +
    xlab("residue") +
    ylab("feature value") +
    ggtitle(f, subtitle = featureDesc$descriptions[featureDesc$.id == f]) +
    scale_fill_manual("type",values = c(plottingCols[["PCP"]], plottingCols[["cis"]],
                                        plottingCols[["revCis"]], plottingCols[["trans"]]))
  
  return(pl)
})

ggsave(filename = "results/SR1vsPCP/_featureDistrManually.pdf", 
       plot = gridExtra::marrangeGrob(pMan, nrow=5, ncol=5, byrow = T), 
       width = 20, height = 20, dpi = "retina")





