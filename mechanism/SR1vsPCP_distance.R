### karat projetc - PCPS mechanism ###
# description:  determine differences between SR1 and PCP (precursors)
#               qualitative data set: Roetschke et al. SciData, EGFR, WT sequences of WT/Mut
# output:       distance between SR1s and PCPs for different features
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
library(eulerr)
library(dbscan)
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
suppressWarnings(dir.create("results/SR1vsPCP/features/distance_approach/"))


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
getFeatureDist = function(PSPdbF, PCPdbF, sr, pcp) {
  
  if (length(sr) > 0 & length(pcp) > 0) {
    # get row indices that belong to each position
    # (precompute to be faster)
    
    srU = unique(sr) %>% sort()
    sridx = lapply(srU, function(k){
      which(sr == k)
    })
    
    pcpU = unique(pcp) %>% sort()
    pcpidx = lapply(pcpU, function(k){
      which(pcp == k)
    })
    
    # calculate pairwise feature similarities for each SR
    if (length(srU) == length(pcpU)) {
      xx = sapply(seq(1,length(srU)), function(i){
        a = na.omit(PSPdbF[sridx[i][[1]]])
        b = na.omit(PCPdbF[pcpidx[i][[1]]])
        if(!!length(a) & !!length(b)) {
          wasserstein1d(a,b)
          # ad_test(a,b, nboots = 500)[2]
        }
      })
    } else {
      print("something seems to be wrong")
      xx = NA
    }
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
  
  # for each feature: get distance between feature distributions of each position
  # iterate features
  pdf(paste0("results/SR1vsPCP/features/distance_approach/",cntDS,"_norm.pdf"), height = 40, width = 15)
  par(mfrow = c(10,length(types)))
  
  meanVars = lapply(features, function(f){
    print(f)
    
    meanVarTypes = lapply(ALL, function(cntType){
      
      meanVarSubs = lapply(cntType, function(cntSubs){
        
        # print(cntSubs$PCP$substrateID[1])
        psp = unlist(cntSubs$PSP[,..f])
        pcp = unlist(cntSubs$PCP[,..f])
        
        xx = getFeatureDist(PSPdbF = psp,
                            PCPdbF = pcp,
                            sr = cntSubs$PSP$residue,
                            pcp = cntSubs$PCP$residue)
        
        return(data.table(t(xx)))
      })
      
      return(meanVarSubs %>% rbindlist(fill = T))
    })
    
    xxx = plyr::ldply(meanVarTypes)
    names(xxx) = c("spliceType", sort(c("P4","P3","P2","P1","P-1")))
    xxx = xxx %>%
      tidyr::gather(residue, Wdist, -spliceType) %>% mutate(
        Wdist = sapply(Wdist, function(d){
          if(is.null(d)){
            return(NA)
          } else {
            return(d)
          }
        }) %>% as.numeric()
      ) %>%
      na.omit()
    
    sapply(types, function(t){
      boxplot(data = xxx %>% filter(spliceType == t), Wdist~residue, main = f, sub = t)
    })
    
    return(xxx %>% group_by(spliceType, residue) %>% summarise(median = median(Wdist)))
    
  })
  
  dev.off()
  names(meanVars) = features
  save(meanVars, file = paste0("results/SR1vsPCP/features/distance_approach/",cntDS,"_norm.RData"))
  
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

getFeatureSimilarity(descriptor = AADescAll,
                     features = names(AADescAll)[names(AADescAll) != "aa"],
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


# ---- load feature distances -----

allDS = c("AADescAll","AAindex")

allFeat = list()
counter = 1
for (cntDS in allDS) {
  print(cntDS)
  
  load(paste0("results/SR1vsPCP/features/distance_approach/",cntDS,"_norm.RData"))
  meanVarsDF = meanVars %>%
    plyr::ldply()
  
  cnt = ggplot(meanVarsDF, aes(x = residue, y = median+1, fill = spliceType)) +
    geom_violin(draw_quantiles = c(0.25,0.5,0.75)) +
    scale_fill_manual(values = c(plottingCols[["cis"]], plottingCols[["revCis"]], plottingCols[["trans"]])) +
    ylab("median Wassterstein distance between PCP/PSP features") +
    ggtitle("PCP/PSP precursor feature distance", subtitle = paste0(cntDS, ", n = ", length(meanVars)))
  
  ggsave(filename = paste0("results/SR1vsPCP/features/distance_approach/_distr",cntDS,"_norm.png"),
         plot = cnt, dpi = "retina", height = 6, width = 6)
  
  meanVarsDF %>%
    group_by(spliceType, residue) %>%
    summarise(n = n(), mean = mean(median, na.rm = T),
              median = median(median, na.rm = T),
              std = sd(median, na.rm = T)) %>%
    print.data.frame()
  
  allFeat[[counter]] = meanVarsDF
  counter = counter+1
}

allFeatDF = plyr::ldply(allFeat) %>%
  tidyr::spread(residue, median)

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

write.csv(allFeatDF, "results/SR1vsPCP/features/distance_approach/allTypes_featureSim_norm.csv", row.names = F)


# ----- identify relevant features -----

residues = c("P4","P3","P2","P1","P-1")
types = c("cis","revCis","trans")

pdf("results/SR1vsPCP/features/distance_approach/_residue-wise.pdf", height = 16, width = 12)
par(mfrow = c(4,3))
for (r in residues[residues != "P1"]) {
  for (t in types) {
    
    plot(x = allFeatDF$P1[allFeatDF$spliceType == t],
         y = allFeatDF[allFeatDF$spliceType == t, r],
         col = plottingCols[t],
         xlab = "distance @ P1",
         ylab = paste0("distance @ ", r),
         pch = 16, main = t)
    abline(a=0,b=1)
    
  }
}
dev.off()

interesting = lapply(types, function(t){
  cnt = allFeatDF[allFeatDF$spliceType == t, ]
  res = lapply(residues, function(x){
    cnt$.id[cnt[,x] >= quantile(cnt[,x], 0.95)]
  })
  names(res) = residues
  return(res)
})
names(interesting) = types


euler(interesting$cis, shape = "ellipse") %>%
  plot(quantities = T, main = "cis")

euler(interesting$revCis, shape = "ellipse") %>%
  plot(quantities = T, main = "revCis")

euler(interesting$trans, shape = "ellipse") %>%
  plot(quantities = T, main = "trans")


# ----- HDBSCAN clustering of relevant features -----
allProps = cbind(AADescAllnorm %>% filter (rownames(AADescAllnorm) != "I") %>% arrange(aa),
                 AAindexnorm %>% arrange(aa))


pdf("results/SR1vsPCP/features/distance_approach/_HDBSCAN_residue-wise.pdf", height = 16, width = 12)
par(mfrow = c(5,3))
for (r in residues) {
  for (t in types) {
    
    cnt = interesting[names(interesting) == t][[1]]
    cl = hdbscan(allProps[, colnames(allProps) %in% cnt[[r]]], minPts = 2)
    plot(cl$hc, col = plottingCols[t],
         main = r, sub = t,
         xlab="HDBSCAN* Hierarchy", labels = gsub("L","I/L",allProps$aa))
    
  }
}
dev.off()


# cluster partition energy features
PartEnFeat = c("GUYH850105",
         "MIYS850101",
         "GUYH850104",
         "GUYH850103",
         "GUYH850101",
         "GUYH850102",
         "ZASB820101")
cl = hdbscan(allProps[, colnames(allProps) %in% PartEnFeat], minPts = 2)
plot(cl$hc,
     main = "partition energy",
     xlab="HDBSCAN* Hierarchy", labels = gsub("L","I/L",allProps$aa))

# cluster charge features
ChargeFeat = c("FAUJ880112", "GGI6", "JGI6", "FAUJ880111", "CHAM830107",
               "JGI7", "GGI7", "CHAM830108", "GGI5", "JGI5",
               "GGI4", "JGI4", "GGI1", "GGI3", "JGI2", "JGI3", "KLEP840101",
               "JGI1", "GGI2", "JGT")
cl = hdbscan(allProps[, colnames(allProps) %in% ChargeFeat], minPts = 2)
plot(cl$hc,
     main = "charge",
     xlab="HDBSCAN* Hierarchy", labels = gsub("L","I/L",allProps$aa))

# cluster polarisability features
polarFeat = c(
  "RDF035p",
  "Mor17e",
  "RADA880108",
  "PLIV810101",
  "EISD860103",
  "Me",
  "CHAM820101",
  "Mp",
  "Se",
  "nHDon",
  "Sp",
  "ROBB790101",
  "JANJ790102"
)

sapply(str_replace_all(polarFeat, "[:punct:]", "_"), function(x){
  if(x %in% names(aaindex)) {
    return(aaindex[names(aaindex) == x][[1]][["D"]] )
  } else if (x %in% dragonNames) {
    return(dragon$Description[dragonNames == x])
  } else {
    return(NA)
  }
})


cl = hdbscan(allProps[, colnames(allProps) %in% polarFeat], minPts = 2)
png("results/SR1vsPCP/polarisabilty.png", height = 5, width = 7, res = 300, units = "in")
plot(cl$hc,
     main = "polarisability",
     sub = paste(polarFeat, collapse = ", "),
     xlab="HDBSCAN* Hierarchy", labels = gsub("L","I/L",allProps$aa))
dev.off()


# cluster PCP/PSP
curiousP1F = c(
  "Hy",  # hydrophilic factor
  "nH",
  "ZIMJ680102", # bulkiness
  "BHAR880101",  # flexibility
  "WOLR810101",  # hydration potential
  "ROBB790101",  # hydration free energy
  "Mp",  # mean atomic polarisability
  "Me",  # mean electronegativity
  "CHAM820101"  # polarisability parameter
)

cl = hdbscan(allProps[, colnames(allProps) %in% curiousP1F], minPts = 2)
png("results/SR1vsPCP/curiousP1.png", height = 6, width = 10, res = 300, units = "in")
plot(cl$hc,
     main = "P1 features",
     sub = paste(curiousP1F, collapse = ", "),
     xlab="HDBSCAN* Hierarchy", labels = gsub("L","I/L",allProps$aa))
dev.off()


curiousP1p = c(
  "nH",
  "nHDon",
  "BHAR880101",  # flexibility
  "EISD860103",  # direction of hydrophobic moment
  "nN"
)

cl = hdbscan(allProps[, colnames(allProps) %in% curiousP1p], minPts = 2)
png("results/SR1vsPCP/curiousP1_.png", height = 6, width = 10, res = 300, units = "in")
plot(cl$hc,
     main = "P1' features",
     sub = paste(curiousP1p, collapse = ", "),
     xlab="HDBSCAN* Hierarchy", labels = gsub("L","I/L",allProps$aa))
dev.off()


