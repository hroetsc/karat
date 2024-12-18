### karat projetc - PCPS mechanism ###
# description:  get intensities / fold-changes for all WT-Mut pairs in all replicates
# input:        WT-Mut pairs, aSPIre quantification results
# output:       fold-changes for all pairs
# author:       HPR

library(dplyr)
library(stringr)

### INPUT ###
WTMut = read.csv("results/WTMut/WT-Mut-pairs.csv", stringsAsFactors = F)

mutfs = list.files("invitroSPI+aSPIre/results/", pattern = "QUANTITIES_raw.RData",
                   recursive = T, full.names = T)
wtfs = list.files("aSPIre_Mut/results/", pattern = "QUANTITIES_raw.RData",
                   recursive = T, full.names = T)

### MAIN PART ###
suppressWarnings(dir.create("results/WTMut/fold-changes/"))

# ----- preprocessing -----
# load quantity tables
mut = lapply(mutfs, function(x){
  load(x)
  return(QUANTITIES %>% select(substrateID,substrateSeq,pepSeq,biological_replicate,digestTime,intensity))
})
mut = plyr::ldply(mut)

wt = lapply(wtfs, function(x){
  load(x)
  return(QUANTITIES %>% select(substrateID,substrateSeq,pepSeq,biological_replicate,digestTime,intensity))
})
wt = plyr::ldply(wt)

Kinetics = rbind(wt, mut)

# merge with WT/Mut pairs
MASTER = inner_join(Kinetics, WTMut) %>%
  filter(!is.na(residue)) %>%
  select(substrateID, substrateSeq, pepSeq, biological_replicate, digestTime, intensity, protein_name, origSubs.1, spliceType, positions, residue, identical)


# ----- split into comparison groups -----

Y = MASTER %>%
  rename(origSubs = origSubs.1) %>%
  group_by(substrateID, pepSeq, positions, biological_replicate, digestTime) %>%
  mutate(technical_replicate = row_number(),
         condition = paste0(protein_name,"-bio",biological_replicate,"-tech",technical_replicate),
         pepPos = paste0(origSubs,"-",positions)) %>%
  ungroup()

# get number of peptides per comparison
Y = Y %>%
  group_by(origSubs, positions) %>%
  mutate(numberPeps = length(unique(protein_name))) %>%
  ungroup()

# make sure that there are indeed only pairs compared to each other
Y = Y %>%
  filter(numberPeps > 1)

# normalise kinetics
Y = Y %>%
  group_by(substrateID, pepSeq, biological_replicate, technical_replicate) %>%
  mutate(intensity = as.numeric(intensity)) %>%
  mutate(intensity = intensity-min(intensity)) %>%
  ungroup()
Y$intensity[Y$digestTime == 0] = 0


# ----- get fold changes -----

origSubs = Y$origSubs %>% unique()

pdf("results/WTMut/fold-changes/stats.pdf", height = 8, width = 15)
par(mfrow = c(2,3))

for (o in origSubs) {
  print(o)
  
  I0 = Y %>%
    filter(origSubs == o & digestTime == 1) %>%
    select(pepPos, condition, intensity) %>%
    tidyr::spread(condition, -pepPos) %>%
    tibble::column_to_rownames("pepPos") %>%
    as.matrix()
  
  I4 = Y %>%
    filter(origSubs == o & digestTime == 4) %>%
    select(pepPos, condition, intensity) %>%
    tidyr::spread(condition, -pepPos) %>%
    tibble::column_to_rownames("pepPos") %>%
    as.matrix()
  
  # filter out peptides that have a lower intensity at 4 than at 0 hours
  # rem = sapply(1:nrow(I0), function(i){
  #   any(I0[i,] < I4[i,])
  # }) %>% which()
  # print(length(rem) / nrow(I4))
  # 
  # I0 = I0[rem,]
  # I4 = I4[rem,]
  
  # set 0 intensities to NA
  I0[I0==0] = NA
  I4[I4==0] = NA
  
  # background correction with 5% quantile
  cnames = str_split_fixed(str_split_fixed(colnames(I0), "-", Inf)[,1], "_", Inf)[,2]
  q = sapply(unique(cnames), function(c){
    quantile(rbind(I0[,cnames == c], I4[,cnames == c]), 0.05, na.rm = T, names = F)
  })

  for(c in unique(cnames)) {
    I0[,cnames == c] = I0[,cnames == c]/q[c]
    I4[,cnames == c] = I4[,cnames == c]/q[c]
  }
  
  # impute missing values
  I0[!is.finite(I0)] = 0
  I0 = apply(I0, 2, function(x){
    x = x+min(x[x!=0])
  })
  
  I4[!is.finite(I4)] = 0
  I4 = apply(I4, 2, function(x){
    x = x+min(x[x!=0])
  })
  
  # fold change + log-transformation
  FC = (log10(I4)-log10(I0))
  D4 = log10(I4)

  D4norm = apply(D4,2,function(x){
    return((x - min(x))/(max(x)-min(x)))
  })
  # density(D4norm) %>% plot(main = o, sub = "distribution of log fold-changes")
  
  # tmp!
  # FC = D4norm
  
  density(FC) %>% plot(main = o, sub = "distribution of log fold-changes")
  save(FC, file = paste0("results/WTMut/fold-changes/",o,".RData"))
  
}
dev.off()

### OUTPUT ###
write.csv(Y, "results/WTMut/WT-Mut-pairs+intensities.csv", row.names = F)

