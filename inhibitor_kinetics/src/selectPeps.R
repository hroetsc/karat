### INHIBITOR KINETICS ###
# description:  select peptides that are downregulated when inhibiting b5
# input:        limma output (DEanalysis_limma.R)
# output:       b5 specific peptides
# author:       HR

library(dplyr)
library(stringr)
library(eulerr)


### INPUT ###
fs_int = list.files("results/b5/", pattern = "intensities.csv", recursive = T, full.names = T)
fs_folds = list.files("results/b5/", pattern = "fold-changes.csv", recursive = T, full.names = T)


### MAIN PART ###
# ----- extract downregulated peptides from all comparisons -----
intensities = list()
fold_changes = list()

for (i in 1:length(fs_int)) {
  
  intensities[[i]] = read.csv(fs_int[i], stringsAsFactors = F)
  names(intensities)[i] = str_extract(basename(fs_int[i]), "b[:digit:]_vs_[:alnum:]+")
  
  fold_changes[[i]] = read.csv(fs_folds[i], stringsAsFactors = F)
  names(fold_changes)[i] = str_extract(basename(fs_folds[i]), "b[:digit:]_vs_[:alnum:]+")
  
}

intensities = plyr::ldply(intensities) %>%
  na.omit()
fold_changes = plyr::ldply(fold_changes) %>%
  na.omit()


# ----- determine overlap -----
ovl = list(all = intensities$pepSeq[intensities$.id == "b5_vs_all"],
           b1 = intensities$pepSeq[intensities$.id == "b5_vs_b1"],
           b2 = intensities$pepSeq[intensities$.id == "b5_vs_b2"],
           no = intensities$pepSeq[intensities$.id == "b5_vs_no"]) %>%
  euler(shape="ellipse")

set.seed(123)
plot(ovl, quantities = T, main = "b5-specific peptides")

png("results/b5/overlapBetweenComp.png", height = 8, width = 10, units = "in", res = 300)
print(plot(ovl, quantities = T, main = "b5-specific peptides"))
dev.off()

# take all peptides and evaluate manually?

# ----- select peptides ------
# were detected in at least 2 comparisons
pepU = unique(intensities$pepSeq)
counts = rep(NA, length(pepU))

for (j in 1:length(pepU)) {
  counts[j] = length(which(intensities$pepSeq == pepU[j]))
}

pepHigh = pepU[counts >= 2]
length(pepHigh)

# select peptides that were significantly downregulated in at least 2 comparisons
intensitiesU = intensities[intensities$pepSeq %in% pepHigh, names(intensities) != ".id"] %>%
  unique()
fold_changesU = fold_changes[fold_changes$pepSeq %in% pepHigh, names(fold_changes) != ".id"] %>%
  unique()

# some basic statistics
table(intensitiesU$spliceType)


### OUTPUT ###
write.csv(intensitiesU, file = "results/b5/selectedPeps_intensities.csv",
          row.names = F)
write.csv(fold_changesU, file = "results/b5/selectedPeps_fold-changes.csv",
          row.names = F)

