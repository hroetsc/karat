### INHIBITOR KINETICS ###
# description:  select peptides that are downregulated when inhibiting b5
# input:        limma output (DEanalysis_limma.R)
# output:       b5 specific peptides
# author:       HR

library(dplyr)
library(stringr)
library(eulerr)


### INPUT ###
fs = list.files("results/b2/", pattern = "DEanalysis_", recursive = T, full.names = T)
fs_int = fs[grepl("_intensities.csv",fs)]
fs_folds = fs[grepl("fold-changes.csv",fs)]


### MAIN PART ###
# ----- extract downregulated peptides from all comparisons -----
intensities = list()
fold_changes = list()

for (i in 1:length(fs_int)) {
  
  intensities[[i]] = read.csv(fs_int[i], stringsAsFactors = F) %>%
    as.data.frame()
  names(intensities)[i] = str_extract(basename(fs_int[i]), "b[:digit:]_vs_[:alnum:]+")
  
  fold_changes[[i]] = read.csv(fs_folds[i], stringsAsFactors = F) %>%
    as.data.frame()
  names(fold_changes)[i] = str_extract(basename(fs_folds[i]), "b[:digit:]_vs_[:alnum:]+")
  
}

intensities = plyr::ldply(intensities)
fold_changes = plyr::ldply(fold_changes)


# ----- determine overlap -----
ovl = list("all other" = intensities$pepSeq[intensities$.id == "b2_vs_all"],
           "ß1" = intensities$pepSeq[intensities$.id == "b2_vs_b1"],
           "ß5" = intensities$pepSeq[intensities$.id == "b2_vs_b5"],
           "no inhibitor" = intensities$pepSeq[intensities$.id == "b2_vs_no"]) %>%
  euler(shape="ellipse")

set.seed(123)
plot(ovl, quantities = T, main = "downregulated upon ß2 inhibition compared to ...")

# png("results/b2/overlapBetweenComp.png", height = 8, width = 10, units = "in", res = 300)
# print(plot(ovl, quantities = T, main = "b2-specific peptides"))
# dev.off()

png("~/Documents/Studium/Fachvertiefung+BA/Praktikumsbericht/plots/overlap_b2.png",
    height = 6, width = 8, units = "in", res = 300)
print(plot(ovl, quantities = T, main = "downregulated upon ß2 inhibition compared to ..."))
dev.off()


# take all peptides and evaluate manually?

