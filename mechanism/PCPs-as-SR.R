### karat projetc - PCPS mechanism ###
# description:  investigate whether PCPs are more often used as SR1s or SR2s
# input:        all data sets
# output:       stats reg usage of PCPs as SRs
# author:       HPR

library(dplyr)
library(stringr)
library(ggplot2)
source("src/invitroSPI_utils.R")

theme_set(theme_classic())

### INPUT ###
load("data/invitroSPI.RData")
load("../../proteinsPCPS/new/data/aSPIre.RData")
proteins = Kinetics
load("data/aSPIre.RData")
polypeps = Kinetics


### MAIN PART ###
suppressWarnings(dir.create("results/SR-as-PCP/"))

# ----- preprocessing -----
nm = intersect(names(ProteasomeDB), intersect(names(proteins), names(polypeps)))
ALL = rbind(ProteasomeDB[,nm] %>% mutate(seqtype = "polypeptides"),
            proteins[,nm] %>% mutate(seqtype = "proteins"),
            polypeps[,nm] %>% mutate(seqtype = "polypeptides")) %>%
  ILredundancy() %>%
  uniquePeptides() %>%
  tidyr::separate_rows(positions, sep = ";")

# ----- extract SRs -----
pos = str_split_fixed(ALL$positions, "_", Inf)
pos = apply(pos,2,as.numeric)

ALL$sr1 = str_sub(ALL$substrateSeq, pos[,1], pos[,2])
ALL$sr2 = str_sub(ALL$substrateSeq, pos[,3], pos[,4])
ALL$P1 = pos[,2]
ALL$P1_ = pos[,3]
ALL$P1_[ALL$spliceType == "PCP"] = pos[ALL$spliceType == "PCP",1]


# ----- exact matches: PCPs detected as SR -----
subIDs = ALL$substrateID %>% unique()
ALL$SrAsPCP = "none"

ALLinfo = lapply(subIDs, function(s){
  
  PSPdb = ALL[ALL$substrateID == s & ALL$spliceType != "PCP", ]
  PCPdb = ALL[ALL$substrateID == s & ALL$spliceType == "PCP", ]
  pcps = PCPdb$pepSeq %>% unique()
  
  k1 = which(PSPdb$sr1 %in% pcps)
  k2 = which(PSPdb$sr2 %in% pcps)
  kb = intersect(k1,k2)
  
  if (length(k1) > 0) {
    PSPdb$SrAsPCP[k1] = "SR1"
  }
  if (length(k1) > 0) {
    PSPdb$SrAsPCP[k2] = "SR2"
  }
  if (length(kb) > 0) {
    PSPdb$SrAsPCP[kb] = "both"
  }
  
  return(PSPdb)
})


ALLinfo = plyr::ldply(ALLinfo)

# stats and plotting

infoinfo = ALLinfo %>%
  mutate(spliceType = ifelse(spliceType %in% c("cis","revCis","type_multi-mapper"), "allcis", spliceType)) %>%
  group_by(seqtype, substrateID, spliceType, SrAsPCP) %>%
  summarise(n = n()) %>%
  tidyr::spread(SrAsPCP, n, fill = 0) %>%
  tidyr::gather(SrAsPCP, n, -seqtype, -substrateID, -spliceType) %>%
  ungroup() %>%
  group_by(seqtype, substrateID, spliceType) %>%
  mutate(freq = n/sum(n))

ggplot(infoinfo, aes(x = SrAsPCP, y = freq, fill = spliceType)) +
  # geom_violin(draw_quantiles = 0.5) +
  geom_boxplot() +
  scale_fill_manual(values = c(plottingCols[["allcis"]], plottingCols[["trans"]])) +
  facet_wrap(~seqtype, scales = "free") +
  ggtitle("PCP usage as SR for PSPs",
          subtitle = "considering exact sequences") +
  ylab("frequency per substrate") + xlab("usage as ...")
ggsave(filename = "results/SR-as-PCP/exactSequences.png", plot = last_plot(), height = 5, width = 9, dpi = "retina")


# ----- get precursors matches -----

ALL$SR1used = NA
ALL$SR2used = NA
ALL$bothused = NA

ALLprecinfo = lapply(subIDs, function(s){
  
  PSPdb = ALL[ALL$substrateID == s & ALL$spliceType != "PCP", ]
  PCPdb = ALL[ALL$substrateID == s & ALL$spliceType == "PCP", ]
  
  k1 = sapply(1:nrow(PSPdb), function(j){
    length(which(grepl(pattern = PSPdb$sr1[j], x = PCPdb$pepSeq) & PCPdb$P1 == PSPdb$P1[j]))
  })
  PSPdb$SR1used = k1
  
  k2 = sapply(1:nrow(PSPdb), function(j){
    length(which(grepl(pattern = PSPdb$sr2[j], x = PCPdb$pepSeq) & PCPdb$P1_ == PSPdb$P1_[j]))
  })
  PSPdb$SR2used = k2
  
  kb = sapply(1:nrow(PSPdb), function(j){
    length(which(grepl(pattern = PSPdb$sr1[j], x = PCPdb$pepSeq) & PCPdb$P1 == PSPdb$P1[j] & grepl(pattern = PSPdb$sr2[j], x = PCPdb$pepSeq) & PCPdb$P1_ == PSPdb$P1_[j]))
  })
  PSPdb$bothused = kb
  
  return(PSPdb)
})


ALLprecinfo = plyr::ldply(ALLprecinfo)

# stats and plotting

infoinfoinfo = ALLprecinfo %>%
  mutate(SrAsPCP = ifelse(SR1used > 0, "SR1", "none"),
         SrAsPCP = ifelse(SR2used > 0, "SR2", SrAsPCP),
         SrAsPCP = ifelse(SR1used > 0 & SR2used > 0, "both", SrAsPCP)) %>%
  select(-SR1used, -SR2used, -bothused) %>%
  mutate(spliceType = ifelse(spliceType %in% c("cis","revCis","type_multi-mapper"), "allcis", spliceType)) %>%
  group_by(seqtype, substrateID, spliceType, SrAsPCP) %>%
  summarise(n = n()) %>%
  tidyr::spread(SrAsPCP, n, fill = 0) %>%
  tidyr::gather(SrAsPCP, n, -seqtype, -substrateID, -spliceType) %>%
  ungroup() %>%
  group_by(seqtype, substrateID, spliceType) %>%
  mutate(freq = n/sum(n))

ggplot(infoinfoinfo, aes(x = SrAsPCP, y = freq, fill = spliceType)) +
  # geom_violin(draw_quantiles = 0.5) +
  geom_boxplot() +
  scale_fill_manual(values = c(plottingCols[["allcis"]], plottingCols[["trans"]])) +
  facet_wrap(~seqtype, scales = "free") +
  ggtitle("PCP usage as SR precursor for PSPs",
          subtitle = "considering precursors") +
  ylab("frequency per substrate") + xlab("usage as precursor for...")
ggsave(filename = "results/SR-as-PCP/precursorSequences.png", plot = last_plot(), height = 5, width = 9, dpi = "retina")


# join all info
MASTERinfo = inner_join(ALLinfo, ALLprecinfo %>% select(-SrAsPCP))
save(MASTERinfo, file = "results/SR-as-PCP/MASTERinfo.RDAta")
