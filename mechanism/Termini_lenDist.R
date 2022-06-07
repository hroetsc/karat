### karat projetc - PCPS mechanism ###
# description:  deduce the minimal splice-reactant length required for the proteasome
# input:        qualitative data set: Roetschke et al. SciData, EGFR, WT sequences of WT/Mut
#               + all theoretically possible peptides
# output:       length distributions of cis-peptides carrying the substrate's termini
# author:       HPR

library(dplyr)
library(stringr)
library(ggplot2)
source("src/plotting_utils.R")
source("src/invitroSPI_utils.R")
source("../brainstorming/src/number-of-products.R")
source("src/aSPIre_plotting.R")

theme_set(theme_classic())

### INPUT ###
load("../../proteinsPCPS/new/data/aSPIre.RData")
proteins = Kinetics
load("data/invitroSPI.RData")
# has to be loaded in this order!!!
load("data/aSPIre.RData")
polypepsQuant = Kinetics

### MAIN PART ###
# ----- data preprocessing -----

polypeps = ProteasomeDB %>%
  ILredundancy() %>%
  filterPepLength() %>%
  disentangleMultimappers.Type() %>%
  disentangleMultimappers.SRlen() %>%
  disentangleMultimappers.IVSeqLen() %>%
  removeMultimappers.SRlen() %>%
  removeMultimappers.IVSeqLen() %>%
  remSynthErrors() %>%
  filter20Sstandard() %>%
  DBcosmetics()

proteins = proteins %>%
  ILredundancy() %>%
  disentangleMultimappers.Type() %>%
  disentangleMultimappers.SRlen() %>%
  disentangleMultimappers.IVSeqLen() %>%
  removeMultimappers.SRlen() %>%
  removeMultimappers.IVSeqLen() %>%
  distinct(substrateID, pepSeq, .keep_all = T)

polypepsQuant = polypepsQuant %>%
  ILredundancy() %>%
  disentangleMultimappers.Type() %>%
  disentangleMultimappers.SRlen() %>%
  disentangleMultimappers.IVSeqLen() %>%
  removeMultimappers.SRlen() %>%
  removeMultimappers.IVSeqLen() %>%
  distinct(substrateID, pepSeq, .keep_all = T)

# ----- fetch cis-spliced carrying the substrate's termini -----

fetchSRs = function(df) {
  
  df_filter = df %>%
    filter(spliceType == "cis") %>%
    mutate(L = nchar(substrateSeq),
           N = nchar(pepSeq)) %>%
    removeMultimappers.Type() %>%
    uniquePeptides()
  
  df_filter$pepSeq %>% unique() %>% length() %>% print()
  
  pos = str_split_fixed(df_filter$positions,coll("_"),Inf)[,c(1:4)]
  pos = apply(pos,2,as.numeric) %>% as.data.frame()
  names(pos) = c("pos1","pos2","pos3","pos4")
  
  MASTER = cbind(df_filter, pos) %>%
    na.omit() %>%
    mutate(contains_term = ifelse(pos4 == L, "C", NA),
           contains_term = ifelse(pos1 == 1, "N", contains_term)) %>%
    filter(contains_term == "N" | contains_term == "C") %>%
    mutate(sr1len = ifelse(contains_term == "N", pos2-pos1+1, NA),
           sr2len = ifelse(contains_term == "C", pos4-pos3+1, NA)) %>%
    select(contains_term,sr1len,sr2len) %>%
    tidyr::gather(sr,len,-contains_term) %>%
    na.omit() %>%
    group_by(contains_term, sr, len) %>%
    summarise(n = n())
  
  return(MASTER)
}

polypeps_SR = fetchSRs(polypeps)
Qpolypeps_SR = fetchSRs(polypepsQuant)
prots_SR = fetchSRs(proteins)

# ----- get length distribution of splice-reactants at termini -----

po = ggplot(polypeps_SR, aes(x=len, y = n, fill = sr)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_continuous(breaks = seq(1,30,2)) +
  scale_fill_manual("", values = c("gray","lightblue"), labels = c("SR1 (N-terminus)","SR2 (C-terminus)")) +
  xlab("splice-reactant length (aa residues)") +
  ylab("counts") +
  ggtitle("polypeptides")

pr = ggplot(prots_SR, aes(x=len, y = n, fill = sr)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_continuous(limits = c(1,30), breaks = seq(1,30,2)) +
  scale_fill_manual("", values = c("gray","lightblue"), labels = c("SR1 (N-terminus)","SR2 (C-terminus)")) +
  xlab("splice-reactant length (aa residues)") +
  ylab("counts") +
  ggtitle("proteins")

ggsave(filename = "results/termini/SRlenAtTerm.png", plot = po, height = 5, width = 5, dpi = "retina")
ggsave(filename = "results/termini/SRlenAtTerm_proteins.png", plot = pr, height = 5, width = 5, dpi = "retina")
# "length distribution of cis splice-reactants", subtitle = "carrying the substrate's N-terminus (SR1) or C-terminus (SR2)"


# ----- investigate kinetics of peptides with SR2 @ C-term -----
suppressWarnings(dir.create("results/termini/singleAASR2/"))

poq = ggplot(Qpolypeps_SR, aes(x=len, y = n, fill = sr)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_continuous(breaks = seq(1,30,2)) +
  scale_fill_manual("", values = c("gray","lightblue"), labels = c("SR1 (N-terminus)","SR2 (C-terminus)")) +
  xlab("splice-reactant length (aa residues)") +
  ylab("counts") +
  ggtitle("polypeptides - quantitative data set")
ggsave(filename = "results/termini/singleAASR2/SRlenAtTerm_quantPolypeps.png", plot = poq, height = 5, width = 5, dpi = "retina")


# determine kinetics of peptides with 1/2 aa SR2
df_filter = Kinetics %>%
  ILredundancy() %>%
  disentangleMultimappers.Type() %>%
  disentangleMultimappers.SRlen() %>%
  disentangleMultimappers.IVSeqLen() %>%
  removeMultimappers.SRlen() %>%
  removeMultimappers.IVSeqLen() %>%
  filter(spliceType == "cis") %>%
  mutate(L = nchar(substrateSeq),
         N = nchar(pepSeq)) %>%
  removeMultimappers.Type %>%
  tidyr::separate_rows(digestTimes, intensities, sep=";") %>%
  rename(digestTime = digestTimes,
         intensity = intensities) %>%
  mutate(intensity = as.numeric(intensity),
         digestTime = as.numeric(digestTime))

pos = str_split_fixed(df_filter$positions,coll("_"),Inf)[,c(1:4)]
pos = apply(pos,2,as.numeric) %>% as.data.frame()
names(pos) = c("pos1","pos2","pos3","pos4")

MASTER = cbind(df_filter, pos) %>%
  na.omit() %>%
  mutate(contains_term = ifelse(pos4 == L, "C", NA),
         contains_term = ifelse(pos1 == 1, "N", contains_term)) %>%
  # filter(contains_term == "N" | contains_term == "C") %>%
  mutate(sr1len = ifelse(contains_term == "N", pos2-pos1+1, NA),
         sr2len = ifelse(contains_term == "C", pos4-pos3+1, NA)) %>%
  select(substrateID, pepSeq, contains_term,sr1len,sr2len) %>%
  tidyr::gather(sr,len,-contains_term, -substrateID, -pepSeq) %>%
  na.omit() %>%
  group_by(contains_term, sr, len) %>%
  mutate(n = n())

MASTER = left_join(MASTER, df_filter)

MASTERk = MASTER %>%
  filter(contains_term == "C" & len < 3)
plotKinetics(MASTERk, outfile = "results/termini/singleAASR2/kinetics_SR2_1-2aa.pdf", meanTech = T, earlyOnly = T, sortByInt = T)
plotKinetics(MASTER  %>% filter(contains_term == "N" & len < 3), outfile = "results/termini/singleAASR2/kinetics_SR1_1-2aa.pdf", meanTech = T, earlyOnly = T, sortByInt = T)


# plot summarised kinetics
MASTER  = MASTER %>%
  mutate(group = ifelse(contains_term == "C" & len < 3, "short_SR2", "all_other"),
         group = ifelse(contains_term == "N" & len < 3, "short_SR1", group))

# mean over bio rep
both = MASTER %>%
  group_by(substrateID,pepSeq,digestTime) %>%
  mutate(intensity = mean(intensity)) %>%
  filter(biological_replicate == 1)

# scale between 0 and 1
both = both %>%
  group_by(substrateID,pepSeq) %>%
  filter(group %in% c("short_SR1","short_SR2")) %>%
  mutate(int_norm = (intensity - min(intensity))/(max(intensity) - min(intensity)))

kinets = ggplot(both, aes(x = digestTime, y = int_norm, col = group)) +
  geom_smooth(method = "loess", se = T)+
  scale_color_manual(values = c("gray","lightblue")) +
  ggtitle("summarised kinetics") +
  ylab("scaled intensity") +
  xlab("time [hrs]")
kinets
ggsave(filename = "results/termini/singleAASR2/kinetics_summarised_shortSRs.png", plot = kinets,
       height = 4, width = 6, dpi = "retina")



## intensity distributions
suma = MASTER %>% filter(digestTime == 4 & intensity > 0) %>% group_by(group) %>% summarise(med = median(log10(intensity+1)))
intdist = ggplot(MASTER %>% filter(digestTime == 4 & intensity > 0), aes(x = log10(intensity+1), col = group)) +
  geom_density() +
  geom_vline(xintercept = suma$med, color =  c("olivedrab4","gray","lightblue"), lty = "dashed") +
  scale_color_manual(values = c("olivedrab4","gray","lightblue")) +
  xlab("log10(MS1 intensity + 1)")+
  ggtitle("intensity distribution of fwd. cis PSPs")

ggsave(filename = "results/termini/singleAASR2/intensity_distr.png", plot = intdist,
       height = 4, width = 6, dpi = "retina")


