### PCPS brainstorming ###
# description:  deduce the minimal splice-reactant length required for the proteasome
# input:        Roetschke et al. SciData Proteasome DB, IDP protein data set
#               processed as in NandCterm.R
#               + all theoretically possible peptides
# output:       length distributions of cis-peptides carrying the substrate's termini
# author:       HPR

library(dplyr)
library(stringr)
library(ggplot2)
source("../../proteinsPCPS/new/src/plotting_utils.R")
source("../../proteinsPCPS/new/src/invitroSPI_utils.R")
source("src/number-of-products.R")

theme_set(theme_classic())

### INPUT ###
load("../../proteinsPCPS/new/data/aSPIre.RData")
polypeps = read.csv("../../invitroSPI/revision/results/ProteasomeDB/ProteasomeDB.csv", stringsAsFactors = F)


### MAIN PART ###
# ----- data preprocessing -----
proteins = Kinetics

proteins = proteins %>%
  distinct(substrateID, pepSeq, .keep_all = T) %>%
  tidyr::separate_rows(spectralAngles, assignedScans, sep=";") %>%
  ILredundancy() %>%
  disentangleMultimappers.Type() %>%
  disentangleMultimappers.SRlen() %>%
  disentangleMultimappers.IVSeqLen()

polypeps = polypeps %>%
  ILredundancy() %>%
  filterPepLength() %>%
  disentangleMultimappers.Type() %>%
  disentangleMultimappers.SRlen() %>%
  disentangleMultimappers.IVSeqLen() %>%
  remSynthErrors() %>%
  filter20Sstandard() %>%
  DBcosmetics()


# ----- fetch cis-spliced carrying the substrate's termini -----

fetchSRs = function(df) {
  
  df_filter = df %>%
    filter(spliceType == "cis") %>%
    mutate(L = nchar(substrateSeq),
           N = nchar(pepSeq)) %>%
    removeMultimappers.Type() %>%
    uniquePeptides()
  
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

prots_SR = fetchSRs(proteins)
polypeps_SR = fetchSRs(polypeps)

# ----- get length distribution of splice-reactants at termini -----


pr = ggplot(prots_SR, aes(x=len, y = n, fill = sr)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_continuous(breaks = seq(1,30,2)) +
  scale_fill_manual("", values = c("gray","lightblue"), labels = c("SR1 (N-terminus)","SR2 (C-terminus)")) +
  xlab("splice-reactant length (aa residues)") +
  ylab("counts") +
  ggtitle("proteins")

po = ggplot(polypeps_SR, aes(x=len, y = n, fill = sr)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_continuous(breaks = seq(1,30,2)) +
  scale_fill_manual("", values = c("gray","lightblue"), labels = c("SR1 (N-terminus)","SR2 (C-terminus)")) +
  xlab("splice-reactant length (aa residues)") +
  ylab("counts") +
  ggtitle("polypeptides")


srs = gridExtra::grid.arrange(po,pr,ncol=2)
ggsave(filename = "results/SRlenAtTerm.png", plot = srs, height = 5, width = 10, dpi = "retina")
# "length distribution of cis splice-reactants", subtitle = "carrying the substrate's N-terminus (SR1) or C-terminus (SR2)"
