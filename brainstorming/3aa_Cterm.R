### PCPS brainstorming ###
# description:  peptides that carry the substrate's C-terminal amino acids
# input:        Roetschke et al. SciData Proteasome DB, IDP protein data set
#               processed as in NandCterm.R
#               + all theoretically possible peptides
# output:       overview of peptide products the substrate's 1-2 C-terminal amino acids
# author:       HPR

library(dplyr)
library(stringr)
library(ggplot2)
source("../../proteinsPCPS/new/src/invitroSPI_utils.R")
source("src/number-of-products.R")

theme_set(theme_classic())

### INPUT ###
load("../../proteinsPCPS/new/data/aSPIre.RData")
polypeps = read.csv("../../invitroSPI/revision/results/ProteasomeDB/ProteasomeDB.csv", stringsAsFactors = F)

# think before running!
load("data/no-theor_polypeps.RData")
# load("data/no-theor_prots.RData")


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

# ----- get peptides that have the substrate's C term as SR -----

getCtermSR = function(df, obs=T) {
  
  if (obs) {
    df_filter = df %>%
      filter(spliceType %in% c("cis","revCis","trans")) %>%
      mutate(L = nchar(substrateSeq),
             N = nchar(pepSeq)) %>%
      removeMultimappers.Type() %>%
      uniquePeptides()
    
  } else {
    
    df_filter = plyr::ldply(df)
    
  }
  
  pos = str_split_fixed(df_filter$positions,coll("_"),Inf)[,c(1:4)]
  pos = apply(pos,2,as.numeric) %>% as.data.frame()
  names(pos) = c("pos1","pos2","pos3","pos4")
  
  MASTER = cbind(df_filter, pos) %>%
    na.omit() %>%
    mutate(psbl = ifelse(spliceType %in% c("cis","trans") & pos1 >= 1 & pos4 <= L, T, F),
           psbl = ifelse(spliceType == "revCis" & pos3 >= 1 & pos2 <= L, T, psbl)) %>%  # check the all possible function!
    filter(psbl) %>%
    mutate(contains_Cterm = ifelse((spliceType %in% c("cis","trans") & pos4 == L) | (spliceType == "revCis" & pos2 == L), T, F),
           srlen = ifelse(spliceType %in% c("cis","trans"), pos4-pos3+1, pos2-pos1+1),
           contraMichaux = ifelse(contains_Cterm & srlen <= 2, T,F))
  
  # MASTER[which(MASTER$contraMichaux), ] %>% View()
  
  abund = MASTER %>%
    group_by(spliceType, L, contraMichaux) %>%
    summarise(n = n()) %>%
    tidyr::spread(contraMichaux,n,fill = 0) %>%
    tidyr::gather(contraMichaux,n, -L, -spliceType) %>%
    ungroup() %>% group_by(spliceType, L) %>%
    mutate(freq = n/sum(n))
  
  if (obs) {
    abund$dataset = "identified"
  } else {
    abund$dataset = "theoretically"
  }
  
  return(abund)
}

obs_polypep_frac = getCtermSR(polypeps)
theor_polypep_frac = getCtermSR(theor_poly, obs = F)

obs_protein_frac = getCtermSR(proteins)

# ----- plotting -----
ALL_polypep = rbind(obs_polypep_frac, theor_polypep_frac) %>%
  filter(contraMichaux == T)

plot_polypep = ggplot(ALL_polypep, aes(x = spliceType, y = freq, fill = dataset)) +
  geom_boxplot() +
  scale_fill_manual("data set", values = c("palegreen","steelblue")) +
  ylab("frequency") +
  ggtitle("substrate's last 2 aa used as SR", subtitle = "polypeptides") +
  theme(axis.text.x = element_text(angle = 90))

ggsave(filename = "results/3aaCterminus_polypeptides.png", plot = plot_polypep,
       height = 6, width = 6, dpi = "retina")


