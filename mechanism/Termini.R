### karat projetc - PCPS mechanism ###
# description:  peptides carrying substrate's N/C term
# input:        qualitative data set: Roetschke et al. SciData, EGFR, WT sequences of WT/Mut
# output:       overview of peptide products carrying N/C term
# author:       HPR


library(dplyr)
library(stringr)
library(ggplot2)
source("src/invitroSPI_utils.R")
source("../brainstorming/src/number-of-products.R")

theme_set(theme_classic())

### INPUT ###
load("data/invitroSPI.RData")

### MAIN PART ###
suppressWarnings(dir.create("results/termini/"))

# ----- data preprocessing -----

polypeps = ProteasomeDB %>%
  ILredundancy() %>%
  filterPepLength() %>%
  disentangleMultimappers.Type() %>%
  disentangleMultimappers.SRlen() %>%
  disentangleMultimappers.IVSeqLen() %>%
  remSynthErrors() %>%
  filter20Sstandard() %>%
  DBcosmetics()


# ----- get number of peptides for a given length -----
subLen_poly = polypeps$substrateSeq %>% nchar() %>% unique() %>% sort()

# polypeptides: 5 to 40 aa (MS limitations)
prodLen_poly = seq(5,40)

getAllPos = function(subLen, prodLen, shortV=F) {
  
  print("forward cis")
  cis = lapply(prodLen, function(N){
    x = lapply(subLen, function(L){
      data.frame(L=L, N=N, number_all=numCis(L,N), number_1term = numCis_1term(L,N), number_bothterm = numCis_bothterm(L,N))
    })
    x = plyr::ldply(x)
    return(x)
  })
  cis_df = plyr::ldply(cis)
  cis_df$spliceType = "cis"
  
  print("reverse cis")
  revcis = lapply(prodLen, function(N){
    x = lapply(subLen, function(L){
      data.frame(L=L, N=N, number_all=numRevCis(L,N), number_1term = numRevCis_1term(L,N), number_bothterm = numRevCis_bothterm(L,N))
    })
    x = plyr::ldply(x)
    return(x)
  })
  revcis_df = plyr::ldply(revcis)
  revcis_df$spliceType = "revCis"
  
  df = list(cis = cis_df,
            revCis = revcis_df)
  return(df)
}

theor_poly = getAllPos(subLen_poly, prodLen_poly)

# ----- extract peptides from actual data sets -----
# categories

extractPeps = function(df, obs=T) {
  
  if (obs) {
    df_filter = df %>%
      mutate(L = nchar(substrateSeq),
             N = nchar(pepSeq)) %>%
      removeMultimappers.Type() %>%
      filter(spliceType %in% c("cis","revCis")) %>%
      uniquePeptides()
    
    pos = str_split_fixed(df_filter$positions,coll("_"),Inf)[,c(1:4)]
    pos = apply(pos,2,as.numeric) %>% as.data.frame()
    names(pos) = c("pos1","pos2","pos3","pos4")
    
    CIS = cbind(df_filter, pos) %>%
      filter(spliceType == "cis") %>%
      na.omit() %>%
      mutate(category = ifelse(pos1 == 1 & pos4 < L, "cat1", NA),
             category = ifelse(pos1 > 1 & pos4 == L, "cat2", category),
             category = ifelse(pos1 == 1 & pos4 == L, "cat3", category),
             category = ifelse(pos1 > 1 & pos4 < L, "cat4", category))
    
    REVCIS = cbind(df_filter, pos) %>%
      filter(spliceType == "revCis") %>%
      na.omit() %>%
      mutate(category = ifelse(pos3 == 1 & pos2 < L, "cat1", NA),
             category = ifelse(pos3 > 1 & pos2 == L, "cat2", category),
             category = ifelse(pos3 == 1 & pos2 == L, "cat3", category),
             category = ifelse(pos3 > 1 & pos2 < L, "cat4", category))
    
    abund = rbind(CIS,REVCIS) %>%
      group_by(L,spliceType,category) %>%
      summarise(n = n()) %>%
      tidyr::spread(category,n,fill = 0) %>%
      tidyr::gather(category,n, -L, -spliceType) %>%
      ungroup() %>% group_by(L,spliceType) %>%
      mutate(freq = n/sum(n))
    
  } else {
    
    ALL = plyr::ldply(df)
    
    abund = ALL %>%
      mutate(number_noterm = number_all-number_1term-number_bothterm) %>%
      rename(cat1 = number_1term,
             cat3 = number_bothterm,
             cat4 = number_noterm) %>%
      mutate(cat2 = cat1) %>%
      select(-number_all) %>%
      tidyr::gather("category","n",cat1,cat2,cat3,cat4, -L, -N, -spliceType) %>%
      select(-.id) %>%
      mutate(category = as.character(category)) %>%
      group_by(L,spliceType,category) %>%
      summarise(n = sum(n)) %>%
      tidyr::spread(category,n, fill = 0) %>%
      tidyr::gather(category,n, -L, -spliceType) %>%
      ungroup() %>% group_by(L,spliceType) %>%
      mutate(freq = n/sum(n))
    
    
  }
  
  return(abund)
}

obs_poly_freq = extractPeps(polypeps)
theor_poly_freq = extractPeps(theor_poly, obs = F)


# ----- plotting -----
# polypeptides
obs_poly_freq$dataset = "identified"
theor_poly_freq$dataset = "theoretically"
ALL_poly = rbind(obs_poly_freq,theor_poly_freq)

plot_poly = ALL_poly %>%
  ggplot(aes(x=category,y=freq,fill=dataset)) +
  geom_boxplot() +
  ylim(0,1) +
  # geom_text(aes(label = L, col=dataset), position=position_jitter(width = .3), size = 2) +
  scale_color_manual("", values = c("black","black")) +
  scale_fill_manual("data set", values = c(plottingCols[["ProteasomeDB"]], plottingCols[["randomDB"]])) +
  scale_x_discrete(labels = c("N term", "C term", "N+C term", "no termini")) +
  ylab("frequency") +
  theme(axis.text.x = element_text(angle = 90)) +
  facet_wrap(~spliceType, scales = "free")

ggsave(filename = "results/termini/NandCterm_frequency.png", plot = plot_poly,
       height = 10, width = 16, units = "cm", dpi = "retina")


ALL_poly %>%
  group_by(spliceType, category, dataset) %>%
  summarise(median = median(freq),
            mean = mean(freq),
            std = sd(freq)) %>%
  print.data.frame()

