### PCPS brainstorming ###
# description:  peptides carrying substrate's N/C term in polypeptide and protein data set
# input:        Roetschke et al. SciData Proteasome DB, IDP protein data set
# output:       overview of peptide products carrying N/C term
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


# ----- get number of peptides for a given length -----
subLen_poly = polypeps$substrateSeq %>% nchar() %>% unique() %>% sort()
subLen_prots = proteins$substrateSeq %>% nchar() %>% unique() %>% sort()

# polypeptides: 5 to 40 aa (MS limitations)
prodLen_poly = seq(5,40)
prodLen_prots = seq(7,30)

getAllPos = function(subLen, prodLen, shortV=F) {
  
  print("forward cis")
  cis = lapply(prodLen, function(N){
    x = lapply(subLen, function(L){
      data.frame(L=L, N=N, positions=getCis(L,N))
    })
    x = plyr::ldply(x)
    return(x)
  })
  cis_df = plyr::ldply(cis)
  cis_df$spliceType = "cis"
  
  print("reverse cis")
  revcis = lapply(prodLen, function(N){
    x = lapply(subLen, function(L){
      data.frame(L=L, N=N, positions=getRevCis(L,N))
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
theor_prot = getAllPos(subLen_prots, prodLen_prots)

save(theor_poly, file = "data/no-theor_polypeps.RData")
save(theor_prot, file = "data/no-theor_prots.RData")


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
      mutate(category = if_else(pos1 == 1 & pos4 < L, "cat1", NA),
             category = if_else(pos1 > 1 & pos4 == L, "cat2", category),
             category = if_else(pos1 == 1 & pos4 == L, "cat3", category),
             category = if_else(pos1 > 1 & pos4 < L, "cat4", category))
    
    REVCIS = cbind(df_filter, pos) %>%
      filter(spliceType == "revCis") %>%
      na.omit() %>%
      mutate(category = if_else(pos3 == 1 & pos2 < L, "cat1", NA),
             category = if_else(pos3 > 1 & pos2 == L, "cat2", category),
             category = if_else(pos3 == 1 & pos2 == L, "cat3", category),
             category = if_else(pos3 > 1 & pos2 < L, "cat4", category))
    
  } else {
    
    cis_pos = str_split_fixed(df$cis$positions,coll("_"),Inf)[,c(1:4)]
    cis_pos = apply(cis_pos,2,as.numeric) %>% as.data.frame()
    names(cis_pos) = c("pos1","pos2","pos3","pos4")
    
    print("managed cis positions")
    
    revcis_pos = str_split_fixed(df$revCis$positions,coll("_"),Inf)[,c(1:4)]
    revcis_pos = apply(revcis_pos,2,as.numeric) %>% as.data.frame()
    names(revcis_pos) = c("pos1","pos2","pos3","pos4")
    
    print("managed rev cis positions")
    
    CIS = cbind(df$cis, cis_pos) %>%
      na.omit() %>%
      filter(pos1 >= 1 & pos4 <= L) %>%  # check the all possible function!
      mutate(category = if_else(pos1 == 1 & pos4 < L, "cat1", NA),
             category = if_else(pos1 > 1 & pos4 == L, "cat2", category),
             category = if_else(pos1 == 1 & pos4 == L, "cat3", category),
             category = if_else(pos1 > 1 & pos4 < L, "cat4", category))
    
    print("managed to get cis df")
    
    REVCIS = cbind(df$revCis, revcis_pos) %>%
      na.omit() %>%
      filter(pos3 >= 1 & pos2 <= L) %>%  # check the all possible function!
      mutate(category = if_else(pos3 == 1 & pos2 < L, "cat1", NA),
             category = if_else(pos3 > 1 & pos2 == L, "cat2", category),
             category = if_else(pos3 == 1 & pos2 == L, "cat3", category),
             category = if_else(pos3 > 1 & pos2 < L, "cat4", category))
    
    print("managed to get rev cis df")
  }
  
  
  abund = rbind(CIS,REVCIS) %>%
    group_by(L,spliceType,category) %>%
    summarise(n = n()) %>%
    tidyr::spread(category,n,fill = 0) %>%
    tidyr::gather(category,n, -L, -spliceType) %>%
    ungroup() %>% group_by(L,spliceType) %>%
    mutate(freq = n/sum(n))
  
  print("managed to summarise")
  
  return(abund)
}

obs_poly_freq = extractPeps(polypeps)
obs_prot_freq = extractPeps(proteins)

theor_poly_freq = extractPeps(theor_poly, obs = F)
theor_prot_freq = extractPeps(theor_prot, obs = F)


# ----- plotting -----
# polypeptides
obs_poly_freq$dataset = "identified"
theor_poly_freq$dataset = "theoretically"
ALL_poly = rbind(obs_poly_freq,theor_poly_freq)

plot_poly = ALL_poly %>%
  ggplot(aes(x=category,y=freq,fill=dataset)) +
  geom_boxplot() +
  # geom_text(aes(label = L, col=dataset), position=position_jitter(width = .3), size = 2) +
  scale_color_manual("", values = c("black","black")) +
  scale_fill_manual("data set", values = c("palegreen","steelblue")) +
  scale_x_discrete(labels = c("N term", "C term", "N+C term", "no termini")) +
  ylab("frequency") +
  theme(axis.text.x = element_text(angle = 90)) +
  facet_wrap(~spliceType)

ggsave(filename = "results/NandCterm_polypeps.png", plot = plot_poly,
       height = 6, width = 12, dpi = "retina")


# proteins
obs_prot_freq$dataset = "identified"
theor_prot_freq$dataset = "theoretically"
ALL_prot = rbind(obs_prot_freq,theor_prot_freq)

plot_prot = ALL_prot %>%
  ggplot(aes(x=category,y=freq,fill=dataset)) +
  geom_boxplot() +
  # geom_text(aes(label = L, col=dataset), position=position_jitter(width = .3), size = 2) +
  scale_color_manual("", values = c("black","black")) +
  scale_fill_manual("data set", values = c("palegreen","steelblue")) +
  scale_x_discrete(labels = c("N term", "C term", "N+C term", "no termini")) +
  ylab("frequency") +
  theme(axis.text.x = element_text(angle = 90)) +
  facet_wrap(~spliceType)

ggsave(filename = "results/NandCterm_proteins_attempt2.png", plot = plot_prot,
       height = 6, width = 12, dpi = "retina")

