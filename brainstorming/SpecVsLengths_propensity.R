### protein PCPS ###
# description:  specificity vs lengths - propensity scale approach
# input:        iaPSIre assigned and quantified peptides - amino acids and lengths extracted
# output:       specificity of the peptides vs SR/intervening sequence length
# author:       HPR

library(dplyr)
library(stringr)
library(seqinr)
library(protr)
library(Peptides)
library(transport)
library(ggplot2)

theme_set(theme_classic())

### INPUT ###
load("results/SpecVsLengths/DATA_proteins.RData")
load("results/SpecVsLengths/rndDATA_proteins.RData")
load("results/SpecVsLengths/DATA_polypeps.RData")
load("results/SpecVsLengths/rndDATA_polypeps.RData")


### MAIN PART ###
# ----- relevant information -----
AA = c("A","D","E","F","G","H","K","L","N","P","Q","R","S","T","V","W","Y", "M", "C", "X")
AAchar = c("P","G","A","V","L","M","F","Y","W","H","R","K","D","E","N","Q","S","T","C", "X")

SR1pos = c("P4"=-3, "P3"=-2, "P2"=-1, "P1"=0,
           "P-1"=1, "P-2"=2, "P-3"=3, "P-4"=4)

SR2pos = c("P-4_"=-4, "P-3_"=-3, "P-2_"=-2, "P-1_"=-1,
           "P1_"=0, "P2_"=1, "P3_"=2, "P4_"=3)

PCPpos = c("P4"=-3, "P3"=-2, "P2"=-1, "P1"=0,
           "P1_"=1, "P2_"=2, "P3_"=3, "P4_"=4)

SRnames = c(names(SR1pos), names(SR2pos))
types = c("cis", "revCis", "trans", "PCP")


# ----- encode splice sites: amino acid composition -----

getAAGroups = function(DBaa) {
  pos_pasted = do.call(paste, c(DBaa[,SRnames], sep=""))
  pos_pasted = str_replace_all(pos_pasted, coll(" "), "X")
  
  # get splice reactants
  pos = str_split_fixed(DBaa$positions,coll("_"),Inf)
  pos = apply(pos,2,as.numeric)
  DBaa$sr1 = str_sub(DBaa$substrateSeq,start = pos[,1], end = pos[,2])
  DBaa$sr2 = str_sub(DBaa$substrateSeq,start = pos[,3], end = pos[,4])
  
  getaapropensities = function(sequences) {
    
    aagroups = sapply(sequences, function(x){
      aaComp(x) %>% unlist() %>% t()
    }) %>% 
      t() %>%
      as.data.frame()
    
    aagroups = aagroups[,c(10:18)]
    names(aagroups) = rownames(aaComp("A")[[1]])
    
    return(aagroups)
  }
  
  aagroups_sites = getaapropensities(pos_pasted)
  aagroups_peps = getaapropensities(DBaa$pepSeq)
  aagroups_sr1 = getaapropensities(DBaa$sr1)
  aagroups_sr2 = getaapropensities(DBaa$sr2)
  
  return(list(sites = aagroups_sites,
              peps = aagroups_peps,
              sr1 = aagroups_sr1,
              sr2 = aagroups_sr2))
}



# ----- iterate lengths and calculate pairwise similarities -----
# matrix multiplication
matmult = function(v1 = "", v2 = ""){
  return(as.numeric(v1) %*% as.numeric(v2))
}
# dot product
dot_product = function(v1 = "", v2 = ""){
  p = matmult(v1, v2)/(sqrt(matmult(v1, v1)) * sqrt(matmult(v2, v2)))
  return(p)
}

getSimilarity = function(k, aagroups) {
  
  if (length(k) > 0) {
    
    Sims = matrix(NA,length(k),length(k))
    pb = txtProgressBar(min = 0, max = length(k), style = 3)
    for (i in 1:length(k)) {
      setTxtProgressBar(pb,i)
      for (j in i:length(k)) {
        
        # Sims[i,j] = dot_product(as.numeric(aagroups[k[i],]), as.numeric(aagroups[k[j],])) %>% as.numeric()
        Sims[i,j] = wasserstein1d(as.numeric(aagroups[k[i],]), as.numeric(aagroups[k[j],]))
      }
    }
    distr = Sims[upper.tri(Sims)] %>% as.numeric()
    
  } else {
    distr = NA
  }
  
  
  return(distr)
}


# ----- apply for SR1 length -----

iterateStuff = function(DBaa, type, aagroups, sample = F) {
  
  # --- SR1
  SR1lengths = DBaa$SR1Len %>% unique() %>% sort()
  SR1_sim = lapply(SR1lengths, function(x){
    print(x)
    k = which(DBaa$SR1Len == x & DBaa$spliceType == type)
    if (length(k) > 1) {
      sim = data.frame(type = type, len = x, gr = "SR1_length", n = length(k),
                       sim = getSimilarity(k, aagroups))
    } else {
      sim = data.frame(type = type, len = x, gr = "SR1_length", n = length(k),
                       sim = NA)
    }
    
    return(sim)
  })
  SR1_sim = plyr::ldply(SR1_sim)
  
  
  # --- SR2
  SR2lengths = DBaa$SR2Len %>% unique() %>% sort()
  SR2_sim = lapply(SR2lengths, function(x){
    print(x)
    k = which(DBaa$SR2Len == x & DBaa$spliceType == type)
    if (length(k) > 1) {
      sim = data.frame(type = type, len = x, gr = "SR2_length", n = length(k),
                       sim = getSimilarity(k, aagroups))
    } else {
      sim = data.frame(type = type, len = x, gr = "SR2_length", n = length(k),
                       sim = NA)
    }
    
    return(sim)
  })
  SR2_sim = plyr::ldply(SR2_sim)
  
  
  # --- shorter SR
  IVlengths = DBaa$IVlen %>% unique() %>% sort()
  IVlen_sim = lapply(IVlengths, function(x){
    print(x)
    k = which(DBaa$IVlen == x & DBaa$spliceType == type)
    if (length(k) > 1) {
      sim = data.frame(type = type, len = x, gr = "IV_length", n = length(k),
                       sim = getSimilarity(k, aagroups))
    } else {
      sim = data.frame(type = type, len = x, gr = "IV_length", n = length(k),
                       sim = NA)
    }
    
    return(sim)
  })
  IVlen_sim = plyr::ldply(IVlen_sim)
  
  out = list(SR1dep = SR1_sim,
             SR2dep = SR2_sim,
             IVdep = IVlen_sim)
  return(out)
}

# ----- break down complexity and plot -----
BreakDownAndPlot = function(DBaa, type, suffix="") {
  
  aagroups = getAAGroups(DBaa)
  
  SR1feat = iterateStuff(DBaa, type, aagroups = aagroups[["sr1"]])
  SR2feat = iterateStuff(DBaa, type, aagroups = aagroups[["sr2"]])
  Pepfeat = iterateStuff(DBaa, type, aagroups = aagroups[["peps"]])
  
  SR1feat_df = plyr::ldply(SR1feat) %>%
    mutate(features_of = "SR1feat")
  SR2feat_df = plyr::ldply(SR2feat) %>%
    mutate(features_of = "SR2feat")
  Pepfeat_df = plyr::ldply(Pepfeat) %>%
    mutate(features_of = "Pepfeat")
  
  ALL = rbind(SR1feat_df,SR2feat_df,Pepfeat_df)
  
  # ALL %>%
  #   group_by(type, gr, len, features_of) %>%
  #   summarise(mu = mean(sim)) %>% View()
  
  sr1lenplot = ALL %>%
    filter(gr == "SR1_length") %>%
    mutate(features_of = factor(features_of, levels = c("SR1feat","SR2feat","Pepfeat"))) %>%
    ggplot(aes(x=factor(len), y=sim, fill = features_of)) +
    geom_text(aes(label = factor(n)), y=-1, size = 2, stat = "count") +
    geom_boxplot() +
    scale_fill_manual(values = c("gray","lightblue","lightgreen")) +
    xlab("length") +
    ylab("divergence") +
    ggtitle("SR1 length") +
    facet_wrap(~features_of, scales = "free")
  
  sr2lenplot = ALL %>%
    filter(gr == "SR2_length") %>%
    mutate(features_of = factor(features_of, levels = c("SR1feat","SR2feat","Pepfeat"))) %>%
    ggplot(aes(x=factor(len), y=sim, fill = features_of)) +
    geom_text(aes(label = factor(n)), y=-1, size = 2, stat = "count") +
    geom_boxplot() +
    scale_fill_manual(values = c("gray","lightblue","lightgreen")) +
    xlab("length") +
    ylab("divergence") +
    ggtitle("SR2 length") +
    facet_wrap(~features_of, scales = "free")
  
  # intervening sequence --> summarise stronger
  ivlenplot = ALL %>%
    filter(gr == "IV_length") %>% 
    mutate(len_round = ifelse(len > 20, round(len, -1), len),
           features_of = factor(features_of, levels = c("SR1feat","SR2feat","Pepfeat"))) %>%
    filter(len_round <= 50) %>%
    ggplot(aes(x=factor(len_round), y=sim, fill = features_of)) +
    geom_text(aes(label = factor(n)), y=-1, size = 2, stat = "count") +
    geom_boxplot() +
    scale_fill_manual(values = c("gray","lightblue","lightgreen")) +
    xlab("length") +
    ylab("divergence") +
    ggtitle("intervening sequence length") +
    facet_wrap(~features_of, scales = "free")
  
  
  lenplot = gridExtra::grid.arrange(sr1lenplot, sr2lenplot, ivlenplot, nrow=3)
  ggsave(filename = paste0("results/SpecVsLengths/",type,"_propensity",suffix,".pdf"),
         plot = lenplot, height = 17, width = 19, dpi = "retina")
  
  
}

BreakDownAndPlot(prots_DBaa, type = "cis", suffix = "_proteins")
BreakDownAndPlot(prots_DBaa, type = "revCis", suffix = "_proteins")
BreakDownAndPlot(prots_DBaa, type = "trans", suffix = "_proteins")

BreakDownAndPlot(polyp_DBaa, type = "cis", suffix = "_polypeps")
BreakDownAndPlot(polyp_DBaa, type = "revCis", suffix = "_polypeps")
BreakDownAndPlot(polyp_DBaa, type = "trans", suffix = "_polypeps")


### OUTPUT ###


