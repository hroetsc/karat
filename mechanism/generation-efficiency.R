### karat projetc - PCPS mechanism ###
# description:  relation between PCP and PSP generation efficiency
# input:        quantitative datasets (polypeptides and proteins)
# output:       relation between PCP and PSP generation efficiency
# author:       HPR

library(dplyr)
library(stringr)
library(ggplot2)
library(grid)
library(gridExtra)
source("src/invitroSPI_utils.R")
source("src/plotting_utils.R")
source("../brainstorming/src/number-of-products.R")

theme_set(theme_classic())

### INPUT ###
load("data/invitroSPI.RData")
load("../../proteinsPCPS/new/data/aSPIre.RData")
proteins = Kinetics
load("data/aSPIre.RData")
polypeps = Kinetics

sample_list_polypeps1 = read.csv("invitroSPI+aSPIre/data/sample_list_aSPIre.csv", stringsAsFactors = F)
sample_list_polypeps2 = read.csv("aSPIre_Mut/data/sample_list_aSPIre.csv", stringsAsFactors = F)
sample_list_proteins = read.csv("../../aSPIre/data/sample_list.csv", stringsAsFactors = F)


### MAIN PART ###
# ----- preprocessing -----
# per raw file / time point
nm = intersect(names(sample_list_polypeps1), intersect(names(sample_list_polypeps2), names(sample_list_proteins)))
sample_list = rbind(sample_list_polypeps1[,nm],
                    sample_list_polypeps2[,nm],
                    sample_list_proteins[,nm]) %>%
  rename(scanTime = digestTime)

nm = intersect(names(proteins), names(polypeps))
quant = rbind(proteins[,nm] %>% mutate(seqtype = "proteins"),
              polypeps[,nm] %>% mutate(seqtype = "polypeptides")) %>%
  tidyr::separate_rows(assignedScans, sep = ";") %>%
  mutate(raw_file = paste0(str_extract_all(assignedScans, "^[:graph:]+(?=_[:digit:]+$)", simplify = T) %>% as.character(), ".raw")) %>%
  left_join(sample_list, by = c("substrateID","raw_file")) %>%
  rename(biological_replicate = biological_replicate.x,
         substrateSeq = substrateSeq.x) %>%
  select(-biological_replicate.y, -substrateSeq.y)
quant$substrateID = paste0(quant$substrateID,"_",quant$scanTime)

ProteasomeDB = ProteasomeDB %>% mutate(seqtype = "polypeptides")
ProteasomeDB$substrateID = paste0(ProteasomeDB$substrateID,"_",ProteasomeDB$digestTime)

nm = intersect(names(ProteasomeDB), names(quant))
# identifications
ALL = rbind(ProteasomeDB[,nm], quant[,nm]) %>%
  ILredundancy() %>%
  disentangleMultimappers.Type() %>%
  removeMultimappers.Type() %>%
  uniquePeptides()

# ----- get theoretical number of peptides/splice sites -----
pepLen_prots = seq(7,30)
pepLen_polypeps = seq(5,40)

getTheoreticalNumberofPepsAndSites = function(DB, pepLen) {
  
  THEORETICAL = DB %>%
    mutate(L = nchar(substrateSeq)) %>%
    group_by(substrateSeq, L, spliceType) %>%
    summarise() %>%
    mutate(N = paste(pepLen, collapse = ";")) %>%
    tidyr::separate_rows(N, sep = ";") %>%
    mutate(L = as.numeric(L),
           N = as.numeric(N))
  
  cis = THEORETICAL %>%
    filter(spliceType == "cis") %>%
    ungroup %>% rowwise() %>%
    mutate(numPeps = numCis(L,N),
           numSites = numCis_sites(L,N))
  
  revCis = THEORETICAL %>%
    filter(spliceType == "revCis") %>%
    ungroup %>% rowwise() %>%
    mutate(numPeps = numRevCis(L,N),
           numSites = numRevCis_sites(L,N))
  
  trans = THEORETICAL %>%
    filter(spliceType == "trans") %>%
    ungroup %>% rowwise() %>%
    mutate(numPeps = numTrans(L,N),
           numSites = numTrans_sites(L,N))
  
  pcp = THEORETICAL %>%
    filter(spliceType == "PCP") %>%
    ungroup %>% rowwise() %>%
    mutate(numPeps = numPCP(L,N),
           numSites = numCleavageSites(L,N))
  
  ALL = rbind(cis,revCis,trans,pcp) %>% as.data.frame()
  ALL$numPeps[ALL$numPeps < 0] = 0
  ALL$numSites[ALL$numSites < 0] = 0
  
  ALL = ALL %>%
    mutate(spliceType = ifelse(spliceType %in% c("cis","revCis"), "allcis", spliceType)) %>%
    group_by(substrateSeq, L, spliceType) %>%
    summarise(numPeps = sum(numPeps),
              numSites = sum(numSites)) %>%
    mutate(dataset = "theoretically")
  
  ALL2 = ALL %>% 
    filter(spliceType != "PCP") %>% 
    group_by(substrateSeq,L) %>% 
    summarise(numPeps = sum(numPeps), numSites = sum(numSites)) %>% 
    mutate(spliceType = "PSP", dataset = "theoretically") %>%
    select(substrateSeq,L,spliceType,numPeps,numSites,dataset)
  
  return(rbind(ALL,ALL2) %>% as.data.frame())
}

THERORETICAL_prot = getTheoreticalNumberofPepsAndSites(DB=ALL[ALL$seqtype == "proteins", ], pepLen = pepLen_prots)
THERORETICAL_polypeps = getTheoreticalNumberofPepsAndSites(DB=ALL[ALL$seqtype == "polypeptides", ], pepLen = pepLen_polypeps)


# ----- get actual number of peptides/splice sites ------

getActualNumberofPepsAndSites = function(DB) {
  
  DB$pasted = NA
  pcpidx = which(DB$spliceType == "PCP")
  pspidx = which(DB$spliceType != "PCP")
  
  pos = str_split_fixed(DB$positions, "_", Inf)
  DB$pasted[pspidx] = paste(pos[pspidx,2],pos[pspidx,3], sep = "_")
  DB$pasted[pcpidx] = pos[pcpidx,2]
  
  # number of pep
  ALL = DB %>%
    mutate(N = nchar(pepSeq),
           L = nchar(substrateSeq),
           spliceType = ifelse(spliceType %in% c("cis","revCis"), "allcis", spliceType)) %>%
    group_by(substrateID, substrateSeq, L, spliceType) %>%
    summarise(numPeps = n(),
              numSites = length(unique(pasted))) %>%
    mutate(dataset = "identified")
  
  ALL2 = ALL %>% 
    filter(spliceType != "PCP") %>% 
    group_by(substrateID, substrateSeq, L) %>% 
    summarise(numPeps = sum(numPeps), numSites = sum(numSites)) %>% 
    mutate(spliceType = "PSP", dataset = "identified") %>%
    select(substrateID,substrateSeq,L,spliceType,numPeps,numSites,dataset)
  
  return(rbind(ALL,ALL2) %>% as.data.frame())
}

OBSERVED_prot = getActualNumberofPepsAndSites(DB=ALL[ALL$seqtype == "proteins", ])
OBSERVED_polypeps = getActualNumberofPepsAndSites(DB=ALL[ALL$seqtype == "polypeptides", ])


# ----- theoretical vs. observed -----
# fraction of observed peptides
# fraction of observed splice sites

# proteins
PROT = OBSERVED_prot %>%
  rename(numPeps_observed = numPeps,
         numSites_observed = numSites) %>%
  select(-dataset) %>%
  left_join(THERORETICAL_prot) %>%
  mutate(pepFrac = numPeps_observed/numPeps,
         siteFrac = numSites_observed/numSites)

# polypeptides
POLYPEPS = OBSERVED_polypeps %>%
  rename(numPeps_observed = numPeps,
         numSites_observed = numSites) %>%
  select(-dataset) %>%
  left_join(THERORETICAL_polypeps) %>%
  mutate(pepFrac = numPeps_observed/numPeps,
         siteFrac = numSites_observed/numSites)

# merge
FRACOBSERVED = rbind(PROT %>% mutate(sequences = "proteins"),
                     POLYPEPS %>% mutate(sequences = "polypeptides"))


# correlate fractions of observed peptide products
pepsobserved = FRACOBSERVED %>%
  select(substrateID, spliceType, sequences, pepFrac) %>%
  tidyr::spread(spliceType, pepFrac, fill = 0)

peps = ggplot(pepsobserved, aes(x = PCP, y = PSP)) +
  geom_point() +
  geom_smooth(method = "loess") +
  xlab("fraction of PCPs observed") +
  ylab("fraction of PSPs observed") +
  ggtitle("fraction of observed peptide products", subtitle = "per substrate and time point") +
  facet_wrap(~sequences, scales = "free")
  

# correlate fractions of observed peptide products
sitesobserved = FRACOBSERVED %>%
  select(substrateID, spliceType, sequences, siteFrac) %>%
  tidyr::spread(spliceType, siteFrac, fill = 0)

sites = ggplot(sitesobserved, aes(x = PCP, y = PSP)) +
  geom_point() +
  geom_smooth(method = "loess") +
  xlab("fraction of cleavage sites observed") +
  ylab("fraction of splice sites observed") +
  ggtitle("fraction of observed cleavage/splice sites", subtitle = "per substrate and time point") +
  facet_wrap(~sequences, scales = "free")

ggsave("results/SR-as-PCP/generation_efficiency.png", plot = grid.arrange(peps, sites, nrow = 2), height = 10, width = 10, dpi = "retina")


# ----- plot mean intensities -----
Q = quant %>%
  select(substrateSeq, pepSeq, biological_replicate, digestTimes, intensities, productType, spliceType, positions, seqtype) %>%
  unique() %>%
  tidyr::separate_rows(digestTimes, intensities, sep = ";") %>%
  mutate(digestTimes = as.numeric(digestTimes),
         intensities = as.numeric(intensities))

Qsum = Q %>%
  mutate(intensity = log10(intensities),
         intensity = ifelse(!is.finite(intensity), NA, intensity)) %>%
  group_by(seqtype, substrateSeq, digestTimes, productType, biological_replicate) %>%
  summarise(mean_amount = mean(intensity, na.rm = T)) %>%
  na.omit() %>%
  tidyr::spread(productType, mean_amount)

pcc_poly = cor(Qsum$PCP[Qsum$seqtype == "polypeptides"], Qsum$PSP[Qsum$seqtype == "polypeptides"]) %>% round(4)
pcc_prots = cor(Qsum$PCP[Qsum$seqtype == "proteins"], Qsum$PSP[Qsum$seqtype == "proteins"]) %>% round(4)

ggplot(Qsum, aes(x = PCP, y = PSP)) +
  geom_point() +
  geom_smooth(method = "glm") +
  xlab("mean log10 MS1 intensity PCPs") +
  ylab("mean log10 MS1 intensity PSPs") +
  ggtitle("correspondence between MS1 intensities", subtitle = paste0("per substrate and time point and bio rep \nPCCs: ", pcc_poly, ", ", pcc_prots)) +
  facet_wrap(~seqtype, scales = "free")

ggsave("results/SR-as-PCP/MS1intensity.png", plot = last_plot(), height = 5, width = 10, dpi = "retina")

# --- check if correlation is stronger for those PSPs that are detected as PCPs
load("results/SR-as-PCP/MASTERinfo.RDAta")

Qfilter = quant %>%
  filter(pepSeq %in% unique(MASTERinfo$pepSeq[MASTERinfo$SR1used > 0 | MASTERinfo$SR2used > 0]) | productType == "PCP") %>%
  select(substrateSeq, pepSeq, biological_replicate, digestTimes, intensities, productType, spliceType, positions, seqtype) %>%
  unique() %>%
  tidyr::separate_rows(digestTimes, intensities, sep = ";") %>%
  mutate(digestTimes = as.numeric(digestTimes),
         intensities = as.numeric(intensities))

Qsumfilter = Qfilter %>%
  mutate(intensity = log10(intensities),
         intensity = ifelse(!is.finite(intensity), NA, intensity)) %>%
  group_by(seqtype, substrateSeq, digestTimes, productType, biological_replicate) %>%
  summarise(mean_amount = mean(intensity, na.rm = T)) %>%
  na.omit() %>%
  tidyr::spread(productType, mean_amount)

pcc_poly = cor(Qsumfilter$PCP[Qsumfilter$seqtype == "polypeptides"], Qsumfilter$PSP[Qsumfilter$seqtype == "polypeptides"]) %>% round(4)
pcc_prots = cor(Qsumfilter$PCP[Qsumfilter$seqtype == "proteins"], Qsumfilter$PSP[Qsumfilter$seqtype == "proteins"]) %>% round(4)

ggplot(Qsumfilter, aes(x = PCP, y = PSP)) +
  geom_point() +
  geom_smooth(method = "glm") +
  xlab("mean log10 MS1 intensity PCPs") +
  ylab("mean log10 MS1 intensity PSPs") +
  ggtitle("correspondence between MS1 intensities", subtitle = paste0("per substrate and time point and bio rep \nPCCs: ", pcc_poly, ", ", pcc_prots)) +
  facet_wrap(~seqtype, scales = "free")


### OUTPUT ###


