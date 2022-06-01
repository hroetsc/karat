### karat projetc - PCPS mechanism ###
# description:  kinetics of reverse cis peptides carrying substrate's termini
# input:        quantitative data set: EGFR, WT sequences of WT/Mut
#               qualitative data set: Roetschke et al. SciData, EGFR, WT sequences of WT/Mut
# output:       kinetics of revCis: transpeptidation vs. condensation
# author:       HPR

library(dplyr)
library(stringr)
library(ggplot2)
library(uwot)
library(twosamples)
library(dgof)
source("src/invitroSPI_utils.R")
source("src/aSPIre_plotting.R")
source("src/plotting_utils.R")

theme_set(theme_classic())

### INPUT ###
load("data/aSPIre.RData")
load("data/invitroSPI.RData")

### MAIN PART ###
suppressWarnings(dir.create("results/termini/kinetics/"))

# ----- preprocessing -----

polypeps = Kinetics %>%
  ILredundancy() %>%
  disentangleMultimappers.Type() %>%
  disentangleMultimappers.SRlen() %>%
  disentangleMultimappers.IVSeqLen()

polypepsQual = ProteasomeDB %>%
  ILredundancy() %>%
  filterPepLength() %>%
  disentangleMultimappers.Type() %>%
  disentangleMultimappers.SRlen() %>%
  disentangleMultimappers.IVSeqLen() %>%
  remSynthErrors() %>%
  filter20Sstandard() %>%
  DBcosmetics()

# ----- extract termini information -----

getTermini = function(polypeps) {
  
  df_filter = polypeps %>%
    mutate(L = nchar(substrateSeq),
           N = nchar(pepSeq)) %>%
    filter(!grepl(";", positions)) %>%
    filter(spliceType %in% c("cis","revCis"))
  
  pos = str_split_fixed(df_filter$positions,coll("_"),Inf)[,c(1:4)]
  pos = apply(pos,2,as.numeric) %>% as.data.frame()
  names(pos) = c("pos1","pos2","pos3","pos4")
  
  REVCIS = cbind(df_filter, pos) %>%
    filter(spliceType == "revCis") %>%
    na.omit() %>%
    mutate(category = ifelse(pos3 == 1 & pos2 < L, "cat1", NA),
           category = ifelse(pos3 > 1 & pos2 == L, "cat2", category),
           category = ifelse(pos3 == 1 & pos2 == L, "cat3", category),
           category = ifelse(pos3 > 1 & pos2 < L, "cat4", category))
  
  transpeptidation = REVCIS[which(REVCIS$category == "cat1"), ]
  condensation = REVCIS[which(REVCIS$category %in% c("cat2","cat3")), ]
  
  if (all(c("digestTimes","intensities") %in% names(REVCIS))) {
    
    transpeptidation = transpeptidation %>%
      tidyr::separate_rows(digestTimes, intensities, sep=";") %>%
      rename(digestTime = digestTimes,
             intensity = intensities) %>%
      mutate(intensity = as.numeric(intensity),
             digestTime = as.numeric(digestTime))
    
    condensation = condensation %>%
      tidyr::separate_rows(digestTimes, intensities, sep=";") %>%
      rename(digestTime = digestTimes,
             intensity = intensities) %>%
      mutate(intensity = as.numeric(intensity),
             digestTime = as.numeric(digestTime))
  }
  
  print("transpeptidation")
  transpeptidation$pepSeq %>% unique() %>% length() %>% print()
  print("condensation")
  condensation$pepSeq %>% unique() %>% length() %>% print()
  
  return(list(transpeptidation = transpeptidation,
              condensation = condensation))
}

polypepsTerm = getTermini(polypeps)
polypepsQualTerm = getTermini(polypepsQual)

transpeptidation = polypepsTerm$transpeptidation
condensation = polypepsTerm$condensation

# ----- plot kinetics -----

plotKinetics(transpeptidation, outfile = "results/termini/kinetics/transpeptidation.pdf",
             meanTech = T, earlyOnly = F, sortByInt = T)
plotKinetics(condensation, outfile = "results/termini/kinetics/condensation.pdf",
             meanTech = T, earlyOnly = F, sortByInt = T)


# ----- compare overall characteristics -----

bothMech = rbind(transpeptidation %>% mutate(mechanism = "transpeptidation"),
                 condensation %>% mutate(mechanism = "condensation")) %>%
  filter(digestTime == 4)

# ----- intensity distributions=
intDistr = ggplot(bothMech, aes(x = mechanism, y = log10(intensity+1), fill = mechanism)) +
  geom_violin(draw_quantiles = c(0.5)) + 
  scale_fill_manual(values = c("blue","red")) +
  ggtitle("intensity distributions after 4 hours")
ggsave(filename = "results/termini/kinetics/Transpep-vs-Condens_intensity.png", plot = intDistr,
       height = 4, width = 4, dpi = "retina")

bothMech %>%
  group_by(mechanism) %>%
  summarise(n = n(),
            median = median(log10(intensity+1)),
            mean = mean(log10(intensity+1)),
            std = sd(log10(intensity+1))) %>%
  print.data.frame()

ad_test(log10(bothMech$intensity[bothMech$mechanism == "transpeptidation"]+1),
        log10(bothMech$intensity[bothMech$mechanism == "condensation"]+1))

# ----- peptide lengths
bothMechQual = rbind(polypepsQualTerm$transpeptidation %>% mutate(mechanism = "transpeptidation"),
                     polypepsQualTerm$condensation %>% mutate(mechanism = "condensation")) %>%
  mutate(N = nchar(pepSeq))

LenDistr = ggplot(bothMechQual, aes(x = mechanism, y = N, fill = mechanism)) +
  geom_violin(draw_quantiles = c(0.5)) + 
  scale_fill_manual(values = c("blue","red")) +
  ggtitle("peptide length") +
  ylab("peptide length (aa residues)")
ggsave(filename = "results/termini/kinetics/Transpep-vs-Condens_length.png", plot = LenDistr,
       height = 4, width = 4, dpi = "retina")

bothMechQual %>%
  group_by(mechanism) %>%
  summarise(n = n(),
            median = median(N),
            mean = mean(N),
            std = sd(N)) %>%
  print.data.frame()

ad_test(bothMech$N[bothMech$mechanism == "transpeptidation"],
        bothMech$N[bothMech$mechanism == "condensation"])

# ----- SR lengths
bothMech_SR = bothMechQual %>%
  disentangleMultimappers.SRlen() %>%
  removeMultimappers.SRlen() %>%
  mutate(sr1len = pos2-pos1+1,
         sr2len = pos4-pos3+1)

sr1 = ggplot(bothMech_SR, aes(x = mechanism, y = sr1len, fill = mechanism)) +
  geom_violin(draw_quantiles = c(0.5)) + 
  scale_fill_manual(values = c("blue","red")) +
  ggtitle("SR1 length") +
  ylab("SR1 length (aa residues)")

sr2 = ggplot(bothMech_SR, aes(x = mechanism, y = sr2len, fill = mechanism)) +
  geom_violin(draw_quantiles = c(0.5)) + 
  scale_fill_manual(values = c("blue","red")) +
  ggtitle("SR2 length") +
  ylab("SR2 length (aa residues)")

srlen = gridExtra::grid.arrange(sr1,sr2,ncol = 2)
ggsave(filename = "results/termini/kinetics/Transpep-vs-Condens_SRlength.png", plot = srlen,
       height = 4, width = 8, dpi = "retina")

# ----- cluster kinetics -----
both = rbind(transpeptidation %>% mutate(mechanism = "transpeptidation"),
             condensation %>% mutate(mechanism = "condensation"))

# mean over bio rep
both = both %>%
  group_by(substrateID,pepSeq,digestTime) %>%
  mutate(intensity = mean(intensity)) %>%
  filter(biological_replicate == 1)

# scale between 0 and 1
both = both %>%
  group_by(substrateID,pepSeq) %>%
  mutate(int_norm = (intensity - min(intensity))/(max(intensity) - min(intensity)))

# get differences between time points
bothChar = both %>%
  group_by(mechanism, substrateID, pepSeq) %>%
  summarise(diff0_4 = int_norm[digestTime == 4] - int_norm[digestTime == 0],
            diff0_1 = int_norm[digestTime == 1] - int_norm[digestTime == 0],
            diff1_2 = int_norm[digestTime == 2] - int_norm[digestTime == 1],
            diff2_3 = int_norm[digestTime == 3] - int_norm[digestTime == 2],
            diff3_4 = int_norm[digestTime == 4] - int_norm[digestTime == 3],
            # int_1 = intensity[digestTime == 1],
            #int_2 = int_norm[digestTime == 2],
            #int_3 = int_norm[digestTime == 3],
            int_4 = intensity[digestTime == 4]) %>%
  as.data.frame()

co = rep("red",nrow(bothChar))
co[bothChar$mechanism == "condensation"] = "blue"

y = rep(0, nrow(bothChar))
y[bothChar$mechanism == "condensation"] = 1

# UMAP embeddings
set.seed(41)
embedding = uwot::umap(bothChar %>% select(-mechanism, -substrateID, -pepSeq),
                       y = y,
                       target_weight = .7,
                       n_neighbors = 5,
                       metric = "cosine",
                       init = "agspectral",
                       ret_extra = T,
                       verbose = T)

png("results/termini/kinetics/Transpep-vs-Condens_embedding.png", height = 7, width = 7, units = "in", res = 300)
plot(x = embedding[,1], y = embedding[,2],
     pch = 16, col = add.alpha(co, .8),
     xlab = "UMAP 1", ylab = "UMAP 2",
     main = "embedding of generation kinetics",
     sub = "transpeptidation (red) + condensation (blue)")
legend("bottomleft",
       legend = c("transpeptidation", "condensation"),
       col = c("red", "blue"),
       pch = rep(16,2), bty = "n", cex = .8)
dev.off()


kinets = ggplot(both, aes(x = digestTime, y = int_norm, col = mechanism)) +
  geom_smooth(method = "loess")+
  scale_color_manual(values = c("blue","red")) +
  ggtitle("summarised kinetics") +
  ylab("scaled intensity") +
  xlab("time [hrs]")
ggsave(filename = "results/termini/kinetics/Transpep-vs-Condens_kinetics.png", plot = kinets,
       height = 4, width = 6, dpi = "retina")

# ----- do bootstrapping on kinetics
iter = 20
sz = .7
peps = both$pepSeq %>% unique()

pdf("results/termini/kinetics/Transpep-vs-Condens_kinetics_bootstrap.pdf", height = 4, width = 6)
for (i in 1:iter) {
  
  pepCnt = peps[sample(length(peps), ceiling(sz*length(peps)))]
  print(ggplot(both[both$pepSeq %in% pepCnt, ], aes(x = digestTime, y = int_norm, col = mechanism)) +
    geom_smooth(method = "loess")+
    scale_color_manual(values = c("blue","red")) +
    ggtitle("summarised kinetics") +
    ylab("scaled intensity") +
    xlab("time [hrs]"))
  
}

dev.off()
