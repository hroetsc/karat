### karat projetc - PCPS mechanism ###
# description:  kinetics of peptides with short SR1s / intervening seq lengths
# input:        quantitative data set: EGFR, WT sequences of WT/Mut
# output:       summarised kinetics
# author:       HPR

library(dplyr)
library(stringr)
library(ggplot2)
library(twosamples)
library(dgof)
source("src/invitroSPI_utils.R")
source("src/aSPIre_plotting.R")
source("src/plotting_utils.R")

theme_set(theme_classic())

### INPUT ###
load("data/aSPIre.RData")

### MAIN PART ###
suppressWarnings(dir.create("results/termini/kinetics/"))

# ----- preprocessing -----

polypeps = Kinetics %>%
  ILredundancy() %>%
  disentangleMultimappers.Type() %>%
  disentangleMultimappers.SRlen() %>%
  disentangleMultimappers.IVSeqLen() %>%
  tidyr::separate_rows(digestTimes, intensities, sep=";") %>%
  rename(digestTime = digestTimes,
         intensity = intensities) %>%
  mutate(intensity = as.numeric(intensity),
         digestTime = as.numeric(digestTime))

# get lengths
pos = str_split_fixed(polypeps$positions,coll("_"),Inf)[,c(1:4)]
pos = apply(pos,2,as.numeric) %>% as.data.frame()
names(pos) = c("pos1","pos2","pos3","pos4")

DB = cbind(polypeps, pos) %>%
  mutate(sr1len = pos2-pos1+1,
         sr2len = pos4-pos3+1,
         ivlen = abs(pos3-pos2)-1)


# ----- normalise kinetics -----
# mean over bio rep
DB = DB %>%
  group_by(substrateID,pepSeq,digestTime) %>%
  mutate(intensity = mean(intensity)) %>%
  filter(biological_replicate == 1)

# scale between 0 and 1
DB = DB %>%
  group_by(substrateID,pepSeq) %>%
  mutate(int_norm = (intensity - min(intensity))/(max(intensity) - min(intensity)))

# ----- short vs long SR1 -----
DB = DB %>%
  mutate(sr1_cat = ifelse(sr1len <= 3, "shortSR1", "longSR1"))

table(DB$spliceType, DB$sr1_cat)

srkin = ggplot(DB[DB$spliceType != "PCP",], aes(x = digestTime, y = int_norm, col = sr1_cat)) +
  geom_smooth(method = "loess")+
  scale_color_manual(values = c("black","darkgray")) +
  ggtitle("summarised kinetics") +
  ylab("scaled intensity") +
  xlab("time [hrs]") +
  facet_wrap(~spliceType, scales = "free")

ggsave(filename = "results/length/qual_kinetics_shortSR1.png", plot = srkin,
       height = 8, width = 12, dpi = "retina")


# ----- short vs long intervening sequence -----
DB = DB %>%
  mutate(iv_cat = ifelse(ivlen <= 1, "shortIVSeq", "longIVSeq"))

table(DB$spliceType, DB$iv_cat)

ivkin = ggplot(DB[DB$spliceType != "PCP",], aes(x = digestTime, y = int_norm, col = iv_cat)) +
  geom_smooth(method = "loess")+
  scale_color_manual(values = c("black","darkgray")) +
  ggtitle("summarised kinetics") +
  ylab("scaled intensity") +
  xlab("time [hrs]") +
  facet_wrap(~spliceType, scales = "free")

ggsave(filename = "results/length/qual_kinetics_shortIVSeq.png", plot = ivkin,
       height = 8, width = 12, dpi = "retina")



### OUTPUT ###


