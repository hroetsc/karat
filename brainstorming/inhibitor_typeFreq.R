### PCPS brainstorming ###
# description:  determine product type frequency in inhibitor data set
# input:        inhibitor kinetics of TSN5 analysed with invitroSPI
# output:       product type frequencies
# author:       HPR

library(dplyr)
library(stringr)
library(ggplot2)
source("../../proteinsPCPS/new/src/invitroSPI_utils.R")
source("../../proteinsPCPS/new/src/plotting_utils.R")

theme_set(theme_classic())

### INPUT ###
sample_list = read.csv("../1_inhibitor_kinetics/invitroSPI/INPUT/sample_list.csv", stringsAsFactors = F)
ProteasomeDB = read.csv("../1_inhibitor_kinetics/invitroSPI/OUTPUT/TSN5/ProteasomeDB.csv", stringsAsFactors = F)

### MAIN PART ###
# ----- preprocessing -----
# join with sample list to extract information
sample_list$digestTime = as.numeric(sample_list$digestTime)
ProteasomeDB$filename = str_extract_all(ProteasomeDB$runID, "F0[:digit:]+.csv", simplify = T)
MASTER = left_join(ProteasomeDB, sample_list)

# get info about inhibitors
MASTER$sample = str_extract_all(MASTER$MSfile, "(?<=J_Liepe_)[:alpha:]{2}5", simplify = F)
MASTER = MASTER %>%
  mutate(condition = ifelse(grepl("A",sample), "no_inhibitor", NA),
         condition = ifelse(grepl("B",sample), "b5_inhibited", condition),
         condition = ifelse(grepl("C",sample), "b2_inhibited", condition),
         condition = ifelse(grepl("D",sample), "b1_inhibited", condition))

# normal preprocessing stuff
MASTER = MASTER %>%
  ILredundancy() %>%
  disentangleMultimappers.Type() %>%
  disentangleMultimappers.SRlen() %>%
  disentangleMultimappers.IVSeqLen() %>%
  remSynthErrors() %>%
  uniquePeptides() %>%
  DBcosmetics()

# ----- get product type frequency -----
abund = MASTER %>%
  group_by(condition, spliceType) %>%
  summarise(n = n()) %>%
  ungroup() %>% group_by(condition) %>%
  mutate(freq = n/sum(n),
         height = cumsum(freq))

abund$spliceType = factor(abund$spliceType, levels = c("PCP","cis","revCis","trans"))
freqplot = ggplot(abund, aes(x = condition, y = freq, fill = spliceType)) +
  geom_bar(position = "stack", stat = "identity") +
  # geom_text(aes(label = n), size = 3, y = abund$height) +
  scale_fill_manual("product type", values = c(plottingCols["PCP"], plottingCols["cis"],
                                               plottingCols["revCis"], plottingCols["trans"])) +
  ylab("frequency") +
  ggtitle("product type frequency upon subunit inhibition") +
  theme(axis.text.x = element_text(angle = 90))

### OUTPUT ###
ggsave(filename = "results/inhibitor_productTypeFreq.png", plot = freqplot,
       height = 4, width = 6, dpi = "retina")


# ----- get lengths -----
pos = str_split_fixed(MASTER$positions,coll("_"),Inf)
pos = apply(pos,2,as.numeric)

LenDB = MASTER %>%
  mutate(pepLen = nchar(pepSeq),
         sr1Len = pos[,2]-pos[,1]+1)

LenDB$condition = factor(LenDB$condition, levels = c("no_inhibitor","b5_inhibited","b2_inhibited","b1_inhibited"))

peplen = ggplot(LenDB, aes(x=spliceType, y=pepLen, fill = condition)) +
  geom_boxplot(alpha = .8)+
  scale_fill_manual("condition", values = c("gray","lightpink","lightblue","seagreen")) +
  ylab("length (aa residues)") +
  ggtitle("peptide length") +
  theme(axis.text.x = element_text(angle = 90))

sr1len = ggplot(LenDB[LenDB$spliceType != "PCP", ], aes(x=spliceType, y=sr1Len, fill = condition)) +
  geom_boxplot(alpha = .8)+
  scale_fill_manual("condition", values = c("gray","lightpink","lightblue","seagreen")) +
  ylab("length (aa residues)") +
  ggtitle("SR1 length") +
  theme(axis.text.x = element_text(angle = 90))

lenplot = gridExtra::grid.arrange(peplen, sr1len, nrow=2)

### OUTPUT ###
ggsave(filename = "results/inhibitor_LenDist.png", plot = lenplot,
       height = 6, width = 8, dpi = "retina")


