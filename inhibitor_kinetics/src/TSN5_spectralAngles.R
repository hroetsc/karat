### protein PCPS ###
# description:  plot spectral angles for TSN5 (no inhibitor)
# input:        spectral angle predictions (from JAC)
# output:       spectral angle distributions for peptides identified with invitroSPI
# author:       HPR


library(dplyr)
library(stringr)
library(ggplot2)
source("../../proteinsPCPS/new/src/invitroSPI_utils.R")
source("../../proteinsPCPS/new/src/plotting_utils.R")

theme_set(theme_classic())

### INPUT ###
TSN5 = read.csv("data/spectralAngles/tsn5_Final.csv", stringsAsFactors = F)


### MAIN PART ###
# no cysteines
# shorter than 13 amino acids
filtering = function(DB) {
  DB = DB[which(nchar(DB$pepSeq) < 13),]
  
  k = which(grepl("C",DB$pepSeq))
  if (length(k) > 0) {
    DB = DB[-k, ]
  }
  
  return(DB)
}

TSN5 = filtering(TSN5) %>%
  ILredundancy() %>%
  disentangleMultimappers.Type()

TSN5$spliceType[TSN5$spliceType == ""] = "PCP"
TSN5$spliceType = factor(TSN5$spliceType, levels = c("PCP","cis","revCis","trans"))

# plotting
pl = ggplot(TSN5, aes(x = spliceType, y = spectralAngle, fill = spliceType)) +
  geom_violin() +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", 
               width = .2,
               position = position_dodge(width = .2)) +
  scale_fill_manual(values =c(plottingCols["PCP"],plottingCols["cis"],plottingCols["revCis"],plottingCols["trans"])) +
  xlab("product type") +
  xlab("peptide products") +
  ylab("Prosit spectral angle") +
  ggtitle("TSN5 - no inhibitor")

pl

ggsave(filename = "data/spectralAngles/TSN5_spectralAngles.png",
       plot = pl, dpi = "retina")


# stats
TSN5 %>%
  group_by(spliceType) %>%
  summarise(n = n(),
            med = median(nchar(pepSeq)))



