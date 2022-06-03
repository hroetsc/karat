### karat projetc - PCPS mechanism ###
# description:  plot kinetics for each peptide carrying the same SR1 / the same SR2
# input:        quantitative data set: EGFR, WT sequences of WT/Mut
# output:       kinetics of peptides carrying the same SRs
# author:       HPR

library(dplyr)
library(stringr)
library(ggplot2)
library(uwot)
library(twosamples)
library(dgof)
library(RColorBrewer)
library(wesanderson)
options(dplyr.summarise.inform = FALSE)

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

Quant = Kinetics %>%
  tidyr::separate_rows(digestTimes, intensities, sep=";") %>%
  rename(digestTime = digestTimes,
         intensity = intensities) %>%
  mutate(intensity = as.numeric(intensity),
         digestTime = as.numeric(digestTime)) %>%
  ILredundancy() %>%
  disentangleMultimappers.Type() %>%
  disentangleMultimappers.SRlen() %>%
  disentangleMultimappers.IVSeqLen() %>%
  removeMultimappers.SRlen() %>%
  removeMultimappers.AA()

extractSRs = function(DB) {
  
  pos = str_split_fixed(DB$positions, "_", Inf)[,c(1:4)]
  pos = apply(pos,2,as.numeric)
  
  DB$sr1 = substr(DB$substrateSeq, pos[,1], pos[,2])
  DB$sr2 = substr(DB$substrateSeq, pos[,3], pos[,4])
  
  return(DB)
}

Quant = extractSRs(Quant)
QuantPSP = Quant[Quant$productType == "PSP",]

# ----- plotting function for each SR -----
requ = c("pepSeq","spliceType","positions","biological_replicate")


plotKineticsEachSR = function(Qtable, SRcol, outfile, meanBio=F, earlyOnly=T, sortByInt=T) {
  
  # ----- preprocessing -----
  if (!all(requ %in% names(Qtable))) {
    print("!!! plotting requires the following information:")
    print(requ)
  }
  
  Qtable = disentangleMultimappers.Type(Qtable)
  Qtable$spliceType[is.na(Qtable$spliceType)] = "PCP"
  
  if ("int_pasted" %in% names(Qtable)) {
    Qtable = Qtable %>%
      tidyr::separate_rows(int_pasted, times_pasted, sep = ";") %>%
      rename(intensity = int_pasted,
             digestTime = times_pasted) %>%
      as.data.frame()
  }
  
  # ----- sort by intensity -----
  if (sortByInt) {
    x = Qtable %>%
      filter(digestTime == max(digestTime)) %>%
      arrange(desc(intensity))
    pepU = unique(x[,SRcol])
    
  } else {
    pepU = Qtable[,SRcol] %>% unique()
  }
  
  pdf(outfile, height = 16, width = 16)
  par(mfrow = c(4,4))
  
  # ----- iterate peptides -----
  subU = Qtable$substrateID %>% unique()
  for (ii in 1:length(subU)) {
    print(subU[ii])
    
    for (i in 1:length(pepU)) {
      
      xx = which(Qtable[,SRcol] == pepU[i] & Qtable$substrateID == subU[ii])
      if (length(xx) > 0) {
        
        cnt = Qtable[xx,] %>%
          select(-pepSeq) %>%
          rename(pepSeq = SRcol)
        cnt$intensity = as.numeric(cnt$intensity)
        cnt$digestTime = as.numeric(cnt$digestTime)
        
        # ----- mean over technical replicates -----
        if(meanBio) {
          tmp = cnt %>%
            group_by(pepSeq,spliceType,positions,digestTime) %>%
            summarise(mean_int = if (all(intensity == 0) & all(digestTime != 0)) 0 else mean(intensity[intensity!=0 | digestTime == 0], na.rm=T),
                      sd_int = if (all(intensity == 0) & all(digestTime != 0)) 0 else sd(intensity[intensity!=0 | digestTime == 0], na.rm=T)) %>%
            mutate(biological_replicate = 1) %>%
            arrange(digestTime) %>%
            suppressMessages()
        } else {
          tmp = cnt %>%
            select(pepSeq,spliceType,positions,biological_replicate,digestTime,intensity) %>%
            mutate(mean_int = intensity,
                   sd_int = 0) %>%
            arrange(digestTime) %>%
            suppressMessages()
        }
        k = which(!is.na(tmp$mean_int))
        
        if (earlyOnly) {
          tmp = tmp[tmp$digestTime < 20,]
        }
        
        # ----- add colours -----
        nCol = length(unique(tmp$positions))
        Cols = rainbow(nCol)
        
        tmp = tmp %>%
          group_by(biological_replicate,digestTime) %>%
          mutate(technical_replicate = row_number(),
                 combo = paste0(biological_replicate,"-",technical_replicate))
        tmp = data.frame(combo = unique(tmp$combo), col = as.character(Cols)) %>%
          right_join(tmp) %>%
          suppressMessages()
        
        # ----- dot plot -----
        plot(x=tmp$digestTime, y=tmp$mean_int,
             pch= 16, cex = 1.5,
             xlab = "time [hrs]", ylab = "intensity",
             ylim = c(0, max(tmp$mean_int)+100),
             col = tmp$col,
             main = paste0(SRcol, " = ", tmp$pepSeq[1], ", ", subU[ii]),
             sub = paste(unique(tmp$positions), collapse = ","))
        
        # ----- add lines -----
        for (gr in unique(tmp$combo)) {
          lines(x=tmp$digestTime[tmp$combo == gr],
                y=tmp$mean_int[tmp$combo == gr],
                col = tmp$col[tmp$combo==gr][1], lwd=1.5)
          
          # add standard dev
          if(meanBio) {
            arrows(tmp$digestTime[tmp$combo == gr], tmp$mean_int[tmp$combo == gr]-tmp$sd_int[tmp$combo == gr],
                   tmp$digestTime[tmp$combo == gr], tmp$mean_int[tmp$combo == gr]+tmp$sd_int[tmp$combo == gr],
                   length=0.03, angle=90, code=3,
                   lty = "solid", lwd = 1, col = tmp$col[tmp$combo==gr][1]) %>%
              suppressWarnings()
          }
          
        }
        
        
      }
      
      
    }
    
    
  }
  
  dev.off()
  
}


plotKineticsEachSR(Qtable = QuantPSP, SRcol = "sr1", outfile = "results/SR1vsSR2/kinetics_allSR1s.pdf",
                   meanBio = T, earlyOnly = T, sortByInt = F)

plotKineticsEachSR(Qtable = QuantPSP, SRcol = "sr2", outfile = "results/SR1vsSR2/kinetics_allSR2s.pdf",
                   meanBio = T, earlyOnly = T, sortByInt = F)



