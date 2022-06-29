### protein PCPS ###
# description:  density of SCS-and PSP-P1 for amino acid groups
# input:        all quantitative data sets:
#               qiSPI: WT-Mut, EGFR
#               aSPIre: TSN5 inhibitor, proteins
# output:       cleavage and splicing strength for different aa groups/residues
# author:       HPR

library(dplyr)
library(stringr)
library(seqinr)
library(Peptides)
library(MASS)
library(gplots)
library(RColorBrewer)

source("src/invitroSPI_utils.R")
source("src/SCS+PSP-P1.R")


### INPUT ###
# proteins
load("../../proteinsPCPS/new/data/aSPIre.RData")
proteinPCPS = Kinetics

# TSN% inhibitor data set
TSN5inhibitor = read.csv("../1_inhibitor_kinetics/invitroSPI+aSPIre/aSPIre_manual/results/TSN5inhibitor/finalKinetics.csv", stringsAsFactors = F)

# EGFR
EGFR = read.csv("/Volumes/DATA16040/DATA/ProteasomeDB/EGFR_kinetics/OUTPUT/KineticsDB.csv", stringsAsFactors = F)

# WT-Mut
fs_wtmut = list.files("/Volumes/DATA16040/DATA/ProteasomeDB/WT-Mut/qiSPI/OUTPUT/", pattern = "finalKinetics.csv", full.names = T, recursive = T)
WT_Mut = lapply(fs_wtmut, read.csv, stringsAsFactors=F)
WT_Mut = plyr::ldply(WT_Mut) %>% as.data.frame()


### MAIN PART ###
# ----- data formatting -----
# proteins
proteinPCPS = proteinPCPS %>%
  tidyr::separate_rows(digestTimes, intensities, sep=";") %>%
  rename(digestTime = digestTimes,
         intensity = intensities) %>%
  mutate(intensity = as.numeric(intensity),
         digestTime = as.numeric(digestTime)) %>%
  filter(digestTime %in% c(3,4)) %>%
  disentangleMultimappers.Type() %>%
  ILredundancy() %>%
  mutate(tool = "aSPIre")

# TSN5: only intensity of no inhibitor condition
TSN5inhibitor = TSN5inhibitor %>%
  tidyr::separate_rows(digestTimes, intensities, sep=";") %>%
  filter(!is.na(substrateID)) %>%
  filter(biological_replicate %in% c("noI_bio1","noI_bio2")) %>%
  rename(digestTime = digestTimes,
         intensity = intensities) %>%
  mutate(intensity = as.numeric(intensity),
         digestTime = as.numeric(digestTime)) %>%
  filter(digestTime == 4) %>%
  disentangleMultimappers.Type() %>%
  ILredundancy() %>%
  mutate(tool = "aSPIre")

# EGFR: only standard proteasome
EGFR = EGFR %>%
  filter(grepl("standard", substrateID)) %>%
  tidyr::separate_rows(digestTimes, intensities, sep=";") %>%
  rename(digestTime = digestTimes,
         intensity = intensities) %>%
  mutate(intensity = as.numeric(intensity),
         digestTime = as.numeric(digestTime)) %>%
  mutate(digestTime = 4) %>%
  disentangleMultimappers.Type() %>%
  ILredundancy() %>%
  mutate(tool = "qiSPI")

# WT/Mut
WT_Mut = WT_Mut %>%
  rename(intensity = tp_4) %>%
  mutate(intensity = as.numeric(intensity)) %>%
  mutate(digestTime = 4) %>%
  disentangleMultimappers.Type() %>%
  ILredundancy() %>%
  mutate(tool = "qiSPI")


# ----- normalise intensities and merge table -----
# merge all
DB = rbind(proteinPCPS %>% dplyr::select(substrateID, substrateSeq, pepSeq, digestTime, intensity, biological_replicate, productType, spliceType, positions, tool),
           TSN5inhibitor %>% dplyr::select(substrateID, substrateSeq, pepSeq, digestTime, intensity, biological_replicate, productType, spliceType, positions, tool),
           EGFR %>% dplyr::select(substrateID, substrateSeq, pepSeq, digestTime, intensity, biological_replicate, productType, spliceType, positions, tool),
           WT_Mut %>% dplyr::select(substrateID, substrateSeq, pepSeq, digestTime, intensity, biological_replicate, productType, spliceType, positions, tool)) %>%
  as.data.frame()

# intensities as mean over biological replicates
DB = DB %>%
  group_by(substrateID, substrateSeq, pepSeq, digestTime, productType, spliceType, positions, tool) %>%
  summarise(intensity = mean(intensity)) %>%
  mutate(biological_replicate = 1)

DB$intensity[DB$tool == "aSPIre"] = (DB$intensity[DB$tool == "aSPIre"] - min(DB$intensity[DB$tool == "aSPIre"])) / (max(DB$intensity[DB$tool == "aSPIre"]) - min(DB$intensity[DB$tool == "aSPIre"]))
DB$intensity[DB$tool == "qiSPI"] = (DB$intensity[DB$tool == "qiSPI"] - min(DB$intensity[DB$tool == "qiSPI"])) / (max(DB$intensity[DB$tool == "qiSPI"]) - min(DB$intensity[DB$tool == "qiSPI"]))

table(DB$spliceType)

################################################################################
AAchar = c("P","G","A","V","L","M","F","Y","W","H","R","K","D","E","N","Q","S","T","C")

getPlots = function(target) {
  
  # ----- retrieve SCS/PSP-P1 ----
  # SCS-and-PSP
  out = SCSandPSP_allSubs(DB, target)
  data = plyr::ldply(out) %>%
    as.data.frame()
  data$target = target
  
  # assign amino acid groups
  data$aa = apply(data,1,function(x){
    substr(DB$substrateSeq[DB$substrateID == x[".id"]][1], start = as.numeric(x["residue"]), stop = as.numeric(x["residue"]))
  })
  
  aagroups = sapply(data$aa, function(x){
    aaComp(x) %>% unlist() %>% t()
  }) %>% 
    t() %>%
    as.data.frame()
  
  aagroups = aagroups[,c(1:9)]
  names(aagroups) = rownames(aaComp("A")[[1]])
  
  # get C and P as special cases
  CP = sapply(data$aa, function(x){
     if (x %in% c("C","P")) 1 else 0
  })
  
  # get all amino acids
  singleaas = sapply(AAchar, function(x){
    sapply(data$aa, function(y){
      if (y == x) 1 else 0
    })
  })
  
  dat = cbind(data, aagroups, CP, singleaas)
  
  # ---- density estimation and plotting -----
  rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
  r <- rf(100)
  
  # plot only SCS- and PSP-P1
  suppressWarnings(dir.create("results/SCS+PSP/"))
  
  cairo_ps(filename = paste0("results/SCS+PSP/SCS+PSP_",target,"_scatter.ps"),
           onefile = T, fallback_resolution = 600, width = 6, height = 6)
  plot(x = dat$scs_mean, y = dat$psp_mean,
       pch = 16,
       xlab = "cleavage strength (%)", ylab = "splicing strength (%)",
       xlim = c(0,50), ylim = c(0,50))
  dev.off()
  
  # plot for amino acid groups
  lim = 100
  
  cairo_ps(filename = paste0("results/SCS+PSP/SCS+PSP_",target,".ps"),
           onefile = T, fallback_resolution = 600, width = 20, height = 24)
  par(mfrow = c(6,5))
  for (j in (ncol(data)+1):ncol(dat)) {
    
    cnt = dat[dat[,j] == 1, ]
    
    dens = try(kde2d(cnt$scs_mean,cnt$psp_mean,n=40,lims=c(0,lim,0,lim)))
    
    if(!is.null(names(dens))) {
      image(log10(dens$z),col=r,
            xlab="cleavage-strength (%)",ylab="splicing strength (%)",main=names(dat)[j],
            sub = paste(cnt$aa %>% unique() %>% sort(), collapse = ", "),
            axes = F)
      axis(1, at = seq(0,1,0.1), labels = seq(0,lim,lim/10))
      axis(2, at = seq(0,1,0.1), labels = seq(0,lim,lim/10))
    } else {
      plot.new()
    }
    
  }
  dev.off()
  
}


################################################################################

### OUTPUT ###
targets = c("P4", "P3", "P2", "P1", "P1_", "P2_", "P3_", "P4_")
sapply(targets, getPlots)

getPlots(target = "P1")
getPlots(target = "P1_")
getPlots(target = "P2")


################################################################################
# get amino acid composition of substrates
AA = c("A","D","E","F","G","H","K","L","N","P","Q","R","S","T","V","W","Y", "M", "C")
AAchar = c("P","G","A","V","I","L","M","F","Y","W","H","R","K","D","E","N","Q","S","T","C")

same_pos = function(tbl) {
  
  if (! all(AA %in% names(tbl))) {
    k = which(! AA %in% names(tbl))
    
    for (i in 1:length(k)) {
      tbl = append(tbl, 0, k[i]-1)
    }
    names(tbl) = AA
  }
  
  return(tbl)
}

substrates = DB %>%
  group_by(substrateID, substrateSeq) %>%
  summarise()

AAcomp = sapply(substrates$substrateSeq,function(x){
  ((strsplit(x,"") %>% unlist() %>% table()) / sum(strsplit(x,"") %>% unlist() %>% table())) %>% same_pos()
}) %>%
  t() %>%
  as.matrix()


image(AAcomp,col=r,
      axes = T)
axis(1, at = seq(0,1,0.1), labels = seq(0,lim,lim/10))
axis(2, at = seq(0,1,0.1), labels = seq(0,lim,lim/10))

