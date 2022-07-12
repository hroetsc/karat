### karat projetc - PCPS mechanismS ###
# description:  density of SCS-and PSP-P1 for amino acid groups
# input:        EGFR+WT quantitative results (aSPIre), quantitative random database
#               protein PCPS quantitative results (aSPIre)
# output:       cleavage and splicing strength for different aa groups/residues
# author:       HPR

library(dplyr)
library(stringr)
library(seqinr)
library(Peptides)
library(MASS)
library(gplots)
library(RColorBrewer)
library(ggplot2)
library(twosamples)
theme_set(theme_classic())

source("src/invitroSPI_utils.R")
source("src/SCS+PSP-P1.R")

rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(100)


### INPUT ###
load("../../proteinsPCPS/new/data/aSPIre.RData")
proteins = Kinetics
load("data/aSPIre.RData")
load("data/randomDB_Quant_aSPIre.RData")

### MAIN PART ###
suppressWarnings(dir.create("results/SCS+PSP/"))

# ----- preprocessing -----
polypeps = Kinetics %>%
  disentangleMultimappers.Type() %>%
  tidyr::separate_rows(digestTimes, intensities, sep=";") %>%
  rename(digestTime = digestTimes,
         intensity = intensities) %>%
  mutate(intensity = as.numeric(intensity),
         digestTime = as.numeric(digestTime)) %>%
  filter(digestTime == 4) %>%
  ILredundancy() %>%
  resolve_multimapper() %>%
  tidyr::separate_rows(positions, sep = ";")

proteins = proteins %>%
  disentangleMultimappers.Type() %>%
  tidyr::separate_rows(digestTimes, intensities, sep=";") %>%
  rename(digestTime = digestTimes,
         intensity = intensities) %>%
  mutate(intensity = as.numeric(intensity),
         digestTime = as.numeric(digestTime)) %>%
  filter(digestTime %in% c(3,4)) %>%
  ILredundancy() %>%
  resolve_multimapper() %>%
  tidyr::separate_rows(positions, sep=";")

randomQuant = randomQuant %>%
  ILredundancy()

################################################################################
AAchar = c("P","G","A","V","L","M","F","Y","W","H","R","K","D","E","N","Q","S","T","C")

getPlots = function(DB, target, suffix, limited = F) {
  
  # ----- retrieve SCS/PSP-P1 ----
  # SCS-and-PSP
  out = SCSandPSP_allSubs(DB = DB, target = target, SR2forSCS = T, meanOverBio = T, Zscale = F, rawVals = F)  # !!!!!
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
  # plot for amino acid groups
  lim = 100
  
  cairo_ps(filename = paste0("results/SCS+PSP/SCS+PSP_",target,suffix,".ps"),
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
  
  # plot scatter plot
  cairo_ps(filename = paste0("results/SCS+PSP/SCS+PSP_scatter_",target,suffix,".ps"),
           onefile = T, fallback_resolution = 600, width = 20, height = 24)
  par(mfrow = c(6,5))
  for (j in (ncol(data)+1):ncol(dat)) {
    cnt = dat[dat[,j] == 1, ]
    
    if(limited) {
      lim = max(cnt$scs_mean,cnt$psp_mean)
    } else {
      lim = 100
    }
    
    plot(x = cnt$scs_mean, y = cnt$psp_mean, pch = 16, cex = .8,
         xlab = "cleavage strength (%)", ylab = "splicing strength (%)",
         xlim = c(0,lim), ylim = c(0,lim),
         main=names(dat)[j],
         sub = paste0(paste(cnt$aa %>% unique() %>% sort(), collapse = ", "), " - n=", nrow(cnt)))
  }
  dev.off()
  
  # ggplot
  allP = list()
  counter = 1
  for (j in (ncol(data)+1):ncol(dat)) {
    cnt = dat[dat[,j] == 1, ]
    
    if(limited) {
      lim = max(cnt$scs_mean,cnt$psp_mean)
    } else {
      lim = 100
    }
    
    # statistical test
    set.seed(1234)
    pval = ad_test(cnt$scs_mean, cnt$psp_mean)[2]
    
    p = ggplot(cnt, aes(x = scs_mean, y = psp_mean)) +
      geom_point() +
      xlab("cleavage strength (%)") +
      ylab("splicing strength (%)") +
      ggtitle(names(dat)[j], subtitle = paste0(paste(cnt$aa %>% unique() %>% sort(), collapse = ", "), " - n=", nrow(cnt), ", p = ", pval)) +
      xlim(0,lim) + ylim(0,lim)
    
    allP[[counter]] = ggExtra::ggMarginal(p, type = "density", size = 8, col = "gray")
    counter = counter+1
  }
  
  ggsave(filename = paste0("results/SCS+PSP/SCS+PSP_scatterDens_",target,suffix,".pdf"), 
         plot = gridExtra::marrangeGrob(allP, nrow=6, ncol=5, byrow = T), 
         width = 20, height = 24, dpi = "retina")
  
  
  names(allP) = names(dat)[(ncol(data)+1):ncol(dat)]
  return(list(allP = allP,
              data = data))
}


################################################################################
nm = intersect(names(proteins), names(polypeps))
ALL = rbind(proteins[,nm], polypeps[,nm])
table(ALL$spliceType)

table(uniquePeptides(ALL)$productType)
unique(ALL$substrateID) %>% length()

# downsample the random databases
x = nrow(Kinetics)*2
randomQuant_down = randomQuant[sample(seq(1,nrow(randomQuant)), x), ] %>%
  mutate(biological_replicate = 1)

### OUTPUT ###
# targets = c("P4", "P3", "P2", "P1", "P1_", "P2_", "P3_", "P4_")
targets = c("P1","P-1","P1_")
sapply(targets, function(t){
  # getPlots(DB = polypeps, target = t, suffix = "_PolypepTrueLim", limited = T)
  getPlots(DB = ALL, target = t, suffix = "_trueLim", limited = T)
  # getPlots(DB = randomQuant_down, target = t, suffix = "_random")
})

P1 = getPlots(DB = ALL, target = "P1", suffix = "_trueLim", limited = T)
Pm1 = getPlots(DB = ALL, target = "P-1", suffix = "_trueLim", limited = T)


# ----- plots for thesis -----
lim = 100
pval = ad_test(P1$data$scs_mean, P1$data$psp_mean)[2]

p = ggplot(P1$data, aes(x = scs_mean, y = psp_mean)) +
  geom_point() +
  xlab("cleavage strength (%)") +
  ylab("splicing strength (%)") +
  ggtitle("", subtitle = paste0("all residues, n=", nrow(P1$data), ", p = ", pval)) +
  xlim(0,lim) + ylim(0,lim)
p = ggExtra::ggMarginal(p, type = "density", size = 8, col = "gray")
p

ggsave(filename = "results/SCS+PSP/_forthesis_P1.png", plot = p, height = 3.5, width = 3.5, dpi = "retina")


pval = ad_test(Pm1$data$scs_mean, Pm1$data$psp_mean)[2]

p = ggplot(Pm1$data, aes(x = scs_mean, y = psp_mean)) +
  geom_point() +
  xlab("cleavage strength (%)") +
  ylab("splicing strength (%)") +
  ggtitle("", subtitle = paste0("all residues, n=", nrow(Pm1$data), ", p = ", pval)) +
  xlim(0,lim) + ylim(0,lim)
p = ggExtra::ggMarginal(p, type = "density", size = 8, col = "gray")
p

ggsave(filename = "results/SCS+PSP/_forthesis_P-1.png", plot = p, height = 3.5, width = 3.5, dpi = "retina")


ggsave(filename = "results/SCS+PSP/_forthesis_P1_P.png", plot = P1$allP$P, height = 3.5, width = 3.5, dpi = "retina")
ggsave(filename = "results/SCS+PSP/_forthesis_P1_G.png", plot = P1$allP$G, height = 3.5, width = 3.5, dpi = "retina")
ggsave(filename = "results/SCS+PSP/_forthesis_P1_M.png", plot = P1$allP$M, height = 3.5, width = 3.5, dpi = "retina")
ggsave(filename = "results/SCS+PSP/_forthesis_P1_F.png", plot = P1$allP$`F`, height = 3.5, width = 3.5, dpi = "retina")
ggsave(filename = "results/SCS+PSP/_forthesis_P1_R.png", plot = P1$allP$R, height = 3.5, width = 3.5, dpi = "retina")
ggsave(filename = "results/SCS+PSP/_forthesis_P1_H.png", plot = P1$allP$H, height = 3.5, width = 3.5, dpi = "retina")


ggsave(filename = "results/SCS+PSP/_forthesis_P-1_P.png", plot = Pm1$allP$P, height = 3.5, width = 3.5, dpi = "retina")
ggsave(filename = "results/SCS+PSP/_forthesis_P-1_G.png", plot = Pm1$allP$G, height = 3.5, width = 3.5, dpi = "retina")


################################################################################
# get amino acid composition of substrates
AA = c("A","D","E","F","G","H","K","L","N","P","Q","R","S","T","V","W","Y", "M", "C")
AAchar = c("P","G","A","V","I","L","M","F","Y","W","H","R","K","D","E","N","Q","S","T","C")
DB = ALL

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

png("results/SCS+PSP/AAcomposition.png", height = 8, width = 10, units = "in", res = 300)
image(AAcomp %>% t(),col=r,
      axes = F)
axis(1, at = seq(0,1,1/(length(AA)-1)), labels = AA)
axis(2, at = seq(0,1,1/(nrow(substrates)-1)), labels = substrates$substrateID)
dev.off()
