### karat projetc - PCPS mechanism ###
# description:  compare specificity of SR1 with SR2
# input:        invitroSPI: Roetschke SciData, EGFR ery, WT substrates --> spliced peptides only
#               random database for qualitative DB 
# output:       sequence entropy of SR1 vs. SR2
# author:       HPR

library(dplyr)
library(stringr)
library(seqinr)
library(ggplot2)
library(RColorBrewer)
source("src/invitroSPI_utils.R")
source("src/plotting_utils.R")

theme_set(theme_classic())
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(100)


### INPUT ###
load("data/invitroSPI.RData")
load("data/randomDB_smart.RData")  # tmp!


### MAIN PART ###
suppressWarnings(dir.create("results/SRspecificity/"))

# ----- preprocessing -----
DB = ProteasomeDB %>%
  filter(productType == "PSP") %>%
  ILredundancy() %>%
  disentangleMultimappers.Type() %>%
  removeMultimappers.Type() %>%
  disentangleMultimappers.AA() %>%
  removeMultimappers.AA() %>%
  disentangleMultimappers.SRlen() %>%
  removeMultimappers.SRlen() %>%
  uniquePeptides()

randomDB = rndDB %>%
  filter(productType == "PSP") %>%
  ILredundancy() %>%
  disentangleMultimappers.Type() %>%
  removeMultimappers.Type() %>%
  disentangleMultimappers.AA() %>%
  removeMultimappers.AA() %>%
  disentangleMultimappers.SRlen() %>%
  removeMultimappers.SRlen() %>%
  uniquePeptides()


# ----- relevant information -----
AA = c("A","D","E","F","G","H","K","L","N","P","Q","R","S","T","V","W","Y", "M", "C")
AAchar = c("P","G","A","V","L","M","F","Y","W","H","R","K","D","E","N","Q","S","T","C")

SR1pos = c("P4"=-3, "P3"=-2, "P2"=-1, "P1"=0,
           "P-1"=1, "P-2"=2, "P-3"=3, "P-4"=4)

SR2pos = c("P-4_"=-4, "P-3_"=-3, "P-2_"=-2, "P-1_"=-1,
           "P1_"=0, "P2_"=1, "P3_"=2, "P4_"=3)

PCPpos = c("P4"=-3, "P3"=-2, "P2"=-1, "P1"=0,
           "P1_"=1, "P2_"=2, "P3_"=3, "P4_"=4)

SRnames = c(names(SR1pos), names(SR2pos))
types = c("cis", "revCis", "trans", "PCP")

# ----- extract amino acids -----
extract_aminoacids = function(tbl = ""){
  
  tbl$spliceType[(tbl$spliceType == "") | (is.na(tbl$spliceType))] = "PCP"
  
  # table with position indices
  pos = str_split_fixed(tbl$positions, coll("_"), Inf) %>% as.data.frame()
  pos = apply(pos, 2, function(x){as.numeric(as.character(x))})
  
  pcp = which(tbl$spliceType == "PCP")
  psp = which(tbl$spliceType != "PCP")
  
  
  # PCPs
  pcpTBL = sapply(PCPpos, function(x){
    substr(tbl$substrateSeq[pcp], start = pos[pcp,2]+x, stop = pos[pcp,2]+x)
  }) %>%
    as.data.frame() %>%
    mutate(spliceType = tbl$spliceType[pcp],
           positions = tbl$positions[pcp],
           pepSeq = tbl$pepSeq[pcp],
           substrateID = tbl$substrateID[pcp],
           substrateSeq = tbl$substrateSeq[pcp])
  
  
  # PSPs
  pspSR1TBL = sapply(SR1pos, function(x){
    substr(tbl$substrateSeq[psp], start = pos[psp,2]+x, stop = pos[psp,2]+x)
  }) %>%
    as.data.frame() %>%
    mutate(spliceType = tbl$spliceType[psp],
           positions = tbl$positions[psp],
           pepSeq = tbl$pepSeq[psp],
           substrateID = tbl$substrateID[psp],
           substrateSeq = tbl$substrateSeq[psp])
  
  
  pspSR2TBL = sapply(SR2pos, function(x){
    substr(tbl$substrateSeq[psp], start = pos[psp,3]+x, stop = pos[psp,3]+x)
  }) %>%
    as.data.frame() %>%
    mutate(spliceType = tbl$spliceType[psp],
           positions = tbl$positions[psp],
           pepSeq = tbl$pepSeq[psp],
           substrateID = tbl$substrateID[psp],
           substrateSeq = tbl$substrateSeq[psp])
  
  # merge all tables
  pspTBL = cbind(pspSR1TBL[,names(SR1pos)], pspSR2TBL)
  
  pcpPlaceholder = matrix("", length(pcp), length(SR2pos)) %>%
    as.data.frame()
  names(pcpPlaceholder) = SRnames[! SRnames %in% names(PCPpos)]
  pcpTBL2 = cbind(cbind(pcpTBL[,names(PCPpos)], pcpPlaceholder)[, c(names(SR1pos), names(SR2pos))],
                  pcpTBL[,which(!names(pcpTBL) %in% names(PCPpos))])
  
  TBL = rbind(pspTBL, pcpTBL2) %>% as.data.frame()
  return(TBL)
}


# ----- get splice-reactants -----

getSRs = function(df) {
  
  pos = str_split_fixed(df$positions, "_", Inf)[,c(1:4)]
  pos = apply(pos,2,as.numeric)
  
  df$sr1 = substr(df$substrateSeq, pos[,1], pos[,2])
  df$sr2 = substr(df$substrateSeq, pos[,3], pos[,4])
  
  return(df)
}

DB = getSRs(DB)
randomDB = getSRs(randomDB)

DBproc = left_join(DB, extract_aminoacids(DB))
randomDBproc = left_join(randomDB, extract_aminoacids(randomDB))

# ----- amino acid frequency + sequence entropy ----

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


getEntropy = function(DB, randomDB) {
  
  getFreq = function(df) {
    
    freq = sapply(c("sr1","sr2"),function(x){
      
      if (x == "sr1") {
        a = substr(df[,x], nchar(df[,x])-2, nchar(df[,x]))
      } else if (x == "sr2") {
        a = substr(df[,x], 1, 3)
      }
      y = table(unlist(strsplit(a,"")))
      z = same_pos(y)
      z = z[match(AAchar,names(z))]
      # print(names(z))
      return(z)
    }) %>%
      as.matrix()
    rownames(freq) = AAchar
    
    Freq = apply(freq,2,function(x){
      x/sum(x)
    })
    return(Freq)
  }
  
  ftrue = getFreq(DB)
  frnd = getFreq(randomDB)
  
  Fs = ftrue[AAchar,]/frnd[AAchar,]
  # Fs = ftrue[AAchar,]
  
  H = apply(Fs,2,function(sr){
    sr = (sr-min(sr, na.rm = T))/(max(sr, na.rm = T) - min(sr, na.rm = T))
    return(-1*sum(sr*log2(sr), na.rm = T))
  })
  
  return(H)
}

# ----- apply to each substrate -----

subs = DBproc$substrateID %>% unique()
entropy = sapply(subs, function(s){
  getEntropy(DBproc[DBproc$substrateID == s, ],
             randomDBproc[randomDBproc$substrateID == s, ])
}) %>%
  t() %>%
  as.data.frame()

Hdf = entropy %>%
  mutate(substrateID = rownames(entropy)) %>%
  tidyr::gather(sr, entropy, -substrateID)

Hdf %>%
  group_by(sr) %>%
  summarise(n = n(),
            mean = mean(entropy),
            median = median(entropy),
            std = sd(entropy)) %>%
  print.data.frame()

# plot
hbox = ggplot(Hdf, aes(x = sr, y = entropy, fill = sr)) +
  geom_boxplot(draw_quantiles = c(0.5)) + 
  scale_fill_manual(values = c("gray","lightblue")) +
  ggtitle("SR sequence entropy") +
  ylab("Shannon Measure of Information [bits]")
hbox

### OUTPUT ###
ggsave("results/SRspecificity/SRspecificity_entropy.png", plot = hbox,
       height = 4, width = 4, dpi = "retina")

save(DBproc, file = "results/SRspecificity/DBproc.RData")
save(randomDBproc, file = "results/SRspecificity/randomDBproc.RData")

