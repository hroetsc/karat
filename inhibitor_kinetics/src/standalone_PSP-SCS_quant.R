### INHIBITOR KINETICS ###
# description:  substrate-specific cleavage and splicing strength
#               as defined in Mishto et al. FI 2019
# input:        quantification results (finalKinetics.csv), selected peptides
# output:       PSP-P1 and SCS-P1
#               a) for each position in a substrate
#               b) PSP-P1 vs. SCS-P1
# author:       HR

library(dplyr)
library(seqinr)
library(stringr)
source("../brainstorming/src/invitroSPI_utils.R")


### INPUT ###
# final Kinetics
DB = read.csv("qiSPI/OUTPUT/TSN5_0+4/finalKinetics.csv", stringsAsFactors = F)
# selected peptides
selectedPeps = read.csv("results/b5/selectedPeps_intensities.csv", stringsAsFactors = F)
load("qiSPI/OUTPUT/TSN5_0+4/filteredResults.RData")

### MAIN PART ###
# ----- preprocessing -----
# taking into account all possible positions of a multi-mapping peptide
# with a unique product type
# considering all biological replicates
nm = c("XA_r1", "XA_r2", "YA_r1", "YA_r2",
       "XB_r1", "XB_r2", "YB_r1", "YB_r2",
       "XC_r1", "XC_r2", "YC_r1", "YC_r2",
       "XD_r1", "XD_r2", "YD_r1", "YD_r2")
names(filteredResults) = nm

res = plyr::ldply(filteredResults) %>%
  as.data.frame()
names(res)[1] = "sample"

# get conditions
cond = data.frame(sample = nm,
                  condition = c("noInh_bio1_tech1", "noInh_bio1_tech2", "noInh_bio2_tech1", "noInh_bio_tech2",
                                "b5_bio1_tech1", "b5_bio1_tech2", "b5_bio2_tech1", "b5_bio2_tech2",
                                "b2_bio1_tech1", "b2_bio1_tech2", "b2_bio2_tech1", "b2_bio2_tech2",
                                "b1_bio1_tech1", "b1_bio1_tech2", "b1_bio2_tech1", "b1_bio2_tech2"),
                  biological_replicate = rep(seq(1,8), each=2))
res = left_join(res, cond) %>% as.data.frame()
names(res)[2] = "pepSeq"

DB = left_join(DB, res) %>%
  as.data.frame()

# only selected peptides
DB = DB[DB$pepSeq %in% selectedPeps$pepSeq, ]
# only bioReps that are not b5
DB = DB[DB$biological_replicate %in% c(1,2,7,8,5,6), ]

DB = DB %>%
  ILredundancy() %>%
  disentangleMultimappers.Type() %>%
  removeMultimappers.Type() %>%
  as.data.frame()


# ----- resolve multi-mappers -----
# assign weight to position multi-mappers
resolve_multimapper = function(ProteasomeDB) {
  
  k = which(str_detect(ProteasomeDB$positions, coll(";")))
  if (length(k) > 0) {
    
    DB_mm = ProteasomeDB[k, ]
    DB = ProteasomeDB[-k, ]
    
    for (r in 1:nrow(DB_mm)) {
      cnt = DB_mm[r, ] %>%
        tidyr::separate_rows(positions, sep=";") %>%
        as.data.frame()
      
      cnt$intensity = cnt$intensity / nrow(cnt)
      
      DB = rbind(DB, cnt)
    }
    
    return(as.data.frame(DB))
    
  } else {
    return(ProteasomeDB)
  }
  
}
  
# ----- compute SCS-P1 and PSP-P1 -----
SCS_and_PSP = function(ProteasomeDB) {
  
  ProteasomeDB = resolve_multimapper(ProteasomeDB)
  ProteasomeDB$productType = toupper(ProteasomeDB$productType)
  
  pcpidx = which(ProteasomeDB$productType == "PCP")
  pspidx = which(ProteasomeDB$productType == "PSP")
  
  res = str_split(ProteasomeDB$positions, "_", simplify = T)[,2] %>%
    as.numeric()
  
  d = data.frame(residue = c(1:nchar(ProteasomeDB$substrateSeq[1])),
                 scs_mean = 0,
                 scs_sd = 0,
                 scs_n = 0,
                 psp_mean = 0,
                 psp_sd = 0,
                 psp_n = 0)
  
  for (i in 1:nrow(d)) {
    
    cntscs = ProteasomeDB[pcpidx,][res[pcpidx] == d$residue[i],] %>%
      dplyr::group_by(biological_replicate) %>%
      dplyr::summarise(int = sum(intensity))
    
    d$scs_mean[i] = paste(cntscs$int, collapse = "_")
    # d$scs_mean[i] = paste(log2(cntscs$int+1), collapse = "_")
    d$scs_n[i] = ProteasomeDB$pepSeq[pcpidx][res[pcpidx] == d$residue[i]] %>%
      unique() %>%
      length()
    
    cntpsp = ProteasomeDB[pspidx,][res[pspidx] == d$residue[i],] %>%
      dplyr::group_by(biological_replicate) %>%
      dplyr::summarise(int = sum(intensity))
    
    d$psp_mean[i] = paste(cntpsp$int, collapse = "_")
    # d$psp_mean[i] = paste(log2(cntpsp$int+1), collapse = "_")
    d$psp_n[i] = ProteasomeDB$pepSeq[pspidx][res[pspidx] == d$residue[i]] %>%
      unique() %>%
      length()
  }
  
  d[is.na(d)] = 0
  d[d == ""] = 0
  
  # remove substrate's C terminus from the SCS
  d$scs_mean[nchar(ProteasomeDB$substrateSeq[1])] = 0
  d$scs_sd[nchar(ProteasomeDB$substrateSeq[1])] = 0
  
  # normalise by sum for each biological replicate
  scs = apply(str_split(d$scs_mean, pattern = "_", simplify = T), 2, as.numeric) %>%
    as.data.frame()
  scs[is.na(scs)] = 0
  
  psp = apply(str_split(d$psp_mean, pattern = "_", simplify = T), 2, as.numeric) %>%
    as.data.frame()
  psp[is.na(psp)] = 0
  
  # tmp!!!
  # normalise by percentage of cleavage/splicing for each residue
  # s = scs+psp
  # scs = scs/s
  # psp = psp/s
  
  scs = sweep(scs, 2, colSums(scs), FUN = "/")
  d$scs_mean = rowMeans(scs) * 100
  d$scs_sd = apply(scs,1,sd) * 100
  
  psp = sweep(psp, 2, colSums(psp), FUN = "/")
  d$psp_mean = rowMeans(psp) * 100
  d$psp_sd = apply(psp,1,sd) * 100
  
  d[is.na(d)] = 0
  
  return(d)
}

# ----- plotting -----
plotPSPandSCS = function(ProteasomeDB, tp, name) {
  
  out = SCS_and_PSP(ProteasomeDB)
  yl = max(out$scs_mean+max(out$scs_sd), out$psp_mean+max(out$psp_mean)) %>% ceiling()
  
  p = data.frame(residue = rep(out$residue, 2),
                 n = c(out$scs_n, out$psp_n),
                 p1 = c(out$scs_mean, -1*out$psp_mean),
                 stdev = c(out$scs_sd, out$psp_sd),
                 col = c(rep(plottingCols["PCP"], nrow(out)),
                         rep(plottingCols["PSP"], nrow(out))))
  
  pdf(paste0("results/quantPSPandSCS1_", name, "_", tp, "hrs.pdf"),
      height = 6, width = 12)
  plot(p1~residue,
       data=p,
       type = "h",
       lwd = 3,
       col = col,
       main = paste0(name, ": ", tp, " hrs"),
       ylim = c(-yl, yl),
       xlab = "substrate residue forming SCS-/PSP-P1",
       ylab = "substrate-specific cleavage/splicing strength (%)",
       sub = "mean +- S.D.",
       axes=F)
  arrows(p$residue, p$p1-p$stdev, p$residue, p$p1+p$stdev,
         length=0.03, angle=90, code=3,
         lty = "solid", lwd = 1, col = p$col)
  xlabel = c("0",
             paste(ProteasomeDB$substrateSeq[1] %>% strsplit("") %>% unlist(), seq(1, nchar(ProteasomeDB$substrateSeq[1])), sep = ""))
  # positions
  text(x = seq(0, nchar(ProteasomeDB$substrateSeq[1])), par("usr")[3]-.3, labels = xlabel,
       srt=90, xpd=T, cex=.5)
  # number of peptides
  text(x = seq(0, nchar(ProteasomeDB$substrateSeq[1])), par("usr")[3]-5,
       labels = c("PCP:",out$scs_n), srt=0, xpd=T, cex=.5)
  text(x = seq(0, nchar(ProteasomeDB$substrateSeq[1])), par("usr")[3]-10,
       labels = c("PSP:",out$psp_n), srt=0, xpd=T, cex=.5)
  axis(2)
  dev.off()
  
  pdf(paste0("results/quantPSPandSCS2_", name, "_", tp, "hrs.pdf"),
      height = 6, width = 6)
  plot(psp_mean~scs_mean, data = out,
       main = paste0(name, ": ", tp, " hrs"),
       pch = 16,
       xlim = c(0, yl),
       ylim = c(0, yl),
       cex = .8,
       ylab = "mean PSP-P1 (%)",xlab = "mean SCS-P1 (%)")
  dev.off()
  
}

# ----- execute -----
# intensities at 4 hrs detected in all replicates except b5 inhibition
# mean over technical+biological replicates
uniquePeps = unique(DB$pepSeq)
for (u in uniquePeps) {
  k = which(DB$pepSeq == u)
  DB$biological_replicate[k] = seq(1,12)
}

DB$intensity = as.numeric(DB$`4`)
plotPSPandSCS(ProteasomeDB = DB[DB$biological_replicate %in% c(1,2),],
              tp = 4, name = "TSN5abs")
# plotPSPandSCS(ProteasomeDB = DB, tp = 4, name = "TSN5perc")


# DB$intensity = log(DB$`4`+1)
# plotPSPandSCS(ProteasomeDB = DB, tp = 4, name = "TSN5log")


# ----- plot candidates -----
cdts = data.frame(V13 = c(10,14),
         V13_2 = c(10,15),
         D26 = c(26,28),
         D26_2 = c(24,28),
         F15 = c(14,17),
         F7 = c(5,8),
         I9 = c(5,10),
         I9_2 = c(9,12),
         S20 = c(17,22),
         S20_2 = c(17,24),
         N29 = c(28,30),
         D11 = c(11,12)) %>% t() %>% as.data.frame()

plotCandidates = function(ProteasomeDB) {
  
  out = SCS_and_PSP(ProteasomeDB)
  yl = max(out$scs_mean+max(out$scs_sd), out$psp_mean+max(out$psp_mean)) %>% ceiling()
  
  p = data.frame(residue = c(out$residue, out$residue+.2),
                 n = c(out$scs_n, out$psp_n),
                 p1 = c(out$scs_mean, out$psp_mean),
                 stdev = c(out$scs_sd, out$psp_sd),
                 col = c(rep(plottingCols["PCP"], nrow(out)),
                         rep(plottingCols["PSP"], nrow(out))))
  
  plot(p1~residue,
       data=p,
       type = "h",
       lwd = 3,
       col = col,
       ylim = c(-10,yl),
       xlab = "substrate residue forming SCS-/PSP-P1",
       ylab = "substrate-specific cleavage/splicing strength (%)",
       sub = "mean +- S.D.",
       axes=F)
  arrows(p$residue, p$p1-p$stdev, p$residue, p$p1+p$stdev,
         length=0.03, angle=90, code=3,
         lty = "solid", lwd = 1, col = p$col)
  xlabel = c("0",
             paste(ProteasomeDB$substrateSeq[1] %>% strsplit("") %>% unlist(), seq(1, nchar(ProteasomeDB$substrateSeq[1])), sep = ""))
  # positions
  text(x = seq(0, nchar(ProteasomeDB$substrateSeq[1])), par("usr")[3]-.3, labels = xlabel,
       srt=90, xpd=T, cex=.5)
  # number of peptides
  text(x = seq(0, nchar(ProteasomeDB$substrateSeq[1])), par("usr")[3]-5,
       labels = c("PCP:",out$scs_n), srt=0, xpd=T, cex=.5)
  text(x = seq(0, nchar(ProteasomeDB$substrateSeq[1])), par("usr")[3]-8,
       labels = c("PSP:",out$psp_n), srt=0, xpd=T, cex=.5)
  
  y = -1
  cnt.intercept = -1
  for (i in 1:nrow(cdts)) {
    segments(x0 = cdts[i,1], y0 = cnt.intercept, x1 = cdts[i,2], y1 = cnt.intercept,
             lwd = .5)
    points(x=str_extract(rownames(cdts)[i], "[:digit:]+"), y=cnt.intercept,
           pch=16, cex=.5)
    cnt.intercept = cnt.intercept-.8
  }
  
  axis(2)
  
}

pdf(paste0("results/b5/candidates.pdf"), height = 6, width = 12)
plotCandidates(DB[DB$biological_replicate %in% c(1,2),]) %>% print()
dev.off()


