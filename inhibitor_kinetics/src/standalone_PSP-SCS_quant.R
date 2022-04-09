### INHIBITOR KINETICS ###
# description:  substrate-specific cleavage and splicing strength
#               as defined in Mishto et al. FI 2019
# input:        quantification results (finalKinetics.csv), selected peptides,
#               intensity table to remove synthesis errors
#               list of peptides to remove (selected manually from plotRawIntens.R output)
# output:       PSP-P1 and SCS-P1 at different positions
#               + candidates plotted over substrate length
# author:       HR

library(dplyr)
library(seqinr)
library(stringr)
source("../brainstorming/src/invitroSPI_utils.R")


### INPUT ###
# final Kinetics
DB = read.csv("qiSPI/OUTPUT/TSN5_0+4/finalKinetics.csv", stringsAsFactors = F)
load("qiSPI/OUTPUT/TSN5_0+4/filteredResults.RData")
# filtered out synthesis errors
INTtable = read.csv("data/intensity-table-4hrs.csv", stringsAsFactors = F)

# selected peptides
selectedPeps = read.csv("results/b5/selectedPeps_intensities.csv", stringsAsFactors = F)
peptidesToRemove = read.table("results/b5/peptides-to-remove.txt")

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
                  condition = c("noInh_bio1_tech1", "noInh_bio1_tech2", "noInh_bio2_tech1", "noInh_bio2_tech2",
                                "b5_bio1_tech1", "b5_bio1_tech2", "b5_bio2_tech1", "b5_bio2_tech2",
                                "b2_bio1_tech1", "b2_bio1_tech2", "b2_bio2_tech1", "b2_bio2_tech2",
                                "b1_bio1_tech1", "b1_bio1_tech2", "b1_bio2_tech1", "b1_bio2_tech2"),
                  biological_replicate = rep(seq(1,8), each=2))
res = left_join(res, cond) %>% as.data.frame()
names(res)[2] = "pepSeq"

DB = left_join(DB, res) %>%
  as.data.frame()

# no synthesis errors
DB = DB[DB$pepSeq %in% INTtable$pepSeq, ]

# intensity > 1e05
no_idx = c(1:4)
pp = INTtable$pepSeq[rowMeans(INTtable[,4+no_idx]) > 1e05]
DB = DB[DB$pepSeq %in% pp, ]

DB = DB %>%
  ILredundancy() %>%
  disentangleMultimappers.Type() %>%
  removeMultimappers.Type() %>%
  as.data.frame()

# remove manually identified peptides from list of selected peptides
selectedPeps = selectedPeps[!gsub("I","L",selectedPeps$pepSeq) %in% gsub("I","L",peptidesToRemove$V1), ]

table(selectedPeps$spliceType)

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

# ---------- relevant info ----------
# positions
SR1pos = c("P4"=-3, "P3"=-2, "P2"=-1, "P1"=0,
           "P-1"=1, "P-2"=2, "P-3"=3, "P-4"=4)

SR2pos = c("P-4_"=-4, "P-3_"=-3, "P-2_"=-2, "P-1_"=-1,
           "P1_"=0, "P2_"=1, "P3_"=2, "P4_"=3)

PCPpos = c("P4"=-3, "P3"=-2, "P2"=-1, "P1"=0,
           "P1_"=1, "P2_"=2, "P3_"=3, "P4_"=4)

SRnames = c(names(SR1pos), names(SR2pos))
types = c("cis", "revCis", "trans", "PCP")

# ---------- extract amino acids ----------
# extract coordinates of different positions in the substrate
extract_coordinates = function(tbl = ""){
  
  tbl$spliceType[(tbl$spliceType == "") | (is.na(tbl$spliceType))] = "PCP"
  
  # table with position indices
  pos = str_split_fixed(tbl$positions, coll("_"), Inf) %>% as.data.frame()
  pos = apply(pos, 2, function(x){as.numeric(as.character(x))})
  
  pcp = which(tbl$spliceType == "PCP")
  psp = which(tbl$spliceType != "PCP")
  
  
  # PCPs
  pcpTBL = sapply(PCPpos, function(x){
    # substr(tbl$substrateSeq[pcp], start = pos[pcp,2]+x, stop = pos[pcp,2]+x)
    pos[pcp,2]+x
  }) %>%
    as.data.frame() %>%
    mutate(spliceType = tbl$spliceType[pcp],
           positions = tbl$positions[pcp],
           pepSeq = tbl$pepSeq[pcp],
           substrateID = tbl$substrateID[pcp],
           substrateSeq = tbl$substrateSeq[pcp])
  
  
  # PSPs
  pspSR1TBL = sapply(SR1pos, function(x){
    # substr(tbl$substrateSeq[psp], start = pos[psp,2]+x, stop = pos[psp,2]+x)
    pos[psp,2]+x
  }) %>%
    as.data.frame() %>%
    mutate(spliceType = tbl$spliceType[psp],
           positions = tbl$positions[psp],
           pepSeq = tbl$pepSeq[psp],
           substrateID = tbl$substrateID[psp],
           substrateSeq = tbl$substrateSeq[psp])
  
  
  pspSR2TBL = sapply(SR2pos, function(x){
    # substr(tbl$substrateSeq[psp], start = pos[psp,3]+x, stop = pos[psp,3]+x)
    pos[psp,3]+x
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

DBaa = extract_coordinates(DB)

DBMaster = inner_join(DB, DBaa, by=c("pepSeq", "positions")) %>%
  unique() %>%
  rename(spliceType = spliceType.x,
         substrateID = substrateID.x,
         substrateSeq = substrateSeq.x)


# ----- compute SCS-P1 and PSP-P1 -----
SCS_and_PSP = function(DBMaster,target) {
  
  DBMaster = resolve_multimapper(DBMaster)
  DBMaster$productType = toupper(DBMaster$productType)
  
  pcpidx = which(DBMaster$productType == "PCP")
  pspidx = which(DBMaster$productType == "PSP")
  
  res = str_split_fixed(DBMaster$positions, coll("_"), n = Inf) %>%
    as.data.frame()
  
  d = data.frame(residue = c(1:nchar(DBMaster$substrateSeq[1])),
                 scs_mean = 0,
                 scs_sd = 0,
                 scs_n = 0,
                 psp_mean = 0,
                 psp_sd = 0,
                 psp_n = 0)
  
  for (i in 1:nrow(d)) {
    
    kpcp = which(DBMaster[,target] == d$residue[i] & DBMaster$productType == "PCP")
    kpsp = which(DBMaster[,target] == d$residue[i] & DBMaster$productType == "PSP")
    # add the C-terminus of spliced peptides as cleaved residue
    if (target == "P1") {
      kpcp = c(kpcp, which(res$V4[pspidx] == d$residue[i]))
    }
    
    cntscs = DBMaster[kpcp,] %>%
      dplyr::group_by(biological_replicate) %>%
      dplyr::summarise(int = sum(intensity))
    
    d$scs_mean[i] = paste(cntscs$int, collapse = "_")
    d$scs_n[i] = DBMaster$pepSeq[kpcp] %>%
      unique() %>%
      length()
    
    cntpsp = DBMaster[kpsp,] %>%
      dplyr::group_by(biological_replicate) %>%
      dplyr::summarise(int = sum(intensity))
    
    d$psp_mean[i] = paste(cntpsp$int, collapse = "_")
    d$psp_n[i] = DBMaster$pepSeq[kpsp] %>%
      unique() %>%
      length()
  }
  
  d[is.na(d)] = 0
  d[d == ""] = 0
  
  # remove substrate's C terminus from the SCS
  d$scs_mean[nchar(DBMaster$substrateSeq[1])] = 0
  d$scs_sd[nchar(DBMaster$substrateSeq[1])] = 0
  
  # normalise by sum for each biological replicate
  scs = apply(str_split(d$scs_mean, pattern = "_", simplify = T), 2, as.numeric) %>%
    as.data.frame()
  scs[is.na(scs)] = 0
  
  psp = apply(str_split(d$psp_mean, pattern = "_", simplify = T), 2, as.numeric) %>%
    as.data.frame()
  psp[is.na(psp)] = 0
  
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
target = "P1"
DBMaster = DB_b5

plotPSPandSCS = function(DBMaster, name, target) {
  
  out = SCS_and_PSP(DBMaster,target)
  yl = max(out$scs_mean+max(out$scs_sd), out$psp_mean+max(out$psp_mean)) %>% ceiling()
  
  p = data.frame(residue = rep(out$residue, 2),
                 n = c(out$scs_n, out$psp_n),
                 p1 = c(out$scs_mean, -1*out$psp_mean),
                 stdev = c(out$scs_sd, out$psp_sd),
                 col = c(rep(plottingCols["PCP"], nrow(out)),
                         rep(plottingCols["PSP"], nrow(out))))
  
  # pdf(paste0("results/b5/PSP-SCS_", name, ".pdf"), height = 6, width = 12)
  png("~/Documents/Studium/Fachvertiefung+BA/Praktikumsbericht/plots/SCS+PSP_b5.png",
      height = 4, width = 10, units = "in", res = 300)
  plot(p1~residue,
       data=p,
       type = "h",
       lwd = 4,
       col = col,
       # main = name,
       ylim = c(-yl, yl),
       # xlab = "substrate residue forming SCS-/PSP-P1",
       xlab = "",
       # ylab = "substrate-specific cleavage/splicing strength (%)",
       ylab = "cleavage/splicing strength (%)",
       # sub = "mean +- S.D.",
       axes=F)
  arrows(p$residue, p$p1-p$stdev, p$residue, p$p1+p$stdev,
         length=0.03, angle=90, code=3,
         lty = "solid", lwd = 1, col = p$col) %>%
    suppressWarnings()
  xlabel = c("",
             paste(DBMaster$substrateSeq[1] %>% strsplit("") %>% unlist(), seq(1, nchar(DBMaster$substrateSeq[1])), sep = ""))
  # positions
  text(x = seq(0, nchar(DBMaster$substrateSeq[1])), par("usr")[3]-.3, labels = xlabel,
       srt=90, xpd=T, cex=.5)
  # number of peptides
  text(x = seq(0, nchar(DBMaster$substrateSeq[1])), par("usr")[3]-15,
       labels = c("non-spliced:",out$scs_n), srt=0, xpd=T, cex=.5)
  text(x = seq(0, nchar(DBMaster$substrateSeq[1])), par("usr")[3]-22,
       labels = c("spliced:",out$psp_n), srt=0, xpd=T, cex=.5)
  axis(2, at = seq(round(-yl,-1),round(yl,-1),20))
  
  legend("topleft",
         legend = c("cleavage strength (SCS-P1)", "splicing strength (PSP-P1)"),
         lty = c("solid","solid"), col = c(plottingCols["PCP"],plottingCols["PSP"]),
         cex = .7, bty = "n")
  dev.off()
  
  
  
  pdf(paste0("results/b5/scatterPSP-SCS", name, ".pdf"),
      height = 6, width = 6)
  plot(psp_mean~scs_mean, data = out,
       main = name,
       pch = 16,
       xlim = c(0, yl),
       ylim = c(0, yl),
       cex = .8,
       ylab = "mean PSP-P1 (%)",xlab = "mean SCS-P1 (%)")
  dev.off()
  
}

# ----- execute -----

# intensities at 4 hrs detected in no inhibitor
# only bioReps that are not b5
DBMaster_noInh = DBMaster[DBMaster$biological_replicate %in% c(1,2),]
DBMaster_noInh$intensity = as.numeric(DBMaster_noInh$`4`)
plotPSPandSCS(DBMaster = DBMaster_noInh, name = "TSN5all", target = "P1")

# plot the same for the b5 intensities
DBMaster_b5inh = DBMaster[DBMaster$biological_replicate %in% c(3,4),]
DBMaster_b5inh$intensity = as.numeric(DBMaster_b5inh$`4`)
plotPSPandSCS(DBMaster = DBMaster_b5inh, name = "TSN5all_b5intens", target = "P1")

# only selected peptides
DB_b5 = DBMaster_noInh[DBMaster_noInh$pepSeq %in% gsub("I","L",selectedPeps$pepSeq), ]
plotPSPandSCS(DBMaster = DB_b5, name = "TSN5_b5_P1", target = "P1")
plotPSPandSCS(DBMaster = DB_b5, name = "TSN5_b5_P2", target = "P2")
plotPSPandSCS(DBMaster = DB_b5, name = "TSN5_b5_P3", target = "P3")




# ----- plot selected peptides -----
ProteasomeDB = DBMaster_noInh

plotCandidates = function(ProteasomeDB, DB_b5) {
  
  out = SCS_and_PSP(ProteasomeDB, target = "P1")
  yl = max(out$scs_mean+max(out$scs_sd), out$psp_mean+max(out$psp_mean)) %>% ceiling()
  
  p = data.frame(residue = c(out$residue, out$residue+.2),
                 n = c(out$scs_n, out$psp_n),
                 p1 = c(out$scs_mean, out$psp_mean),
                 stdev = c(out$scs_sd, out$psp_sd),
                 col = c(rep(plottingCols["PCP"], nrow(out)),
                         rep(plottingCols["PSP"], nrow(out))))
  
  png("~/Documents/Studium/Fachvertiefung+BA/Praktikumsbericht/plots/SCS+PSP_all.png",
      height = 6, width = 10, units = "in", res = 300)
  plot(p1~residue,
       data=p,
       type = "h",
       lwd = 4,
       col = col,
       ylim = c(-50,yl),
       xlab = "",
       ylab = "cleavage/splicing strength (%)",
       axes=F)
  arrows(p$residue, p$p1-p$stdev, p$residue, p$p1+p$stdev,
         length=0.03, angle=90, code=3,
         lty = "solid", lwd = 1, col = p$col) %>%
    suppressWarnings()
  xlabel = c("",
             paste(ProteasomeDB$substrateSeq[1] %>% strsplit("") %>% unlist(), seq(1, nchar(ProteasomeDB$substrateSeq[1])), sep = ""))
  # positions
  text(x = seq(0, nchar(ProteasomeDB$substrateSeq[1])), par("usr")[3]-5, labels = xlabel,
       srt=90, xpd=T, cex=.5)
  # number of peptides
  text(x = seq(0, nchar(ProteasomeDB$substrateSeq[1])), par("usr")[3]-15,
       labels = c("non-spliced:",out$scs_n), srt=0, xpd=T, cex=.5)
  text(x = seq(0, nchar(ProteasomeDB$substrateSeq[1])), par("usr")[3]-20,
       labels = c("spliced:",out$psp_n), srt=0, xpd=T, cex=.5)
  
  # extract
  y = -1
  cnt.intercept = -1
  
  UDB_b5 = DB_b5[-which(duplicated(DB_b5$pepSeq)), ]
  UDB_b5 = UDB_b5[!str_detect(UDB_b5$positions, ";"), ]
  
  out_b5 = SCS_and_PSP(DB_b5, target = "P1")
  resi = out_b5$residue[out_b5$psp_mean > 0 | out_b5$scs_mean > 0]
  for (r in resi) {
    kk = which(UDB_b5$P1 == r)
    
    for (i in 1:length(kk)) {
      pos = str_split(UDB_b5$positions[kk[i]], coll("_"), simplify = T)
      
      segments(x0 = as.numeric(pos[1])+.1, y0 = cnt.intercept, x1 = as.numeric(pos[2])+.1, y1 = cnt.intercept,
               lwd = .5, col = if(length(pos) == 4) {plottingCols[["PSP"]]} else {plottingCols[["PCP"]]})
      points(x=r+.1, y=cnt.intercept, pch=16, cex=.3)
      cnt.intercept = cnt.intercept-.5
      
    }
    
  }
  
  axis(2, at = seq(0, round(yl, -1), by=20))
  
  legend("topleft",
         legend = c("cleavage strength (SCS-P1)", "splicing strength (PSP-P1)"),
         lty = c("solid","solid"), col = c(plottingCols["PCP"],plottingCols["PSP"]),
         cex = .7, bty = "n")
  dev.off()
}

pdf(paste0("results/b5/candidates.pdf"), height = 6, width = 12)
plotCandidates(DBMaster_noInh, DB_b5) %>% print()
dev.off()

pdf(paste0("results/b5/candidates_b5intens.pdf"), height = 6, width = 12)
plotCandidates(DBMaster_b5inh, DB_b5) %>% print()
dev.off()

### OUTPUT ###
save(selectedPeps, file = "results/b5/selectedPeps_final.RData")




nchar(DB_b5$pepSeq) %>% density() %>% plot()
nchar(DBMaster_noInh$pepSeq) %>% density() %>% lines()


