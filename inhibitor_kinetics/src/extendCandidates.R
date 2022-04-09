### INHIBITOR KINETICS ###
# description:  extend selected candidates C-terminally and check if those peptides are being detected
# input:        candidate peptides, DB (no synthesis errors, int > 1e05)
# output:       -
# author:       HR


library(dplyr)
library(stringr)
source("../brainstorming/src/invitroSPI_utils.R")


### INPUT ###
DB = read.csv("qiSPI/OUTPUT/TSN5_0+4/finalKinetics.csv", stringsAsFactors = F)
INTtable = read.csv("data/intensity-table-4hrs.csv", stringsAsFactors = F)
selectedPeps = read.csv("results/b5/selectedPeps_intensities.csv", stringsAsFactors = F)

### MAIN PART ###
# ----- preprocessing -----
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
  disentangleMultimappers.AA() %>%
  removeMultimappers.AA()
  as.data.frame()


# ----- candidates -----
cdts = data.frame(S4 = c(1,4),  # PSP and PCP
                  F7 = c(1,7),  # only PSP
                  F7_2 = c(1,8),
                  F7_3 = c(1,9),
                  F7_4 = c(1,10),
                  F7_5 = c(1,11),
                  S20 = c(16,20),
                  A24 = c(14,24),
                  N29 = c(25,29),  # nothing
                  L22 = c(14,22)) %>% t() %>% as.data.frame()  # only PCP
names(cdts) = c("start","end")

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

# ----- plot candidates -----
# from plot_SCS_and_PSPS.R
out = SCS_and_PSP(DB_b5, target = "P1")
yl = max(out$scs_mean+max(out$scs_sd), out$psp_mean+max(out$psp_mean)) %>% ceiling()

p = data.frame(residue = c(out$residue, out$residue+.2),
               n = c(out$scs_n, out$psp_n),
               p1 = c(out$scs_mean, out$psp_mean),
               stdev = c(out$scs_sd, out$psp_sd),
               col = c(rep(plottingCols["PCP"], nrow(out)),
                       rep(plottingCols["PSP"], nrow(out))))

png("~/Documents/Studium/Fachvertiefung+BA/project_proposal/plots/candidates.png",
    height = 6, width = 10, units = "in", res = 300)
plot(p1~residue,
     data=p,
     type = "h",
     lwd = 3,
     col = col,
     ylim = c(-10,yl),
     xlab = "",
     ylab = "cleavage/splicing strength (%)",
     axes=F)
suppressWarnings(arrows(p$residue, p$p1-p$stdev, p$residue, p$p1+p$stdev,
       length=0.03, angle=90, code=3,
       lty = "solid", lwd = 1, col = p$col))
xlabel = c("",
           paste(ProteasomeDB$substrateSeq[1] %>% strsplit("") %>% unlist(), seq(1, nchar(ProteasomeDB$substrateSeq[1])), sep = ""))
# positions
text(x = seq(0, nchar(ProteasomeDB$substrateSeq[1])), par("usr")[3]-5, labels = xlabel,
     srt=90, xpd=T, cex=.5)
# number of peptides
text(x = seq(0, nchar(ProteasomeDB$substrateSeq[1])), par("usr")[3]-8,
     labels = c("non-spliced:",out$scs_n), srt=0, xpd=T, cex=.5)
text(x = seq(0, nchar(ProteasomeDB$substrateSeq[1])), par("usr")[3]-10,
     labels = c("spliced:",out$psp_n), srt=0, xpd=T, cex=.5)

# extract
y = -1
cnt.intercept = -1

for (i in 1:nrow(cdts)) {
  segments(x0 = cdts[i,1]+.1, y0 = cnt.intercept, x1 = cdts[i,2]+.1, y1 = cnt.intercept,
           lwd = .5)
  points(x=as.numeric(str_extract(rownames(cdts)[i], "[:digit:]+"))+.1, y=cnt.intercept,
         pch=16, cex=.5)
  cnt.intercept = cnt.intercept-.8
}

legend("topright",
       legend = c("cleavage strength (SCS-P1)", "splicing strength (PSP-P1)"),
       lty = c("solid","solid"), col = c(plottingCols["PCP"],plottingCols["PSP"]),
       cex = .9, bty = "n")
axis(2, at = seq(0, round(yl, -1), by=20))
dev.off()

# ----- check for precursors -----
# when extending the candidates towards the C-term
# how many peptides do we find that have the C-terminal aa as P1
DBU = DBMaster[-which(duplicated(DBMaster$pepSeq)), ]
pos = str_split_fixed(DBU$positions, coll("_"), Inf)[,c(1:4)]
pos = apply(pos,2,as.numeric)

cdts$CtermP1 = NA
comp = list()
for (i in 1:nrow(cdts)) {
  
  hit = F
  count = 1
  while(!hit) {
    
    k = which(DBU$P1 == cdts$end[i]+count & pos[,1]==cdts$start[i])
    if (length(k) > 0) {
      hit = T
    }
    
    count = count+1
  }
  
  comp[[i]] = DBU[k,]
  names(comp)[i] = rownames(cdts)[i]
  cdts$CtermP1[i] = cdts$end[i]+count-1
}

# for all candidates we find the C-terminal 1aa precursor
cdts

### OUTPUT ###
save(cdts, file = "results/b5/selected-candidates.RData")
save(comp, file = "results/b5/selected-candidates-CtermPrec.RData")

