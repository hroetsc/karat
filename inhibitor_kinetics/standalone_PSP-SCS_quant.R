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
source("../../brainstorming/src/invitroSPI_utils.R")
source("src/SCS+PSP-P1.R")

### INPUT ###
finalKinetics = read.csv("aSPIre_manual/results/TSN5inhibitor/finalKinetics.csv", stringsAsFactors = F)
INTtable = read.csv("data/intensity-table-4hrs.csv", stringsAsFactors = F)

# selected peptides
selectedPeps = read.csv("results/b5/selectedPeps_intensities.csv", stringsAsFactors = F)
peptidesToRemove = read.table("results/b5/peptides-to-remove", header = F, sep = " ")

### MAIN PART ###
# ----- preprocessing & formatting -----
DATA = finalKinetics

DATA$spliceType[is.na(DATA$spliceType)] = "PCP"
DATA = na.omit(DATA)
DATA = DATA %>%
  disentangleMultimappers.Type()

# add full replicate info
DATA$condition = paste0(DATA$biological_replicate,
                        "_tech", DATA$technical_replicate)
DATA = DATA[order(DATA$pepSeq), ]


# no synthesis errors
DATA = DATA[gsub("I","L",DATA$pepSeq) %in% gsub("I","L",DATA$pepSeq), ]
# intensity > 1e05
no_idx = c(13:16)
pp = INTtable$pepSeq[rowMeans(INTtable[,3+no_idx]) > 1e05]
DATA = DATA[gsub("I","L",DATA$pepSeq) %in% gsub("I","L",pp), ]

DATA = DATA %>%
  ILredundancy() %>%
  disentangleMultimappers.Type() %>%
  removeMultimappers.Type() %>%
  as.data.frame()

# remove manually identified peptides from list of selected peptides
selectedPeps = selectedPeps[!gsub("I","L",selectedPeps$pepSeq) %in% gsub("I","L",peptidesToRemove$V2), ]
table(selectedPeps$spliceType)

# disentangle time points
DB = DATA %>%
  tidyr::separate_rows(digestTimes, intensities, sep = ";") %>%
  mutate(intensity = as.numeric(intensities),
         digestTime = as.numeric(digestTimes)) %>%
  filter(digestTime == 4 & grepl("noI", biological_replicate))


# ----- plotting -----
plotPSPandSCS = function(DB, name, target) {
  
  out = SCSandPSP_allSubs(DB,target)[[1]]
  yl = max(out$scs_mean+max(out$scs_sd), out$psp_mean+max(out$psp_mean)) %>% ceiling()
  
  p = data.frame(residue = rep(out$residue, 2),
                 n = c(out$scs_n, out$psp_n),
                 p1 = c(out$scs_mean, -1*out$psp_mean),
                 stdev = c(out$scs_sd, out$psp_sd),
                 col = c(rep(plottingCols["PCP"], nrow(out)),
                         rep(plottingCols["PSP"], nrow(out))))
  
  png(paste0("results/b5/",name,"_SCS+PSP.png"),
      height = 6, width = 12, units = "in", res = 300)
  plot(p1~residue,
       data=p,
       type = "h",
       lwd = 4,
       col = col,
       ylim = c(-80, 80),
       xlab = "",
       ylab = "cleavage/splicing strength (%)",
       axes=F)
  arrows(p$residue, p$p1-p$stdev, p$residue, p$p1+p$stdev,
         length=0.03, angle=90, code=3,
         lty = "solid", lwd = 1, col = p$col) %>%
    suppressWarnings()
  
  # positions
  xlabel = c("", paste(DB$substrateSeq[1] %>% strsplit("") %>% unlist(), seq(1, nchar(DB$substrateSeq[1])), sep = ""))
  text(x = seq(0, nchar(DB$substrateSeq[1])), par("usr")[3]-.3, labels = xlabel,
       srt=90, xpd=T, cex=1)
  # number of peptides
  text(x = seq(0, nchar(DB$substrateSeq[1])), par("usr")[3]-15,
       labels = c("n.-spliced:",out$scs_n), srt=0, xpd=T, cex=.8)
  text(x = seq(0, nchar(DB$substrateSeq[1])), par("usr")[3]-22,
       labels = c("spliced:",out$psp_n), srt=0, xpd=T, cex=.8)
  axis(2, at = seq(round(-yl,-1),round(yl,-1),20))
  
  legend("topleft",
         legend = c("cleavage strength (SCS-P1)", "splicing strength (PSP-P1)"),
         lty = c("solid","solid"), col = c(plottingCols["PCP"],plottingCols["PSP"]),
         cex = 1, bty = "n", horiz = T)
  dev.off()
  
}

# ----- execute -----

# intensities at 4 hrs detected in no inhibitor
plotPSPandSCS(DB = DB, name = "TSN5all", target = "P1")
plotPSPandSCS(DB = DB[gsub("I","L",DB$pepSeq) %in% gsub("I","L",selectedPeps$pepSeq), ], name = "TSN5b5specific", target = "P1")


# ----- plot selected peptides -----
cdts = c(4,7,20,24,22)
p1_cdts = c(4,7,8,9,10,11,20,24,22,29)

cdts = data.frame(start = c(1,1,1,1,1,1,16,14,14,25),
                  end = c(4,7,8,9,10,11,20,24,22,29))

plotCandidates = function(DB, name, synthetic=F) {
  
  out = SCSandPSP_allSubs(DB, target = "P1")[[1]]
  yl = max(out$scs_mean+max(out$scs_sd), out$psp_mean+max(out$psp_mean)) %>% ceiling()
  
  p = data.frame(residue = c(out$residue, out$residue+.2),
                 n = c(out$scs_n, out$psp_n),
                 p1 = c(out$scs_mean, out$psp_mean),
                 stdev = c(out$scs_sd, out$psp_sd),
                 col = c(rep(plottingCols["PCP"], nrow(out)),
                         rep(plottingCols["PSP"], nrow(out))))
  
  png(paste0("results/b5/",name,"_SCS+PSP+b5peptides.png"),
      height = 6, width = 12, units = "in", res = 300)
  plot(p1~residue,
       data=p,
       type = "h",
       lwd = 4,
       col = col,
       ylim = if (!synthetic) c(-80,80) else c(-30,80),
       xlab = "",
       ylab = "cleavage/splicing strength (%)",
       axes=F)
  arrows(p$residue, p$p1-p$stdev, p$residue, p$p1+p$stdev,
         length=0.03, angle=90, code=3,
         lty = "solid", lwd = 1, col = p$col) %>%
    suppressWarnings()
  
  # positions
  fonts = rep(1,nchar(DB$substrateSeq[1]))
  fonts[p1_cdts] = 2
  fonts = c(1,fonts)
  
  xlabel = c("", paste(DB$substrateSeq[1] %>% strsplit("") %>% unlist(), seq(1, nchar(DB$substrateSeq[1])), sep = ""))
  text(x = seq(0, nchar(DB$substrateSeq[1])), par("usr")[3]-.3, labels = xlabel,
       srt=90, xpd=T, cex=1, font = fonts)
  # number of peptides
  text(x = seq(0, nchar(DB$substrateSeq[1])), par("usr")[3]-15,
       labels = c("n.-spliced:",out$scs_n), srt=0, xpd=T, cex=.8)
  text(x = seq(0, nchar(DB$substrateSeq[1])), par("usr")[3]-22,
       labels = c("spliced:",out$psp_n), srt=0, xpd=T, cex=.8)
  
  legend("topleft",
         legend = c("cleavage strength (SCS-P1)", "splicing strength (PSP-P1)"),
         lty = c("solid","solid"), col = c(plottingCols["PCP"],plottingCols["PSP"]),
         cex = 1, bty = "n", horiz = T)
  
  # extract
  y = -1
  cnt.intercept = -1
  
  if (!synthetic) {
    DB_b5 = DB[gsub("I","L",DB$pepSeq) %in% gsub("I","L",selectedPeps$pepSeq), ]
    UDB_b5 = DB_b5[-which(duplicated(DB_b5$pepSeq)), ] %>%
      extract_coordinates()
    UDB_b5 = UDB_b5[!str_detect(UDB_b5$positions, ";"), ]
    
    out_b5 = SCSandPSP_allSubs(DB_b5, target = "P1")[[1]]
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
    
  } else {
    for(r in 1:nrow(cdts)) {
      segments(x0 = as.numeric(cdts$start[r])+.1, y0 = cnt.intercept, x1 = as.numeric(cdts$end[r])+.1, y1 = cnt.intercept,
               lwd = 1, col = "black")
      points(x=as.numeric(cdts$end[r])+.1, y=cnt.intercept, pch=16, cex=1)
      cnt.intercept = cnt.intercept-1
    }
  }
  
  
  axis(2, at = seq(0, round(yl, -1), by=20))
  
  dev.off()
}

plotCandidates(DB, name = "TSN5all")
plotCandidates(DB = DB[gsub("I","L",DB$pepSeq) %in% gsub("I","L",selectedPeps$pepSeq), ], name = "TSN5b5specific", synthetic = T)

TSN5all = SCSandPSP_allSubs(DB, target = "P1")[[1]]
TSN5b5specific = SCSandPSP_allSubs(DB[gsub("I","L",DB$pepSeq) %in% gsub("I","L",selectedPeps$pepSeq), ], target = "P1")[[1]]

### OUTPUT ###
save(selectedPeps, file = "results/b5/selectedPeps_final.RData")
write.csv(selectedPeps, file = "results/selectedPeps_final.csv", row.names = F)

save(TSN5all, file = "data/SCS+PSP-P1_TS5all.RData")
save(TSN5b5specific, file = "data/SCS+PSP-P1_TSN5b5specific.RData")

