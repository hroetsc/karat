### INHIBITOR KINETICS ###
# description:  plot coverage map of selected peptides
# input:        b5 specific peptides
# output:       coverage map
# author:       HR


library(dplyr)
library(stringr)

source("../brainstorming/src/invitroSPI_utils.R")

### INPUT ###
pepList = read.csv("results/b5/selectedPeps_intensities.csv", stringsAsFactors = F)
table(pepList$spliceType)

### MAIN PART ###
# ----- get mean intensities (no inhibitor) -----
names(pepList)
meanInt = log(rowMeans(pepList[, c(5:8)]))

# ----- plot coverage -----
plotCoverage = function(ProteasomeDB) {
  
  k = which(str_detect(ProteasomeDB$positions, ";"))
  if (length(k) > 0) {
    print("removing all multi-mappers")
    ProteasomeDB = ProteasomeDB[-k,]
  }
  
  pcp = which(ProteasomeDB$spliceType == "PCP")
  psp = which(ProteasomeDB$spliceType != "PCP")
  
  pcp_pos = ProteasomeDB$positions[pcp]
  psp_pos = ProteasomeDB$positions[psp]
  S = ProteasomeDB$substrateSeq[1]
  
  # --- PCPs ---
  pcp_pos = str_split(pcp_pos, "_", simplify = T)
  pcp_pos = apply(pcp_pos,2,as.numeric) %>%
    as.data.frame()
  pcp_pos = pcp_pos[order(pcp_pos$V1, decreasing = F), ]
  
  pdf(paste0("results/b5/coverage_map.pdf"),
      height = 12, width = 16)
  
  par(mfrow = c(2,1))
  y0 =0.1
  plot(NULL, ylim=c(0,nrow(pcp_pos)*0.05),
       ylab = "", xlab="substrate",
       xlim = c(0, nchar(S)),
       axes=F)
  
  for (i in 1:nrow(pcp_pos)) {
    
    if (length(which(pcp_pos$V1[1:i] <= pcp_pos$V1[i] & 
                     pcp_pos$V2[1:i] >= pcp_pos$V1[i])) > 1) {
      
      cnt.intercept = length(which(pcp_pos$V1[1:i] <= pcp_pos$V1[i] & 
                                     pcp_pos$V2[1:i] >= pcp_pos$V1[i]))*0.1
      
    } else {cnt.intercept = y0}
    
    segments(x0 = pcp_pos$V1[i], y0 = cnt.intercept,
             x1 = pcp_pos$V2[i], y1 = cnt.intercept,
             col = plottingCols["PCP"],
             lwd = meanInt[pcp][i]/max(meanInt[pcp]))
  }
  
  xlabel = c("0",
             paste(S %>% strsplit("") %>% unlist(), seq(1, nchar(S)), sep = ""))
  text(x = seq(0, nchar(S)), par("usr")[3]-.3, labels = xlabel,
       srt=90, xpd=T, cex=.5)
  axis(2, labels = NA)
  
  
  # --- PSPs ---
  psp_pos = str_split(psp_pos, "_", simplify = T)
  psp_pos = apply(psp_pos,2,as.numeric) %>%
    as.data.frame()
  psppos = psp_pos[order(psp_pos$V1, decreasing = F), ]
  
  psp_pos = rbind(as.matrix(psppos[,c(1:2)]), as.matrix(psppos[,c(3:4)])) %>%
    as.data.frame()
  
  y0 =0.1
  plot(NULL, ylim=c(0,nrow(psp_pos)*0.05),
       ylab = "", xlab="substrate",
       xlim = c(0, nchar(S)),
       axes=F)
  
  for (i in 1:nrow(psp_pos)) {
    
    if (length(which(psp_pos$V1[1:i] <= psp_pos$V1[i] & 
                     psp_pos$V2[1:i] >= psp_pos$V1[i])) > 1) {
      
      cnt.intercept = length(which(psp_pos$V1[1:i] <= psp_pos$V1[i] & 
                                     psp_pos$V2[1:i] >= psp_pos$V1[i]))*0.1
      
    } else {cnt.intercept = y0}
    
    segments(x0 = psp_pos$V1[i], y0 = cnt.intercept,
             x1 = psp_pos$V2[i], y1 = cnt.intercept,
             col = plottingCols["PSP"],
             lwd = meanInt[psp][i]/max(meanInt[psp]))
  }
  
  xlabel = c("0",
             paste(S %>% strsplit("") %>% unlist(), seq(1, nchar(S)), sep = ""))
  text(x = seq(0, nchar(S)), par("usr")[3]-.3, labels = xlabel,
       srt=90, xpd=T, cex=.5)
  axis(2, labels = NA)
  
  dev.off()
}

plotCoverage(ProteasomeDB = pepList)
