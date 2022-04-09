### INHIBITOR KINETICS ###
# description:  plot raw intensities of selected peptides for 0,2,4 hrs
#               peptides downregulated in b1 and b2 inhibition (irrgardless of significance)
# input:        selectedPeps, qiSPI output for 0+4 and 0+2 hrs
# output:       plots for selected peptides
# author:       HR


library(dplyr)
library(stringr)
source("../brainstorming/src/invitroSPI_utils.R")

# subunit = "b2"

### INPUT ###
# selected peptides
# fs_peps = list.files(paste0("results/",subunit), pattern = "_statistics.csv", full.names = T, recursive = T)

# 0+4
DB4 = read.csv("qiSPI/OUTPUT/TSN5_0+4/finalKinetics.csv", stringsAsFactors = F)
load("qiSPI/OUTPUT/TSN5_0+4/filteredResults.RData")
filteredResults4 = filteredResults

# 0+2
DB2 = read.csv("qiSPI/OUTPUT/TSN5_0+2/finalKinetics.csv", stringsAsFactors = F)
load("qiSPI/OUTPUT/TSN5_0+2/filteredResults.RData")
filteredResults2 = filteredResults

### MAIN PART ###
# ----- get peptide sequences -----
# # read in all sequences with p-values and fold changes
# selectedPeps = lapply(fs_peps, read.csv, stringsAsFactors=F) %>%
#   plyr::ldply() %>%
#   as.data.frame()
# 
# # select those that have a negatice fold change irregardless of magnitude/p-val
# selectedPeps = selectedPeps[selectedPeps$estimate<0, ]
# selectedPeps = selectedPeps[-which(duplicated(selectedPeps$X)), ]

# assign positions and splice types
DB = DB4 %>%
  ILredundancy() %>%
  disentangleMultimappers.Type()
DB$spliceType[is.na(DB$spliceType)] = "PCP"

# ----- preprocessing -----
getIntensities = function(filteredResults) {
  
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
  
  return(res)
}

int2hrs = getIntensities(filteredResults2)
int4hrs = getIntensities(filteredResults4)

# # only selected peptides
# int2hrs = int2hrs[int2hrs$pepSeq %in% selectedPeps$gene, ]
# int4hrs = int4hrs[int4hrs$pepSeq %in% selectedPeps$gene, ]

# ----- plotting -----
pepSeqs = selectedPeps$gene
pepSeqs = unique(DB4$pepSeq)

grps = unique(str_extract_all(int4hrs$condition, "^[:alnum:]+_bio[:digit:]", simplify = T)) %>%
  as.character() %>%
  sort()
cols = c("black", "lightblue", "purple","gray")
matching = data.frame(group = grps, cols = rep(cols,each=2))


# pdf(paste0("results/",subunit,"/rawIntensities-for-downregulated-peps.pdf"), height = 16, width = 16)
pdf(paste0("results/rawIntensities-all.pdf"), height = 16, width = 16)
par(mfrow = c(4,4))
for (i in 1:length(pepSeqs)) {
  
  # get intensities
  cnt4 = int4hrs[int4hrs$pepSeq == pepSeqs[i],] %>%
    mutate(group = str_extract_all(condition, "^[:alnum:]+_bio[:digit:]", simplify = T)) %>%
    left_join(matching)
  
  cnt2 = int2hrs[int2hrs$pepSeq == pepSeqs[i],]
  if (nrow(cnt2) > 0) {
    cnt4$`2` = cnt2$`2`[match(cnt4$condition,cnt2$condition)]
  } else {
    cnt4$`2` = NA
  }
  
  # mean over technical replicates
  names(cnt4)
  tmp = cnt4 %>%
    tidyr::gather(time, intensity, -sample, -pepSeq, -types, -condition, -biological_replicate, -group, -cols) %>% 
    group_by(time,group,cols) %>%
    summarise(mean_int = mean(intensity),
              sd_int = sd(intensity))
  k = which(!is.na(tmp$mean_int))
  
  # I guess we need to normalise the data too...
  
  # plot it
  plot(x = tmp$time, y = tmp$mean_int,
       pch= 16, cex=.8,
       col = tmp$cols,
       xlab = "time [hrs]", ylab = "intensity",
       main = paste0(pepSeqs[i]," - ",DB4$spliceType[DB4$pepSeq == gsub("I","L",pepSeqs[i])][1]),
       sub = paste0("positions: ",DB4$positions[DB4$pepSeq == gsub("I","L",pepSeqs[i])][1]))
  for (g in unique(tmp$group)) {
    lines(x = as.numeric(tmp$time[k][tmp$group == g]), y = tmp$mean_int[k][tmp$group == g], col = tmp$cols[tmp$group==g][1])
    arrows(as.numeric(tmp$time[tmp$group == g]), tmp$mean_int[tmp$group == g]-tmp$sd_int[tmp$group == g],
           as.numeric(tmp$time[tmp$group == g]), tmp$mean_int[tmp$group == g]+tmp$sd_int[tmp$group == g],
           length=0.03, angle=90, code=3,
           lty = "solid", lwd = 1, col = tmp$cols[tmp$group==g][1])
  }
  
  legend("topleft", cex=.8,
         legend = c("no", "ß1", "ß2", "ß5"),
         lty = rep("solid",4), col = c("gray","black", "lightblue", "purple"))
  
}

dev.off()


