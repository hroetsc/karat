### INHIBITOR KINETICS ###
# description:  parse qiSPI output for DE analysis
# input:        qiSPI filtered resulst (from 9_filterKinetics.R), unfiltered results,
#               substrate degradation determined using XICs
# output:       table with intensities for all peptides
# author:       HR

library(dplyr)
library(stringr)
source("qiSPI/SOURCE/qiSPI_utils.R")
source("../brainstorming/src/invitroSPI_utils.R")
source("src/data_utils.R")

bioReps = 8
techReps = 2

### INPUT ###
load("qiSPI/OUTPUT/TSN5_0+4/filteredResults.RData")
load("qiSPI/OUTPUT/TSN5_0+4/quantity_rep.RData")
sample_list = read.csv("qiSPI/INPUT/sample_list.csv", stringsAsFactors = F)
load("data/substrate-degradation.RData")

### MAIN PART ###
# ----- preprocessing & formatting -----
# map names of filteredResults
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

# get positions and splice type
# substrateSeq, productType, spliceType, pepSeq
DATA = res %>%
  mutate(productType = toupper(types),
         spliceType = NA,
         substrateSeq = sample_list$substrateSeq[1],
         pepSeq = sequence) %>%
  mapping() %>%
  disentangleMultimappers.Type()

DATA = DATA[order(DATA$pepSeq), ]

# some stats
DATA$spliceType[is.na(DATA$spliceType)] = "PCP"
table(DATA$spliceType[DATA$condition == "noInh_bio1_tech1"])

# ----- fetch intensities -----
NuniquePeps = DATA$pepSeq %>% unique() %>% length()

getIntensityTable = function(target) {
  
  INTtable = matrix(NA, nrow = nrow(DATA), ncol = 4+bioReps*techReps) %>%
    as.data.frame()
  names(INTtable) = c("substrateSeq", "spliceType", "positions", "pepSeq",
                      cond$condition %>% unique())
  intIDX = which(names(INTtable) %in% unique(cond$condition))
  
  DATA = DATA[order(DATA$pepSeq, DATA$biological_replicate), ]
  k = 1
  counter = 1
  while (k <= NuniquePeps) {
    
    cnt = DATA[c(counter):c(counter+bioReps*techReps-1), ]
    intensity = cnt[,c(target)] %>% as.numeric()
    
    f = all(cnt$biological_replicate == rep(seq(1,bioReps), each=techReps))
    if(!f) {print("SOMETHING IS WRONG")}
    
    # sequence info
    INTtable$substrateSeq[k] = cnt$substrateSeq[1]
    INTtable$spliceType[k] = cnt$spliceType[1]
    INTtable$positions[k] = cnt$positions[1]
    INTtable$pepSeq[k] = cnt$pepSeq[1]
    
    # intensities
    INTtable[k, intIDX] = intensity[match(cnt$condition, names(INTtable)[intIDX])]
    
    k = k+1
    counter = counter+techReps*bioReps
  }
  
  INTtable = INTtable[1:c(k-1),]
  
  # sanity check
  NuniquePeps == INTtable$pepSeq %>% unique() %>% length()
  apply(INTtable, 2, function(x){any(is.na(x))})
  
  INTtable$spliceType[is.na(INTtable$spliceType)] = "PCP"
  
  # remove substrate
  # k2 = which(INTtable$pepSeq == INTtable$substrateSeq)
  # if (length(k2) > 0) {
  #   INTtable = INTtable[-k2,]
  # }
  
  return(INTtable)
}

INT_0hrs = getIntensityTable(target = "0")
INT_4hrs = getIntensityTable(target = "4")

intIDX = which(names(INT_0hrs) %in% unique(cond$condition))
I0 = as.matrix(INT_0hrs[,intIDX])
I4 = as.matrix(INT_4hrs[,intIDX])

rownames(I0) = INT_0hrs$pepSeq
rownames(I4) = INT_4hrs$pepSeq


# ----- remove synthesis errors -----
# for no inhibitor: remove everything that has a lower intensity at 4 hrs than at 0 hrs
# mean over technical and biological replicates

no_idx = c(1:4)
b1_idx = c(13:16)
b2_idx = c(9:12)
b5_idx = c(5:8)

rem = which(rowMeans(I4[,no_idx]) < rowMeans(I0[,no_idx]))
length(rem) / nrow(I4)

I0 = I0[-rem,]
I4 = I4[-rem,]
INT_0hrs = INT_0hrs[-rem,]
INT_4hrs = INT_4hrs[-rem,]

table(INT_4hrs$spliceType)

# ---- statistics -----

apply(INT_0hrs[,intIDX], 2, summary)
apply(INT_4hrs[,intIDX], 2, summary)

diagnostics = function(M) {
  
  density(na.omit(M[,1])) %>% plot()
  apply(M, 2, function(x){
    lines(density(x[is.finite(x)]), col="red")
  })
  
  apply(M,2,function(x){
    x[!is.finite(x)] = 0
    y = sort(x)
    return(y)
  }) %>%
    pairs(pch = 16, lower.panel = NULL, cex = .2, alpha=.8) %>%
    print()
  
}


# ----- get substrate degradation -----
# from Distiller XIC
# how?!
# sk = which(INT_0hrs$substrateSeq == INT_0hrs$pepSeq)
# subsDeg = I4[sk,] - I0[sk,]
# subsDeg
# subsDeg/max(I4) * 100


degraded$degraded_norm = degraded$degraded - min(degraded$degraded) +1

colnames(I4)
cnames = str_extract_all(colnames(I4), "^[:alnum:]+_bio[:digit:]", simplify = T)
dg = degraded$degraded_norm[match(cnames, degraded$condition)]
DEGR = matrix(rep(dg, each=nrow(I4)), nrow = nrow(I4))

# ----- data normalisation -----
## normalise per mol degraded substrate
I4_perSubs = I4/DEGR
I0_perSubs = I0/DEGR

FC_perSubs = ((I4_perSubs-min(I4_perSubs)+1)/(I0_perSubs-min(I0_perSubs)+1)) - 1
FC_perSubs = FC_perSubs - min(FC_perSubs, na.rm = T)
FC_perSubs = log(FC_perSubs+1)


## quantile normalisation over all time points/ technical replicates
## normalise 4 hrs intensities by the substrate degradation

## set 0 intensities to NA
I0[I0==0] = NA
I4[I4==0] = NA

## plot 5% quantiles
quantiles0 = apply(I0, 2, function(x){
  quantile(x, 0.05, na.rm=T)
})
quantiles4 = apply(I4, 2, function(x){
  quantile(x, 0.05, na.rm=T)
})

plot(log10(quantiles0), pch=16, col="black", ylim = c(0,5),
     main = "log10 5% quantiles of intensities",
     sub = "black: 0 hrs, green: 4 hrs")
points(log10(quantiles4), pch=16, col="limegreen")
abline(v = c(4.5, 8.5, 12.5), lty="dashed")


## background correction using the 5% quantile
q_no = quantile(rbind(I0[,no_idx], I4[,no_idx]), 0.05, na.rm = T)
q_b1 = quantile(rbind(I0[,b1_idx], I4[,b1_idx]), 0.05, na.rm = T)
q_b2 = quantile(rbind(I0[,b2_idx], I4[,b2_idx]), 0.05, na.rm = T)
q_b5 = quantile(rbind(I0[,b5_idx], I4[,b5_idx]), 0.05, na.rm = T)

I0[,no_idx] = I0[,no_idx]/q_no
I4[,no_idx] = I4[,no_idx]/q_no

I0[,b1_idx] = I0[,b1_idx]/q_b1
I4[,b1_idx] = I4[,b1_idx]/q_b1

I0[,b2_idx] = I0[,b2_idx]/q_b2
I4[,b2_idx] = I4[,b2_idx]/q_b2

I0[,b5_idx] = I0[,b5_idx]/q_b5
I4[,b5_idx] = I4[,b5_idx]/q_b5

# plot(dg~quantiles4)

## plot 5% quantiles
quantiles0 = apply(I0, 2, function(x){
  quantile(x, 0.05, na.rm=T)
})
quantiles4 = apply(I4, 2, function(x){
  quantile(x, 0.05, na.rm=T)
})

plot(log10(quantiles0), pch=16, col="black", ylim = c(-1,3),
     main = "log10 5% quantiles of intensities",
     sub = "black: 0 hrs, green: 4 hrs")
points(log10(quantiles4), pch=16, col="limegreen")
abline(v = c(4.5, 8.5, 12.5), lty="dashed")


# I0 = apply(I0, 2, function(x){
#   x = x/quantile(x, 0.05, na.rm=T)
#   return(x)
# })
# 
# I4 = apply(I4, 2, function(x){
#   x = x/quantile(x, 0.05, na.rm=T)
#   return(x)
# })

pdf("results/parse-intensities_diagnostics.pdf", height = 10, width = 10)
diagnostics(I0)
diagnostics(I4)

## impute missing values
I0[!is.finite(I0)] = 0
I0 = apply(I0, 2, function(x){
  x = x+min(x[x!=0])
})

I4[!is.finite(I4)] = 0
I4 = apply(I4, 2, function(x){
  x = x+min(x[x!=0])
})

## fold change + log-transformation
FC = I4/I0 - 1
FC = FC - min(FC, na.rm = T)
FC = log(FC+1)

diagnostics(FC)

## double log
FCl = log(FC+1)
diagnostics(FCl)

dev.off()


### OUTPUT ###
save(I0, file = "data/intensities-0hrs.RData")
write.csv(INT_0hrs, "data/intensity-table-0hrs.csv", row.names = F)

save(I4, file = "data/intensities-4hrs.RData")
write.csv(INT_4hrs, "data/intensity-table-4hrs.csv", row.names = F)

save(FC, file = "data/fold-changes.RData")
save(FCl, file = "data/fold-changes_double-log.RData")


# ----- some thoughts ----
## square root
FCsq = sqrt(FC)
diagnostics(FCsq)

## median polishing
# not really working
mp = medpolish(FC)
FCm = medpolish(FC)$residuals
diagnostics(FCm)

## quantile normalisation
# destroys p-values
# FCq = quantileNorm(FC)
# diagnostics(FCq)

# subtract column medians
# FCs = apply(FC, 2, function(x){
#   x = x-median(x)
# }) 
# diagnostics(FCs)

# ----- plot double-log transformed fold changes -----
library(ggplot2)
library(stringr)

FCl_tbl = tidyr::gather(as.data.frame(FCl))
FCl_tbl$condition = str_extract_all(FCl_tbl$key, "^[:alnum:]+", simplify = T)

s = FCl_tbl %>%
  dplyr::group_by(condition) %>%
  dplyr::summarise(n = dplyr::n(),
                   mean = mean(value),
                   median = median(value),
                   std = sd(value))
print.data.frame(s)

theme_set(theme_classic())
vp = ggplot(FCl_tbl, aes(x = condition, y=value, fill = condition)) +
  geom_violin(size = .1,  alpha = .8) +
  stat_summary(fun = mean, fun.min = mean, fun.max = mean,
               geom = "crossbar", 
               width = .3,
               position = position_dodge(width = 4),
               col = "black") +
  scale_fill_manual(values = c("black", "lightblue", "purple","gray")) +
  ylab("double log-transformed fold changes") +
  theme(legend.position = "none")
vp
ggsave(filename = "~/Documents/Studium/Fachvertiefung+BA/Praktikumsbericht/plots/intensityDistr.ps",
       plot = vp, device = cairo_ps,
       dpi = "retina", height = 3*5, width = 4*5, units = "cm")





