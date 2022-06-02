### INHIBITOR KINETICS ###
# description:  parse aSPIre output for DE analysis
# input:        aSPIre filtered intensities (synthesis errors removed) + final kinetics
# output:       table with intensities for all peptides
# author:       HR

library(dplyr)
library(stringr)
source("../../brainstorming/src/invitroSPI_utils.R")


### INPUT ###
finalKinetics = read.csv("aSPIre_manual/results/TSN5inhibitor/finalKinetics.csv", stringsAsFactors = F)
load("aSPIre_manual/results/TSN5inhibitor/QUANTITIES_filtered.Data")
sample_list = read.csv("invitroSPI/INPUT/sample_list.csv", stringsAsFactors = F)


### MAIN PART ###
# ----- preprocessing & formatting -----
DATA = QUANTfiltered

DATA$spliceType[is.na(DATA$spliceType)] = "PCP"
DATA = na.omit(DATA)
DATA = DATA %>%
  disentangleMultimappers.Type()

KU = DATA %>%
  distinct(substrateID, pepSeq, .keep_all = T)
table(KU$substrateID,KU$spliceType)

# add full replicate info
DATA$condition = paste0(DATA$biological_replicate, "_tech", DATA$technical_replicate)
DATA = DATA[order(DATA$pepSeq), ]


# ----- fetch intensities -----
NuniquePeps = DATA$pepSeq %>% unique() %>% length()
print(NuniquePeps)

I0 = DATA %>%
  tidyr::separate_rows(times_pasted, int_pasted, sep = ";") %>%
  mutate(digestTime = as.numeric(times_pasted),
         intensity = as.numeric(int_pasted)) %>%
  filter(digestTime == 0) %>%
  select(pepSeq, condition, intensity) %>%
  tidyr::spread(condition, -pepSeq) %>%
  tibble::column_to_rownames("pepSeq") %>%
  as.matrix()

INT_4hrs = DATA %>%
  tidyr::separate_rows(times_pasted, int_pasted, sep = ";") %>%
  mutate(digestTime = as.numeric(times_pasted),
         intensity = as.numeric(int_pasted)) %>%
  filter(digestTime == 4) %>%
  select(pepSeq, spliceType, positions, condition, intensity) %>%
  tidyr::spread(condition, intensity)

I4 = INT_4hrs %>%
  select(-spliceType, -positions) %>%
  tibble::column_to_rownames("pepSeq") %>%
  as.matrix()


# ----- remove synthesis errors -----
# for no inhibitor: remove everything that has a lower intensity at 4 hrs than at 0 hrs
# mean over technical and biological replicates
colnames(I0) == colnames(I4)

no_idx = c(13:16)
b1_idx = c(1:4)
b2_idx = c(5:8)
b5_idx = c(9:12)

rem = which(rowMeans(I4[,no_idx]) < rowMeans(I0[,no_idx]))
length(rem) / nrow(I4)

I0 = I0[-rem,]
I4 = I4[-rem,]
INT_4hrs = INT_4hrs[-rem,]

table(INT_4hrs$spliceType)

# ---- statistics -----

apply(I0, 2, summary)
apply(I4, 2, summary)

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


## quantile normalisation over all time points/ technical replicates
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

plot(log10(quantiles0), pch=16, col="black", ylim = c(0,8),
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

plot(log10(quantiles0), pch=16, col="black", ylim = c(-1,5),
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
save(I4, file = "data/intensities-4hrs.RData")

save(FC, file = "data/fold-changes.RData")
save(FCl, file = "data/fold-changes_double-log.RData")

write.csv(INT_4hrs, "data/intensity-table-4hrs.csv", row.names = F)

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
ggsave(filename = "results/intensityDistr.ps",
       plot = vp, device = cairo_ps,
       dpi = "retina", height = 3*3, width = 4*3, units = "cm")





