### INHIBITOR KINETICS ###
# description:  get substrate amount
# input:        XICs for all replicates, 0 and 4 hrs, charge +3 and +4
# output:       substrate degradation
# author:       HR

library(dplyr)
library(stringr)
library(ggplot2)
library(zoo)

### INPUT ###
fs = list.files("data/substrate_degradation", pattern = ".txt",
                full.names = T, recursive = T)
fs


### MAIN PART ###
# ----- preprocessing & formatting -----
# get conditions
cond = data.frame(sample =  c("XA", "YA", "XB", "YB",
                              "XC", "YC","XD", "YD"),
                  condition = c("noInh_bio1", "noInh_bio2",
                                "b5_bio1", "b5_bio2",
                                "b2_bio1", "b2_bio2",
                                "b1_bio1", "b1_bio2"),
                  biological_replicate = seq(1,8))


# ------ extract XICs -----
XICs = lapply(fs, read.table, header=F)
names(XICs) = fs
XICs = plyr::ldply(XICs)

names(XICs) = c("filename", "RT", "intensity")
XICs$sampleName = str_extract_all(basename(XICs$filename), "(?<=XIC_)[:alpha:]{2}[:digit:]{1}",simplify = T)
XICs$sample = str_extract_all(XICs$sampleName, "[:alpha:]{2}", simplify = T)
XICs$digestTime = sapply(XICs$sampleName, function(x){
  if (grepl("1", x)) {
    0
  } else {
    4
  }
})

MASTER = left_join(XICs,cond)
MASTER$techRep = str_extract_all(basename(MASTER$filename),"(?<=rep)[:digit:]{1}",simplify = T)
MASTER$charge = str_extract_all(basename(MASTER$filename),"(?<=charge)[:digit:]{1}",simplify = T)


# ----- plot XICs ------
MASTER$digestTime = factor(MASTER$digestTime, levels = c("0","4"))
theme_set(theme_classic())

xic_plot = ggplot(MASTER, aes(x = RT, y=intensity, col=digestTime)) +
  geom_line(aes(y=rollmean(intensity, 10, na.pad = T))) +
  scale_color_manual("digestion time", values = c("midnightblue", "mediumseagreen")) +
  ylim(c(0,1.5e10)) +
  xlim(c(2700,3800)) +
  ggtitle("XICs") + 
  xlab("RT [s]") + 
  ylab("intensity")

xic_plot
xic_grid = xic_plot + facet_wrap(~condition, ncol = 2, nrow=4)
xic_grid

xic_grid_large = xic_plot + facet_wrap(condition~charge, ncol = 4, nrow = 4)
xic_grid_large

ggsave(filename = "results/XICs_substrate.png", plot = xic_grid,
       device = "png", height = 10, width = 6, dpi = "retina")
ggsave(filename = "results/XICs_substrate_large.png", plot = xic_grid_large,
       device = "png", height = 10, width = 12, dpi = "retina")


# ------ get substrate degradation -----
# mean over technical replicates for both charges
tmp = MASTER %>%
  dplyr::group_by(condition, digestTime, RT, charge) %>%
  dplyr::summarise(mu = mean(intensity))

# define range of integration
integrated = tmp %>%
  dplyr::filter(RT>2700 & RT<3800 & mu > 1e05) %>%
  dplyr::group_by(condition, digestTime) %>%
  dplyr::summarise(AUC = sum(mu))

# get actual degradation
degraded = integrated %>%
  dplyr::group_by(condition) %>%
  dplyr::summarise(degraded = -1*diff(AUC))

degraded
# cond = unique(integrated$condition)
# sapply(cond, function(x) {
#   integrated$AUC[integrated$condition == x & integrated$digestTime == 4] - integrated$AUC[integrated$condition == x & integrated$digestTime == 0]
# })

### OUTPUT ###
save(integrated, file = "data/substrate-intensities.RData")
save(degraded, file = "data/substrate-degradation.RData")


