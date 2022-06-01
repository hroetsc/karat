### protein PCPS ###
# description:  validate old and new approach to generate random databases
# input:        random database from invitroSPI manuscript, random database with new approach
# output:       random DB comparison
# author:       HPR


library(dplyr)
library(stringr)
library(ggplot2)
source("src/invitroSPI_utils.R")
source("src/plotting_utils.R")

theme_set(theme_classic())

### INPUT ###
# load("../../invitroSPI/revision/data/wholeDB_random/wholeDB_random_AAMultimapper.RData")
load("data/randomDB.RData")
# old_randomDB = rbind(randomDB.whole_processed_AA$cis,randomDB.whole_processed_AA$revCis, randomDB.whole_processed_AA$trans)
new_randomDB = rndDB
load("data/randomDB_smart.RData")
smart_randomDB = rndDB

### MAIN PART ###
# ----- preprocessing -----
Subs = intersect(old_randomDB$substrateID %>% unique(),
                 new_randomDB$substrateID %>% unique())

old_randomDB2 = old_randomDB[old_randomDB$substrateID %in% Subs, ]
new_randomDB2 = new_randomDB[new_randomDB$substrateID %in% Subs, ] %>%
  disentangleMultimappers.Type()

old_randomDB2$productType %>% unique()
table(old_randomDB2$spliceType)
table(new_randomDB2$spliceType)


# ----- violin plot -----
splittedViolinPlot_comp = function(data) {
  
  data$db = factor(data$db, levels = c("PCPcomb", "samplePos"))
  
  # print stats
  s = data %>%
    dplyr::group_by(type, db) %>%
    dplyr::summarise(n = dplyr::n(),
                     mean = mean(value),
                     median = median(value),
                     std = sd(value))
  print.data.frame(s)
  
  theme_set(theme_classic())
  
  # plotting
  sc = ggplot(data = data, aes(x = type, y = value, fill = db)) +
    geom_split_violin(size = .25) +
    stat_summary(fun = median, fun.min = median, fun.max = median,
                 geom = "crossbar", 
                 width = .2,
                 position = position_dodge(width = .2)) +
    xlab("product type") +
    scale_fill_manual("database",
                      values = c(plottingCols[["proteins"]], plottingCols[["polypeptides"]]),
                      drop = F) +
    theme(axis.text.x = element_text(angle = 90))
  
  # add labs and title downstream
  
  return(sc)
}


# ----- peptide length -----

PepLength = function(old, new) {
  
  old = removeMultimappers.Type(old)
  new = removeMultimappers.Type(new)
  
  pl = data.frame(value = nchar(old$pepSeq),
                  substrateID = old$substrateID,
                  type = old$spliceType,
                  db = "PCPcomb")
  pl = rbind(pl,
             data.frame(value = nchar(new$pepSeq),
                        substrateID = new$substrateID,
                        type = new$spliceType,
                        db = "samplePos")) %>%
    na.omit()
  
  pl$type = factor(pl$type, levels = c("PCP","cis", "revCis", "trans"))
  
  pepLen = splittedViolinPlot_comp(data = pl) +
    ylab("peptide product length (aa residues)") +
    ggtitle("peptide length distribution")
  
  ggsave(filename = paste0("results/validation/_randomDBValidation_compareApproaches.png"),
         plot = pepLen, dpi = "retina",
         height = 6, width = 6)
}

# PepLength(old_randomDB2, new_randomDB2[new_randomDB2$spliceType != "PCP", ])
PepLength(new_randomDB, smart_randomDB)



### OUTPUT ###

