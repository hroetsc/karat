

library(stringr)
library(dplyr)

load("results/SRspecificity/DBproc.RData")

DBproc$sr1[which(str_detect(DBproc$sr1, "G$") & nchar(DBproc$sr1) == 5)] %>% unique()
DBproc$sr1[which(str_detect(DBproc$sr1, "P$") & nchar(DBproc$sr1) == 5)] %>% unique()


DBproc$sr1[which(str_detect(DBproc$sr1, "F$") & nchar(DBproc$sr1) == 5)] %>% unique()
DBproc$sr1[which(str_detect(DBproc$sr1, "Y$") & nchar(DBproc$sr1) == 5)] %>% unique()

DBproc$sr1 %>% unique() %>% length()
DBproc$sr1 %>% unique() %>% length()
