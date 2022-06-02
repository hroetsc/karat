### INHIBITOR KINETICS ###
# description:  integrated cleavage/splicing strength for each synthetic peptide candidate
# input:        synthetic peptides, overall and b5-specific cleavage/splicing strength
# output:       SCS + PSP-P1 for each candidate --> for b5 (LLVY) and all other subunits
# author:       HR

library(dplyr)
library(stringr)
library(readxl)
library(ggplot2)
source("../../brainstorming/src/invitroSPI_utils.R")

theme_set(theme_classic())

### INPUT ###
load("data/SCS+PSP-P1_TS5all.RData")
load("data/SCS+PSP-P1_TSN5b5specific.RData")
synPeps = read_excel("../../3_competitor_assays/syntheticPeps_220428.xlsx", range = "A1:E11")

### MAIN PART ###
# ----- get SCS and PSP for each candidate -----
b5 = synPeps %>%
  mutate(start = str_extract_all(position_in_substrate, "^[:digit:]+", simplify = T) %>% as.numeric(),
         end = str_extract_all(position_in_substrate, "[:digit:]+$", simplify = T) %>% as.numeric())
b5 = b5 %>%
  mutate(SCS_P1 = sapply(b5$end, function(x){
    TSN5b5specific$scs_mean[TSN5b5specific$residue == x]
  }),
  PSP_P1 = sapply(b5$end, function(x){
    TSN5b5specific$psp_mean[TSN5b5specific$residue == x]
  }),
  SCS_integr = apply(b5,1,function(x){
    TSN5b5specific$scs_mean[which(TSN5b5specific$residue == as.numeric(x[["start"]])):which(TSN5b5specific$residue == as.numeric(x[["end"]]))] %>% sum()
  }),
  PSP_integr = apply(b5,1,function(x){
    TSN5b5specific$psp_mean[which(TSN5b5specific$residue == as.numeric(x[["start"]])):which(TSN5b5specific$residue == as.numeric(x[["end"]]))] %>% sum()
  }))


all = synPeps %>%
  mutate(start = str_extract_all(position_in_substrate, "^[:digit:]+", simplify = T) %>% as.numeric(),
         end = str_extract_all(position_in_substrate, "[:digit:]+$", simplify = T) %>% as.numeric())
all = all %>%
  mutate(SCS_P1 = sapply(all$end, function(x){
    TSN5all$scs_mean[TSN5all$residue == x]
  }),
  PSP_P1 = sapply(all$end, function(x){
    TSN5all$psp_mean[TSN5all$residue == x]
  }),
  SCS_integr = apply(all,1,function(x){
    TSN5all$scs_mean[which(TSN5all$residue == as.numeric(x[["start"]])):which(TSN5all$residue == as.numeric(x[["end"]]))] %>% sum()
  }),
  PSP_integr = apply(all,1,function(x){
    TSN5all$psp_mean[which(TSN5all$residue == as.numeric(x[["start"]])):which(TSN5all$residue == as.numeric(x[["end"]]))] %>% sum()
  }))


# ----- visualisation -----

MASTER = rbind(b5 %>% mutate(strength = "b5_specific"),
               all %>% mutate(strength = "all_peptides")) %>%
  mutate(P1 = paste0(substr(sequence,nchar(sequence),nchar(sequence)), end)) %>%
  select(substrateID,strength,site,start,end,P1,SCS_P1,PSP_P1,SCS_integr,PSP_integr) %>%
  mutate(Delta_P1 = PSP_P1-SCS_P1,
         Delta_integr = PSP_integr-SCS_integr,
         PSP_P1 = -1*PSP_P1,
         PSP_integr = -1*PSP_integr)
save(MASTER, file = "results/CANDIDATES_scs+psp.RData")

MASTER = MASTER %>%
  tidyr::gather(type,value, -substrateID,-site,-start,-end,-strength,-P1)

p1 = MASTER %>%
  filter(type %in% c("SCS_P1","PSP_P1")) %>%
  ggplot(aes(x=reorder(P1,end),y=value, group = strength, alpha=strength, fill=type)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c(plottingCols[["PSP"]],plottingCols[["PCP"]]),
                    labels = c("splicing strength", "cleavage strength")) +
  scale_alpha_manual(values=c(0.5, 1)) +
  ylim(c(-65,65)) +
  ylab("cleavage/splicing strength (%)") +
  xlab("P1 of synthetic peptide") +
  ggtitle("SCS-/PSP-P1")

integr = MASTER %>%
  filter(type %in% c("SCS_integr","PSP_integr")) %>%
  ggplot(aes(x=reorder(P1,end),y=value, group = strength, alpha=strength, fill=type)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c(plottingCols[["PSP"]],plottingCols[["PCP"]]),
                    labels = c("splicing strength", "cleavage strength")) +
  scale_alpha_manual(values=c(0.5, 1)) +
  ylim(c(-65,65)) +
  ylab("cleavage/splicing strength (%)") +
  xlab("P1 of synthetic peptide") +
  ggtitle("integrated SCS/PSP")

delta = MASTER %>%
  filter(type %in% c("Delta_integr")) %>%
  ggplot(aes(x=reorder(P1,end),y=value, group = strength, alpha=strength, fill=type)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("orchid4","hotpink")) +
  scale_alpha_manual(values=c(0.4, 0.8)) +
  ylim(c(-30,30)) +
  ylab("splicing - cleavage strength (%)") +
  xlab("P1 of synthetic peptide") +
  ggtitle("difference between PSP and SCS")
delta

both = gridExtra::grid.arrange(p1,integr,delta,ncol = 1)

ggsave(filename = "results/CANDIDATES_scs+psp.png", plot = both,
       height = 15, width = 10, dpi = "retina")

