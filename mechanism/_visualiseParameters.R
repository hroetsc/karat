### karat projetc - PCPS mechanism ###
# description:  visualise stiff and informative parameters for PCPs/PSPs
# input:        DATA.RData, vector of stiff + informative parameters, posteriors
# output:       parameter visualisation
# author:       HPR

library(dplyr)
library(stringr)
library(eulerr)
library(ggplot2)
source("src/invitroSPI_utils.R")

theme_set(theme_classic())
AAchar_here = c("P","G","C","M","A","V","I","L","F","Y","W","H","R","K","D","E","N","Q","S","T","X")
AAchar_here_sorted = sort(AAchar_here)


### INPUT ###
load("results/Bayesian_ProteaSMM/PLOTS/LOV/0622_PCPposteriors_stiff+informative.RData")
pcp_parameters = parameters
load("results/Bayesian_ProteaSMM/PLOTS/LOV/0623_PSPposteriors_stiff+informative.RData")
psp_parameters = parameters


### MAIN PART ###
suppressWarnings(dir.create("results/Bayesian_ProteaSMM/PLOTS/_informativeParams/"))

# ----- overlap between stiff+informative parameters -----
png("results/Bayesian_ProteaSMM/PLOTS/_informativeParams/parameterOverlap.png",
    height = 4, width = 5, units = "in", res = 300)
euler(list(psp = colnames(psp_parameters),
           pcp = colnames(pcp_parameters))) %>% plot(quantities = T)
dev.off()

# ----- get parameter distributions -----
PCPparam = tidyr::gather(pcp_parameters %>% as.data.frame()) %>%
  mutate(reactant = "PCP")
PSPparam = tidyr::gather(psp_parameters %>% as.data.frame()) %>%
  mutate(reactant = "PSP")

bothR = rbind(PCPparam,PSPparam)
bothR$position = str_split_fixed(bothR$key,";",Inf)[,1]
bothR$aa = str_split_fixed(bothR$key,";",Inf)[,2]

# reorder amino acids
bothR$aa = factor(bothR$aa, levels = AAchar_here)


# ----- plot per position -----
# plot per position
allP = list()
counter = 1
pos = c("P6","P5","P4", "P3", "P2", "P1", "P-6", "P-5", "P-4", "P-3","P-2", "P-1")
for (p in 1:length(pos)) {
  cntR = bothR[bothR$position == as.character(pos[p]), ]
  
  if (nrow(cntR) > 0) {
    cntP = ggplot(cntR, aes(x = aa, y = value, fill = reactant)) +
      geom_boxplot(outlier.shape = NA, coef = 0, position = position_dodge(0.5), aes(alpha = .8)) +
      scale_fill_manual(values = c(plottingCols[["PCP"]], plottingCols[["PSP"]])) +
      scale_x_discrete(limits = AAchar_here) +
      ylim(c(0,1)) +
      xlab("amino acid") +
      ylab("posterior parameter distribution") +
      ggtitle(pos[p])
    
    allP[[counter]] = cntP
    counter = counter + 1
  }
}

ggsave(filename = "results/Bayesian_ProteaSMM/PLOTS/_informativeParams/PARAMETERS.pdf", 
       plot = gridExtra::marrangeGrob(allP, nrow=6, ncol=2, byrow = F), 
       width = 15, height = 27, dpi = "retina")


# ----- plot per amino acid ----

allPaa = list()
counter = 1
bothR$position = factor(bothR$position, levels = c("P6","P5","P4", "P3", "P2", "P1",
                                                   "P-1","P-2","P-3","P-4","P-5","P-6"))
for (a in 1:length(AAchar_here)) {
  
  cntR = bothR[bothR$aa == as.character(AAchar_here[a]), ]
  
  if (nrow(cntR) > 0) {
    cntP = ggplot(cntR, aes(x = position, y = value, fill = reactant)) +
      geom_boxplot(outlier.shape = NA, coef = 0, position = position_dodge(0.5), aes(alpha = .8)) +
      scale_fill_manual(values = c(plottingCols[["PCP"]], plottingCols[["PSP"]])) +
      scale_x_discrete(limits = c("P6","P5","P4", "P3", "P2", "P1",
                                  "P-1","P-2","P-3","P-4","P-5","P-6")) +
      ylim(c(0,1)) +
      xlab("position") +
      ylab("posterior parameter distribution") +
      ggtitle(AAchar_here[a])
    
    allPaa[[counter]] = cntP
    counter = counter + 1
  }
  
}

ggsave(filename = "results/Bayesian_ProteaSMM/PLOTS/_informativeParams/PARAMETERS_aawise.pdf", 
       plot = gridExtra::marrangeGrob(allPaa, nrow=5, ncol=4, byrow = T), 
       width = 20, height = 25, dpi = "retina")

# ----- number of informative parameters per position -----



