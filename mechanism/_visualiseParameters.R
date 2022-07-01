### karat projetc - PCPS mechanism ###
# description:  visualise stiff and informative parameters for PCPs/PSPs
# input:        DATA.RData, vector of stiff + informative parameters, posteriors
# output:       parameter visualisation
# author:       HPR

library(dplyr)
library(stringr)
library(eulerr)
library(ggplot2)
library(ggseqlogo)
library(RColorBrewer)
source("src/invitroSPI_utils.R")

rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
rcol <- rf(100)
theme_set(theme_classic())

AAchar_here = c("P","G","C","M","A","V","I","L","F","Y","W","H","R","K","D","E","N","Q","S","T","X")
AAchar_here_sorted = sort(AAchar_here)

# create colour scheme
cs = make_col_scheme(chars = c("P","G","C",# special case
                               "M", "A","V","I","L", # hydrophobic
                               "F","Y","W",  # aromatic
                               "H","R","K",  # basic
                               "D","E",  # acidic
                               "N","Q","S","T"),  # nucleophilic/polar
                     cols = c(rep("deeppink",3),
                              rep("orange", 5),
                              rep("mediumpurple", 3),
                              rep("seagreen3",3),
                              rep("firebrick",2),
                              rep("dodgerblue",4)),
                     groups = c(rep("special",3),
                                rep("hydrophobic", 5),
                                rep("aromatic", 3),
                                rep("basic",3),
                                rep("acidic",2),
                                rep("polar",4)))


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


# ----- heatmaps -----
interesting_residues = c("P6","P5","P4", "P3", "P2", "P1", "P-1","P-2","P-3","P-4","P-5","P-6")


getLogos = function(prodType) {
  png(paste0("results/Bayesian_ProteaSMM/PLOTS/heatMaps",prodType,".png"), height = 8, width = 10, res = 300, units = "in")
  # get frequency matrix
  freqMat = bothR %>%
    filter(reactant == prodType) %>%
    group_by(position, aa) %>%
    summarise(med = median(value))
  
  c = log(median(exp(freqMat$med)))
  freqMat$med = freqMat$med - c
  
  DFt = freqMat %>%
    tidyr::spread(position,med) %>%
    tibble::column_to_rownames("aa")
  
  
  # mising positions / aa
  lack = interesting_residues[!interesting_residues %in% colnames(DFt)]
  if (length(lack) > 0) {
    DFt = cbind(DFt, matrix(NA, nrow(DFt), length(lack)))
    colnames(DFt)[(ncol(DFt)-length(lack)+1):ncol(DFt)] = lack
  }
  lack = AAchar_here[!AAchar_here %in% rownames(DFt)]
  if (length(lack) > 0) {
    DFt = rbind(as.matrix(DFt), matrix(NA, length(lack), ncol(DFt)))
    rownames(DFt)[(nrow(DFt)-length(lack)+1):nrow(DFt)] = lack
  }
  
  
  DFt = as.matrix(DFt[,interesting_residues])
  DFt = DFt[AAchar_here, ]
  
  image(DFt %>% t(), axes = F, col = rcol,
        main = "median posteriors",
        sub = "low: blue, high: red")
  axis(2, at = seq(0,1,1/(length(rownames(DFt))-1)), labels = rownames(DFt))
  axis(1, at = seq(0,1,1/(length(colnames(DFt))-1)), labels = colnames(DFt))
  
  
  # ----- logoplots -----
  
  l1 = ggseqlogo(DFt, method = "custom", seq_type = "aa", font = "helvetica_bold",col_scheme = cs) +
    scale_x_discrete(limits = colnames(DFt)) +
    xlab("position") +
    ylab("median weight") +
    ggtitle(prodType)
  
  l2 = ggseqlogo(DFt-min(DFt,na.rm = T), seq_type = "aa", font = "helvetica_bold",col_scheme = cs) +
    scale_x_discrete(limits = colnames(DFt)) +
    xlab("position") +
    ylab("entropy [bits]") +
    ggtitle(prodType)
  
  dev.off()
  
  return(list(logo1 = l1, logo2 = l2))
}

pcpLogo = getLogos("PCP")
pspLogo = getLogos("PSP")

ggsave(filename = "results/Bayesian_ProteaSMM/PLOTS/logoPlots_1.png",
       plot = gridExtra::grid.arrange(pcpLogo$logo1, pspLogo$logo1),
       height = 8, width = 6, dpi = "retina")
ggsave(filename = "results/Bayesian_ProteaSMM/PLOTS/logoPlots_2.png",
       plot = gridExtra::grid.arrange(pcpLogo$logo2, pspLogo$logo2),
       height = 8, width = 6, dpi = "retina")


# ----- logoplots on entire posterior -----
# not just stiff+informative parameters

load("results/Bayesian_ProteaSMM/PLOTS/LOV/0623_PSP_DeltaH.RData")
pspDeltaH = DeltaH
load("results/Bayesian_ProteaSMM/PLOTS/LOV/0622_PCP_DeltaH.RData")
pcpDeltaH = DeltaH


PCPparam = data.frame(key = names(pcpDeltaH),
                      value = pcpDeltaH) %>%
  mutate(reactant = "PCP")
PSPparam = data.frame(key = names(pspDeltaH),
                      value = pspDeltaH) %>%
  mutate(reactant = "PSP")

bothR = rbind(PCPparam,PSPparam)
bothR$position = str_split_fixed(bothR$key,";",Inf)[,1]
bothR$aa = str_split_fixed(bothR$key,";",Inf)[,2]

# reorder amino acids
bothR$aa = factor(bothR$aa, levels = AAchar_here)

getEntropyLogos = function(prodType) {
  
  # get frequency matrix
  freqMat = bothR %>%
    filter(reactant == prodType) %>%
    group_by(position, aa) %>%
    summarise(med = median(value))
  
  DFt = freqMat %>%
    tidyr::spread(position,med) %>%
    tibble::column_to_rownames("aa")
  DFt = as.matrix(DFt[,interesting_residues])
  DFt = DFt[AAchar_here, ]
  
  # get entropy
  l = ggseqlogo(DFt, method = "custom", seq_type = "aa", font = "helvetica_bold",col_scheme = cs) +
    scale_x_discrete(limits = colnames(DFt)) +
    xlab("position") +
    ylab("Delta H (prior - posterior)") +
    ggtitle(prodType)
  
  return(l)
}

lpsp = getEntropyLogos("PSP")
lpcp = getEntropyLogos("PCP")

ggsave(filename = "results/Bayesian_ProteaSMM/PLOTS/logoPlots_infoGain.png",
       plot = gridExtra::grid.arrange(lpcp, lpsp),
       height = 8, width = 6, dpi = "retina")


