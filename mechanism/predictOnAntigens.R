### karat projetc - PCPS mechanism ###
# description:  predict cleavage/splice sites of model antigens
# input:        antigen sequences, stiff and informative parameters
# output:       predicted cleavage/splicing strengths across antigen sequence
# author:       HPR

library(dplyr)
library(stringr)
source("src/invitroSPI_utils.R")
source("src/_extract-aa.R")
source("src/SCS+PSP-P1.R")

AAchar_here = c("P","G","C","M","A","V","I","L","F","Y","W","H","R","K","D","E","N","Q","S","T","X")
AAchar_here_sorted = sort(AAchar_here)

SRpos = c("P6"=-5,"P5"=-4,"P4"=-3, "P3"=-2, "P2"=-1, "P1"=0,
          "P-1"=1, "P-2"=2, "P-3"=3, "P-4"=4, "P-5"=5, "P-6"=6)

interesting_residues = names(SRpos)
allCombos = tidyr::crossing(interesting_residues,AAchar_here_sorted)
allCombos = do.call(paste, c(allCombos[c("interesting_residues","AAchar_here_sorted")], sep=";"))

pseudo = 1e-05

### INPUT ###
# ovalbumin antigens
antigen = "LELPFASGTMSMLVLLPDEVSGLEQLESIINFEKLTEWTSSNVMEERKIKVYLPRMKMEE"
antigen = "SMLVLLPDEVSGLEQLESIINFEKLTEWTSSNVMEERKIK"
antigenName = "siinfekl"

# gp100
# antigen = "DWLGVSRQLRTKAWNRQLYPEWTEAQRLDC"
# antigen = "VGATKVPRNQDWLGVSRQLRTKAWNRQLYPEWTEAQRLDCWRGGQVSLKV"
antigen = "RRGSRSYVPLAHSSSAFTITDQVPFSVSVSQLRALDGGNKHFLRNQPLTF"

antigen = "AHSSSAFTITDQVPFSVSVSQLRALDGGNK"
antigenName = "gp100_201-230"

antigen = "AHSSSAFTIMDQVPFSVSVSQLRALDGGNK"
antigenName = "gp100_201-230_mutated"

antigen = "ALGVNAENPPAYISSVAYGRQVYLKLSTNSHSTKVKAAFD"
antigenName = "LLO"

antigen = "MAHHHHHHMSDNGPQNQRNAPRITFGGPSDSTGSNQNGERSGARSKQRRPQGLPNNTASWFTALTQHGKEDLKFPRGQGVPINTNSSPDDQIGYYRRATRRIRGGDGKMKDLSPRWYFYYLGTGPEAGLPYGANKDGIIWVATEGALNTPKDHIGTRNPANNAAIVLQLPQGTTLPKGFYAEGSRGGSQASSRSSSRSRNSSRNSTPGSSRGTSPARMAGNGGDAALALLLLDRLNQLESKMSGKGQQQQGQTVTKKSAAEASKKPRQKRTATKAYNVTQAFGRRGPEQTQGNFGDQELIRQGTDYKHWPQIAQFAPSASAFFGMSRIGMEVTPSGTWLTYTGAIKLDDKDPNFKDQVILLNKHIDAYKTFPPTEPKKDKKKKADETQALPQRQKKQQTVTLLPAADLDDFSKQLQQSMSSADSTQA"
antigenName = "SARSCoV2-nucleoprotein"

load("results/Bayesian_ProteaSMM/PLOTS/LOV/0622_PCPposteriors_stiff+informative.RData")
pcp_parameters = parameters
load("results/Bayesian_ProteaSMM/PLOTS/LOV/0623_PSPposteriors_stiff+informative.RData")
psp_parameters = parameters


### MAIN PART ###
suppressWarnings(dir.create("results/Bayesian_ProteaSMM/predictOnAntigens/"))

# ----- get substrate matrix -----
# for only one SR
getSubstrateCounts = function(Seq, paramNames) {
  
  pos = data.frame(substrateSeq = Seq, residue = c(1:nchar(Seq)))
  SRTBL = sapply(SRpos, function(x){
    substr(pos$substrateSeq, start = pos$residue+x, stop = pos$residue+x)
  }) %>%
    as.data.frame()
  SRTBL[SRTBL == ""] = "X"
  
  master = matrix(0, nrow = nrow(SRTBL), ncol = length(allCombos))
  colnames(master) = allCombos
  
  for (j in 1:nrow(SRTBL)) {
    cntN = paste(colnames(SRTBL),SRTBL[j,],sep = ";") %>% as.character()
    master[j,cntN] = 1
  }
  
  master = master[, paramNames]
  
  return(master)
}

Xpsp = getSubstrateCounts(Seq = antigen, paramNames = colnames(psp_parameters))
Xpcp = getSubstrateCounts(Seq = antigen, paramNames = colnames(pcp_parameters))

# ----- predict SCS- and PSP-P1 -----

predictSCSandPSP = function(X, p) {
  
  X = X[,colnames(p)]
  
  tsim = apply(p,1,function(j){
    X%*%log(j)
  })
  
  t = (exp(tsim) - pseudo)*100
  
  tsims_mean = apply(t,1,mean,na.rm = T)
  tsims_sd = apply(t,1,sd,na.rm = T)
 
  return(list(mean = tsims_mean,
              sd = tsims_sd)) 
}

pspPred = predictSCSandPSP(X = Xpsp, p = psp_parameters)
pcpPred = predictSCSandPSP(X = Xpcp, p = pcp_parameters)

# ----- plot predictions ------
# SARS-Cov2 nucleoprotein for thesis
x = readxl::read_excel("data/SARS-CoV2_validation/SARS-COV2-Nucleoprotein_PCP-epitopes_270522.xlsx", sheet = 2)
sarspeps = x$`Peptide Sequence` %>% na.omit() %>% unique()
loci = str_locate_all(antigen, pattern = sarspeps) %>% plyr::ldply() %>% as.data.frame() %>% arrange(start)

p = data.frame(residue = c(c(1:nchar(antigen)), c(1:nchar(antigen))+.2),
               p1 = c(pcpPred$mean, pspPred$mean),
               stdev = c(pcpPred$sd, pspPred$sd),
               col = c(rep(plottingCols["PCP"], nchar(antigen)),
                       rep(plottingCols["PSP"], nchar(antigen))))
yl = max(p$p1+p$stdev) %>% ceiling()

pdf(paste0("results/Bayesian_ProteaSMM/predictOnAntigens/",antigenName,".pdf"), height = 6, width = 15)

# par(mfrow = c(2,2))
start = 1
ends = c(100,200,295,nchar(antigen))

for (j in ends) {
  cntloci = loci[intersect(which(loci$start >= start),which(loci$end <= j)),]
  
  plot(p1~residue, data=p,
       type = "h", lwd = 1, col = col,
       xlim = c(start,j),
       ylim = c(-12, yl),
       main = antigenName,
       xlab = "",
       ylab = "predicted cleavage/splicing strength (%)",
       axes=F)
  arrows(p$residue, p$p1-p$stdev, p$residue, p$p1+p$stdev,
         length=0.03, angle=90, code=3,
         lty = "solid", lwd = .3, col = p$col) %>%
    suppressWarnings()
  
  # positions
  # xlabel = paste(antigen %>% strsplit("") %>% unlist(), seq(1, nchar(antigen)), sep = "")[start:j]
  # xlabel = seq(1, nchar(antigen))[start:j]
  # text(x = seq(start, j), par("usr")[3]-.3, labels = xlabel,
  #      # srt=90,
  #      xpd=T, cex=.5)
  axis(1)
  axis(2, at = seq(0,round(yl,-1),20))
  
  legend("topleft",
         legend = c("cleavage strength (SCS-P1)", "splicing strength (PSP-P1)"),
         lty = c("solid","solid"), col = c(plottingCols["PCP"],plottingCols["PSP"]),
         lwd = c(1,1),
         cex = .6, bty = "n", horiz = F)
  
  y = -1
  cnt.intercept = -1
  
  for (i in 1:nrow(cntloci)) {
    
    segments(x0 = as.numeric(cntloci$start[i])+.1, y0 = cnt.intercept, x1 = as.numeric(cntloci$end[i])+.1, y1 = cnt.intercept,
             lwd = .5, col = "black")
    cnt.intercept = cnt.intercept-.5
    
  }
  
  start = j+1
  
}


dev.off()




# SIINFEKL for thesis
yl = 50
png(paste0("results/Bayesian_ProteaSMM/predictOnAntigens/",antigenName,".png"), height = 5, width = 10, units = "in", res = 300)
plot(p1~residue, data=p,
     type = "h", lwd = 4, col = col,
     ylim = c(-yl, yl),
     main = antigenName,
     xlab = "",
     ylab = "predicted cleavage/splicing strength (%)",
     axes=F)
arrows(p$residue, p$p1-p$stdev, p$residue, p$p1+p$stdev,
       length=0.03, angle=90, code=3,
       lty = "solid", lwd = 1, col = p$col) %>%
  suppressWarnings()
# positions
# xlabel = c("", paste(antigen %>% strsplit("") %>% unlist(), seq(1, nchar(antigen)), sep = ""))
# text(x = seq(0, nchar(antigen)), par("usr")[3]-.3, labels = xlabel,
#      srt=90, xpd=T, cex=1)

xlabel = c("", paste(antigen %>% strsplit("") %>% unlist()))
font = rep(1,length(xlabel))
font[str_locate(antigen,"SIINFEKL")[1]:str_locate(antigen,"SIINFEKL")[2]] = 2
text(x = seq(1, nchar(antigen)+1), par("usr")[3]-.3, labels = paste(antigen %>% strsplit("") %>% unlist()),
     xpd = T, cex = 1, font = font)
axis(2, at = seq(round(-yl,-1),round(yl,-1),20))

abline(v = str_locate(antigen,"SIINFEKL")[1]-.5, lty = "dashed")
abline(v = str_locate(antigen,"SIINFEKL")[2]+.5, lty = "dashed")

legend("topleft",
       legend = c("cleavage strength (SCS-P1)", "splicing strength (PSP-P1)"),
       lty = c("solid","solid"), col = c(plottingCols["PCP"],plottingCols["PSP"]),
       lwd = c(4,4),
       cex = .8, bty = "n", horiz = F)
dev.off()


### OUTPUT ###

