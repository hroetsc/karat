### karat projetc - PCPS mechanism ###
# description:  correlate predicted cleavage/splicing strength to peptide MS1 intensity
# input:        quantitative polypeptide dataset, posteriors
# output:       predicted peptide score, correlation with MS1 intensity
# author:       HPR

library(dplyr)
library(stringr)
source("src/invitroSPI_utils.R")
source("src/SCS+PSP-P1.R")
source("src/_extract-aa.R")


### INPUT ###
load("data/aSPIre.RData")
load("results/Bayesian_ProteaSMM/PLOTS/LOV/0622_PCPposteriors_stiff+informative.RData")
load("data/ProteaSMM/PCP_SR1extfeat_P1/stiff_informative_params.RData")


### MAIN PART ###
# ----- get median scores -----
paramNames = colnames(parameters)
scores = data.frame(position = str_split_fixed(paramNames,";",Inf)[,1],
                    aa = str_split_fixed(paramNames,";",Inf)[,2],
                    score = apply(parameters,2,median))


# ----- preprocessing -----
# actually, PCPs need to be filtered out

Quant = Kinetics %>%
  disentangleMultimappers.Type() %>%
  tidyr::separate_rows(digestTimes, intensities, sep=";") %>%
  rename(digestTime = digestTimes,
         intensity = intensities) %>%
  mutate(intensity = as.numeric(intensity),
         digestTime = as.numeric(digestTime)) %>%
  tidyr::separate_rows(positions, sep = ";")


# ----- extract amino acids -----

DB = uniquePeptides(Quant) %>% select(-biological_replicate, -digestTime, -intensity)
DBaa = extract_aminoacids(DB, onlyValidSeq = F)

DBall = left_join(DB,DBaa)

# keep only the amino acids that occur in the actual peptide (?)
# + the intervening sequence length? is this valid?

DBall = DBaa %>%
  tidyr::gather(position, aa, -spliceType, -positions, -pepSeq, -substrateID, -substrateSeq)
DBall$aa[DBall$aa == ""] = "X"

# join with scores
DBall = left_join(DBall, scores) %>%
  na.omit()
DBallScores = DBall %>%
  group_by(substrateID, spliceType, pepSeq) %>%
  summarise(overallScore = sum(score, na.rm = T))


# ----- correlate with intensity fold change -----

QuantAll = inner_join(Quant, DBallScores)
ALLfc = QuantAll %>%
  group_by(substrateID, pepSeq, biological_replicate) %>%
  summarise(FC = (intensity[digestTime == 4]+1e-05)/(intensity[digestTime == 1]+1e-05) - 1)

QuantAll = left_join(QuantAll, ALLfc)
QuantAll = QuantAll[is.finite(QuantAll$FC), ]

k = which(QuantAll$digestTime == 4 & QuantAll$spliceType == "PCP")
cor(log10(QuantAll$FC[k]), QuantAll$overallScore[k])
plot(log10(QuantAll$FC[k]), QuantAll$overallScore[k], pch = 16)



