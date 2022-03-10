### INHIBITOR KINETICS ###
# description:  parse qiSPI output for DE analysis
# input:        qiSPI output
# output:       table with 4 hrs intensities for all peptides
#               at all conditions / bio reps
# author:       HR

library(dplyr)

bioReps = 8

### INPUT ###
finalKinetics = read.csv("qiSPI/OUTPUT/TSN5_0+4/finalKinetics.csv",
                         stringsAsFactors = F)


### MAIN PART ###
# ----- stats -----
NuniquePeps = finalKinetics$pepSeq %>% unique() %>% length()
paste0("unique peptides identified: ", NuniquePeps) %>% print()

# ----- map conditions to bioReps -----

mping = c("noInh_rep1" = 1, "noInh_rep2" = 2,
          "b1_rep1" = 7, "b1_rep2" = 8,
          "b2_rep1" = 5, "b2_rep2" = 6,
          "b5_rep1" = 3, "b5_rep2" = 4)

finalKinetics$condition = sapply(finalKinetics$biological_replicate, function(x){
  names(mping)[mping == x]
})

# ----- summarise intensities -----
# 4 hrs intensity
plot(density(log10(finalKinetics$tp_4)))

# difference between 0 and 4 hrs intensity
finalKinetics$diff = finalKinetics$tp_4 - finalKinetics$tp_0

plot(density(log10(finalKinetics$diff[finalKinetics$diff>0])))
lines(density(log10(abs(finalKinetics$diff[finalKinetics$diff<0]))),
      col = "red")

summary(finalKinetics$diff)
quantile(finalKinetics$diff, c(0.05, 0.5, 0.9, 0.95))

# ----- data formatting -----
INTtable = data.frame(substrateID = rep(NA, nrow(finalKinetics)),
                      substrateSeq = NA,
                      spliceType = NA,
                      positions = NA,
                      pepSeq = NA,
                      noInh_rep1=NA,noInh_rep2=NA,b1_rep1=NA,b1_rep2=NA,
                      b2_rep1=NA,b2_rep2=NA,b5_rep1=NA,b5_rep2=NA)

k = 1
counter = 1
while (k <= NuniquePeps) {
  
  cnt = finalKinetics[c(counter):c(counter+bioReps-1), ]
  f = all(cnt$biological_replicate == seq(1,bioReps))
  if(!f) {print("SOMETHING IS WRONG")}
  
  # sequence info
  INTtable$substrateID[k] = cnt$substrateID[1]
  INTtable$substrateSeq[k] = cnt$substrateSeq[1]
  INTtable$spliceType[k] = cnt$spliceType[1]
  INTtable$positions[k] = cnt$positions[1]
  INTtable$pepSeq[k] = cnt$pepSeq[1]
  
  # intensities
  # INTtable[k, c(6:13)] = cnt$diff[match(cnt$condition, names(INTtable)[6:13])]
  # INTtable[k, c(6:13)] = cnt$tp_4[match(cnt$condition, names(INTtable)[6:13])]
  INTtable[k, c(6:13)] = log2(cnt$tp_4[match(cnt$condition, names(INTtable)[6:13])] / (cnt$tp_0[match(cnt$condition, names(INTtable)[6:13])]+1))
  
  k = k+1
  counter = counter+bioReps
}

INTtable = INTtable[1:c(k-1),]

# sanity check
NuniquePeps == INTtable$pepSeq %>% unique() %>% length()
apply(INTtable, 2, function(x){any(is.na(x))})

INTtable$spliceType[is.na(INTtable$spliceType)] = "PCP"

# remove substrate
k2 = which(INTtable$pepSeq == INTtable$substrateSeq)
if (length(k2) > 0) {
  INTtable = INTtable[-k2,]
}

# ----- statistics -----

log(INTtable$noInh_rep1+1) %>% density() %>% plot()
log(INTtable$noInh_rep2+1) %>% density() %>% lines(col = "gray")
log(INTtable$b1_rep1+1) %>% density() %>% lines(col = "brown")
log(INTtable$b1_rep2+1) %>% density() %>% lines(col = "orange")
log(INTtable$b2_rep1+1) %>% density() %>% lines(col = "darkblue")
log(INTtable$b2_rep2+1) %>% density() %>% lines(col = "lightblue")
log(INTtable$b5_rep1+1) %>% density() %>% lines(col = "purple")
log(INTtable$b5_rep2+1) %>% density() %>% lines(col = "hotpink")

intIdx = seq(6,13)

k3 = which(rowSums(INTtable[,intIdx]) == 0)
INTtable = INTtable[-k3, ]



### OUTPUT ###
write.csv(INTtable, "data/intensity_table.csv", row.names = F)

