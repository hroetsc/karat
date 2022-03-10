




### b5 against all ###
{
  # create design matrix
  sample_info = data.frame(condition = c(rep("b5_inactive",4), rep("b5_active",12)),
                           replicate = c(c(c("r1","r2","r3","r4"),paste0("r", seq(1,12)))))
  design = model.matrix(~condition, sample_info)
  
  # modify intensity table
  names(QUANT)
  L = QUANT[, c(5:8,1:4,9:16)]
  names(L) = paste0("Intensity ", paste0(sample_info$condition,"_",sample_info$replicate))
}

### b5 against b1 ###
{
  # create design matrix
  sample_info = data.frame(condition = c(rep("b5",4), rep("b1",4)),
                           replicate = rep(c("r1", "r2", "r3", "r4"),2))
  design = model.matrix(~condition, sample_info)
  
  # modify intensity table
  names(QUANT)
  L = QUANT[, c(5:8,13:16)]
  names(L) = paste0("Intensity ", paste0(sample_info$condition,"_",sample_info$replicate))
}

### b5 against b2 ###
{
  # create design matrix
  sample_info = data.frame(condition = c(rep("b5",4), rep("b2",4)),
                           replicate = rep(c("r1", "r2", "r3", "r4"),2))
  design = model.matrix(~condition, sample_info)
  
  # modify intensity table
  names(QUANT)
  L = QUANT[, c(5:8,9:12)]
  names(L) = paste0("Intensity ", paste0(sample_info$condition,"_",sample_info$replicate))
}

### b5 against no inhibitor ###
{
  # create design matrix
  sample_info = data.frame(condition = c(rep("b5",4),rep("noInh",4)),
                           replicate = rep(c("r1", "r2", "r3", "r4"),2))
  design = model.matrix(~condition, sample_info)
  
  # modify intensity table
  names(QUANT)
  L = QUANT[, c(5:8,1:4)]
  names(L) = paste0("Intensity ", paste0(sample_info$condition,"_",sample_info$replicate))
}
















# background correction
# 5% quantile

# transformation
# log2 or VSN (ln)


Qnorm = apply(Q,2,function(x){
  x = x-quantile(x,0.05)
})

Qnorm[Qnorm<=0] = NA
Qnorm %>%
  as.data.frame() %>%
  vis_miss(sort_miss = FALSE, cluster = T)

Qimputed <- Qnorm %>%
  as.data.frame() %>%
  mutate_at(1:ncol(Qnorm),~replace(., . == 0, min(.[.>0], na.rm = TRUE))) %>%
  as.data.frame()
rownames(Qimputed) <- rownames(Q)

Qtransform = log(Qimputed+1)



Qtransform = log2(Q+1)
Qcorr = backgroundCorrect(Qtransform, method="normexp", offset=50)

# Qquantnorm = quantileNorm(Qztr)
tk = medpolish(Qcorr)
Qmedianpolished = tk$residuals

rownames(Qmedianpolished) = INTtable$pepSeq
QUANT = Qmedianpolished %>% as.data.frame()



Qtransform = log2(Qquantnorm+1)
Qmedianpolished = medpolish(Qtransform)$residuals

rownames(Qmedianpolished) = INTtable$pepSeq
QUANT = Qmedianpolished %>% as.data.frame()

## 1. log-transform intensities
# use pseudocounts
Q = apply(INTtable[,intIdx],2,function(x){
  y = log2(x+1)
  return(y)
}) %>%
  as.data.frame()
rownames(Q) = INTtable$pepSeq

## 2. impute missing values
# visualise 0 intensities (missing values)
Q %>%
  replace_with_na_all(condition = ~.x == 0) %>%
  vis_miss(sort_miss = FALSE, cluster = T)

# impute with minimal intensity
Q.minimal <- Q %>%
  mutate_at(1:ncol(Q),~replace(., . == 0, min(.[.>0], na.rm = TRUE))) %>%
  as.data.frame()
rownames(Q.minimal) <- rownames(Q)

# impute with Bayesian PCA
{
  in_matrix <- as.matrix(Q)
  in_matrix[in_matrix==0] <- NA
  # remove rows where all variables are missing
  keep_prot_NA <- rowSums(is.na(in_matrix)) < (ncol(in_matrix)-1)
  in_matrix <- in_matrix[keep_prot_NA,]
  pcaMethods::checkData(in_matrix, verbose=TRUE)
  # BPCA impute
  resBPCA <- in_matrix %>% 
    pcaMethods::pca(method="bpca", center=TRUE, scale = "pareto", nPcs=6)
  plot(resBPCA)
  Q.BPCA <- pcaMethods::completeObs(resBPCA)
  print(resBPCA)
}
{
  x11()
  par(mfrow = c(3,1))
  plotDensities(Q,col = rep(c('red', 'black', "blue", "green"), each = 2), 
                main = 'before missing value imputation', legend = F)
  plotDensities(Q.minimal,col = rep(c('red', 'black', "blue", "green"), each = 2), 
                main = 'missing value imputation with minimum intensity', legend = F)
  plotDensities(Q.BPCA,col = rep(c('red', 'black', "blue", "green"), each = 2), 
                main = 'missing value imputation using BPCA', legend = F)
  dev.copy2pdf(file = "results/missingness.pdf")
}

# Qnorm = quantileNorm(Q.minimal)
# rownames(Qnorm) = rownames(Q.minimal)
# Qnorm = medpolish(Qnorm)$residuals

# # Tukey's median polish
# Qnorm = matrix(NA, nrow(Q), ncol(Q))
# 
# QTk = medpolish(Q.BPCA[,c(1:2)], eps = 0.01, maxiter = 10, trace.iter = TRUE, na.rm = FALSE)
# plot(QTk)
# Qnorm[,c(1:2)] = QTk$overall + outer(QTk$row,QTk$col, "+") + QTk$residuals

# decide with what to continue
# QUANT = Q.BPCA %>% as.data.frame()
# QUANT = Q.minimal %>% as.data.frame()
# QUANT = Q

Qnorm = Q.minimal
Qnorm[,c(1:2)] = quantileNorm(Q.minimal[,c(1:2)])
Qnorm[,c(3:4)] = quantileNorm(Q.minimal[,c(3:4)])
Qnorm[,c(5:6)] = quantileNorm(Q.minimal[,c(5:6)])
Qnorm[,c(7:8)] = quantileNorm(Q.minimal[,c(7:8)])

Qnorm = medpolish(Qnorm)$residuals

QUANT = Qnorm %>% as.data.frame()














# subtract peptide-wise minimal intensity
INTtable[,intIdx] = apply(INTtable[,intIdx],1,function(x){
  k = which.min(x)
  y = x - x[k]
  return(y)
})


# in case of 0/neg values
{
  # set negative values to NA/0
  INTtable[,intIdx] = apply(INTtable[,intIdx], 2, function(x){
    sapply(x, function(y){
      if (y < 0) {
        NA
      } else {y}
    })
  })
  
  k = apply(INTtable[,intIdx],1,function(x){
    all(is.na(x))
  }) %>%
    which()
  
  if (length(k) > 0){
    INTtable = INTtable[-k,]
  }
  
  # impute missing values with minimal intensity across all peptides for a given sample
  # should be 0
  minInt = apply(INTtable[,intIdx], 2, function(x){
    min(x, na.rm = T)
  })
  
  minInt
  # only execute if minInt is 0!
  INTtable[,intIdx] = apply(INTtable[,intIdx], 2, function(x){
    sapply(x, function(y){
      if (is.na(y)) {
        return(0)
      } else {
        return(y)
      }
    })
  })
  
}

# in case value is already intensity difference
# # subtract peptide-wise smallest value
# INTtable[,intIdx] = apply(INTtable[,intIdx],2,function(x){
#   k = which.min(x)
#   y = x - x[k]
#   return(y)
# })


# some ideas
{
  # # subtract sample-wise smallest value
  # INTtable[,intIdx] = apply(INTtable[,intIdx],2,function(x){
  #   k = which.min(x)
  #   y = x - x[k]
  #   return(y)
  # })
  
  # # subtract peptide-wise smallest value
  # INTtable[,intIdx] = apply(INTtable[,intIdx],1,function(x){
  #   k = which.min(x)
  #   y = x - x[k]
  #   return(y)
  # })
}
















k_b1 = which(res$estimate<0 & -log10(res$p.value) > -log10(0.1))
k_b2 = which(res$estimate<0 & -log10(res$p.value) > -log10(0.1))
k_noInh = which(res$estimate<0 & -log10(res$p.value) > -log10(0.1))

# ----- select peptides -----
library(eulerr)

klst = list(b1=k_b1, b2=k_b2, noInh=k_noInh)
euler(klst) %>% plot(quantities=T)

kk = intersect(k_b1, k_b2)
INTtable[kk, ]

write.csv(INTtable[kk, ], "data/intensity_table_selected.csv", row.names = F)

kk2 = c(k_b1, k_b2, k_noInh) %>% unique()
write.csv(INTtable[kk2, ], "data/intensity_table_selected2.csv", row.names = F)

kk3 = which(res$estimate<0 & -log10(res$p.value) > -log10(pval.threshold))
write.csv(INTtable[kk3, ], "data/intensity_table_selected3.csv", row.names = F)

# ----- ranking ----
res2 = res %>%
  rename(pepSeq = gene)

master = left_join(INTtable, res2)

# select peptides that are not 0 in any condition but b5
# remove those that are 0 at no b5
names(master)
j1 = apply(master[, c(6:11)], 1, function(x){
  any(x == 0)
}) %>%
  which()

master = master[-j1, ]

# select downregulated peptides
master = master[master$estimate < 0, ]

plot(density(master$qvalues))
which(master$qvalues < 0.2)

plot(density(master$p.value))

jj = which(master$p.value < 0.05)
master[jj,] %>% View()

write.csv(master[jj, ], "data/intensity_table_selectedMaster.csv", row.names = F)


