### karat projetc - PCPS mechanism ###
# description:  compare specificity of SR1 with SR2
# input:        invitroSPI: Roetschke SciData, EGFR ery, WT substrates --> spliced peptides only
#               random database for qualitative DB 
# output:       sequence divergence of SR1 vs. SR2
# author:       HPR

library(dplyr)
library(stringr)
library(seqinr)
library(protr)
library(Peptides)
library(transport)
library(ggplot2)
library(uwot)
library(class)
library(dbscan)
library(text2vec)

theme_set(theme_classic())


### INPUT ###
load("results/SRspecificity/DBproc.RData")
load("results/SRspecificity/randomDBproc.RData")


### MAIN PART ###
# ----- fetch AA index database -----
data("aaindex")
aaindices = matrix(ncol = 20, nrow = length(aaindex))
rownames(aaindices) = names(aaindex)
colnames(aaindices) = aaindex[[1]][["I"]] %>% names() %>% a()

for (k in 1:length(aaindex)) {
  cnt = aaindex[[k]][["I"]]
  aaindices[k, ] = cnt
}

aaindices = aaindices[, order(colnames(aaindices))] %>%
  na.omit() %>%
  t() %>%
  as.data.frame()

# scale between 0 and 1?

aaindices$AA = rownames(aaindices)


# ----- aa propensities of splice-reactants -----

getAAprop = function(df) {
  
  SR1feat = list()
  SR2feat = list()
  
  pb = txtProgressBar(min = 0, max = nrow(df), style = 3)
  for (i in 1:nrow(df)) {
    setTxtProgressBar(pb, i)
    
    x1 = data.frame(AA = unlist(strsplit(df$sr1[i],"")) %>% as.character()) %>%
      left_join(aaindices) %>% 
      select(-AA) %>% suppressMessages()
    SR1feat[[i]] = colMeans(x1)
    
    x2 = data.frame(AA = unlist(strsplit(df$sr2[i],"")) %>% as.character()) %>%
      left_join(aaindices) %>% 
      select(-AA) %>% suppressMessages()
    SR2feat[[i]] = colMeans(x2)
    
  }
  
  SR1 = plyr::ldply(SR1feat)
  SR2 = plyr::ldply(SR2feat)
  
  return(list(SR1 = SR1, SR2 = SR2))
}

DB_AAprop = getAAprop(DBproc)

# ----- determine similar amino acid indices -----
# matrix multiplication
matmult = function(v1, v2){
  return(as.numeric(v1) %*% as.numeric(v2))
}
# dot product
dot_product = function(v1, v2){
  p = matmult(v1, v2)/(sqrt(matmult(v1, v1)) * sqrt(matmult(v2, v2)))
  return(p)
}


AAind = aaindices %>% select(-AA) %>% t() %>% as.data.frame()

# UMAP
embedding = uwot::umap(AAind,
                       n_neighbors = 10,
                       min_dist = .3,
                       metric = "euclidean",
                       scale = "z",
                       init = "agspectral",
                       n_components = 5,
                       ret_extra = T,
                       verbose = T)

# plotting
plot(x = embedding[,1], y = embedding[,2],
     pch = 16,
     xlab = "UMAP 1", ylab = "UMAP 2",
     main = "embedding of AA index database")


# HDBSCAN clustering
cl = hdbscan(t(AAind), minPts = 1)
cl
plot(cl)
plot(cl$hc, main="HDBSCAN* Hierarchy", labels = names(AAind))


# get pairwise similarities between amino acid indices
Sims = list()
for (a in 1:(ncol(aaindices)-1)) {
  Sims[[a]] = dot_product(DB_AAprop$SR1[,a], DB_AAprop$SR2[,a])
}

names(Sims) = colnames(aaindices)[1:(ncol(aaindices)-1)]
Sims = unlist(Sims)
summary(Sims)

# most dissimilar amino acid indices between SR1 and SR2
plot(density(Sims))
extracted = Sims[Sims < 0.5]

descriptions = sapply(names(extracted), function(x){
  aaindex[names(aaindex) %in% x][[1]][["D"]]
})

extractedInd = data.frame(AAindex = names(extracted),
                          desc = descriptions,
                          sr1sr2_sim = extracted)


cl = hdbscan(t(AAind[rownames(AAind) %in% extractedInd$AAindex, ]), minPts = 2)
plot(cl$hc, main="HDBSCAN* Hierarchy", labels = names(AAind))

cl = hdbscan(t(AAind[rownames(AAind) %in% extractedInd$AAindex[grep("hydrophobicity",extractedInd$desc)], ]), minPts = 2)
plot(cl$hc, main="HDBSCAN* Hierarchy", labels = names(AAind))


# ----- determine difference of SR1 and SR2 for each AA index -----

getSim = function(v) {
  
  Sims = list()
  
  pb = txtProgressBar(min = 0, max = nrow(v), style = 3)
  v = as.matrix(v)
  
  for (i in 1:nrow(v)) {
    setTxtProgressBar(pb, i)
    Sims[[i]] = sapply(seq(i,nrow(v)), function(j){
      dist(v[i,], v[j,])
    })
  }
  
  dist(v, upper = T)
  
  return(Sims)
}


SR1 = DB_AAprop$SR1
SR1 = SR1[,colnames(SR1) %in% extractedInd$AAindex]

SR2 = DB_AAprop$SR2
SR2 = SR2[,colnames(SR2) %in% extractedInd$AAindex]




### OUTPUT ###
