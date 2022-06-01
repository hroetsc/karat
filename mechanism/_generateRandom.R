### protein PCPS ###
# description:  generate random databases for appropriate comparison of
#               identified peptides from proteins
# input:        identified sequences
# output:       random DBs, incl quality control
# author:       HPR

library(dplyr)
library(stringr)
library(seqinr)
library(Biostrings)

source("src/invitroSPI_utils.R")
source("src/_generateRandom_qc.R")

# factor of which the random DB should be larger than the actual DB
rndSize = 10
# determined by Mascot server settings/ Prosit
minLength = 5
maxLength = 40

#### INPUT ###
load("data/invitroSPI.RData")
load("data/aSPIre.RData")

### MAIN PART ###
# DB = Kinetics %>% uniquePeptides()
DB = ProteasomeDB
S = DB$substrateID %>% unique()
# S = c("TSN2","TSN89","TSN5")
tps = "all"


# ----- 2) function to retrieve all spliced peptides for a given length -----
genPCP = function(subSeq = "") {
  print("all possible PCPs")
  
  subEnc = strsplit(subSeq, "") %>%
    unlist()
  
  t = length(subEnc)^2  # theoretical amount of PCPs
  pcp.products = rep(NA, t)
  pos = rep(NA, t)
  
  counter = 1
  
  for (o in 1:length(subEnc)) {
    
    for (p in 1:length(subEnc)) {
      
      if (! p < o) {
        cnt.pcp = paste(subEnc[o:p], collapse = "")
        pcp.products[counter] = cnt.pcp
        pos[counter] = paste(o, p, sep = "_")
        
        counter = counter + 1
        
      }
    }
  }
  
  res = data.frame(product = pcp.products,
                   positions = pos,
                   type = "PCP") %>% na.omit()
  
  res$product = as.character(res$product)
  
  return(res)
  
}



# ----- 3) generate the random DB ------
randomDB.allSubs = list()
for (s in 1:length(S)) {
  
  print("----------------------------------")
  print("generate random DB for substrate:")
  print(S[s])
  print("----------------------------------")
  
  
  randomDB.allTPs = list()
  for (t in 1:length(tps)) {
    
    print("----------------------------------")
    print("time point:")
    print(tps[t])
    print("----------------------------------")
    
    if (tps == "all") {
      cntDB = DB[DB$substrateID == S[s], ]
    } else {
      cntDB = DB[DB$substrateID == S[s] & DB$digestTime == tps[t], ]
    }
    
    
    if (nrow(cntDB) > 0) {
      
      subSeq = cntDB$substrateSeq[1]
      pcpALL = genPCP(subSeq)
      
      pcpALL$pos1 = str_split_fixed(pcpALL$positions, coll("_"), Inf)[,1] %>% as.numeric()
      pcpALL$pos2 = str_split_fixed(pcpALL$positions, coll("_"), Inf)[,2] %>% as.numeric()
      
      # maximum length cutoff
      if (any(nchar(pcpALL$product) > maxLength)) {
        pcpALL = pcpALL[-which(nchar(pcpALL$product) > maxLength), ]
      }
      
      # generate all combinations of PCPs and cut too long sequences
      # to improve sampling success
      reallyall = crossing(pcpALL %>% select(product,pos1,pos2),
                           pcpALL %>% select(product,pos1,pos2),
                           .name_repair = "unique")
      colnames(reallyall) = c("SR1","pos1","pos2","SR2","pos3","pos4")
      reallyall$product = do.call(paste, c(reallyall[c("SR1","SR2")], sep=""))
      reallyall$n = nchar(reallyall$product)
      if (max(reallyall$n) > maxLength) {
        reallyall = reallyall[reallyall$n <= maxLength, ]
      }
      reallyall = reallyall[reallyall$n >= minLength, ]
      
      # get number of products
      Ns = c(
        "cis" = length(which(cntDB$spliceType == "cis")) * rndSize,
        "revCis" = length(which(cntDB$spliceType == "revCis")) * rndSize,
        "trans" = length(which(cntDB$spliceType == "trans")) * rndSize,
        "PCP" = length(which(cntDB$spliceType == "PCP")) * rndSize
      )
      
      counter = c( "cis" = 0, "revCis" = 0, "trans" = 0,
                   "PCP" = length(which(cntDB$spliceType == "PCP")) * rndSize)
      
      pspAll = matrix(NA, ncol = 3, nrow = 5*(Ns["cis"]+Ns["revCis"]+Ns["trans"])) %>% as.data.frame()
      colnames(pspAll) = c("product","positions","type")
      
      if (nrow(reallyall) <= Ns["cis"]+Ns["revCis"]+Ns["trans"]){
        randomVec = sample(nrow(reallyall), nrow(reallyall), replace = F)
      } else {
        randomVec = sample(nrow(reallyall), Ns["cis"]+Ns["revCis"]+Ns["trans"], replace = F)
      }
      
      i = 1
      while (counter != Ns & i<=length(randomVec)) {
        
        # sample two PCPs and create PSP
        x = reallyall[randomVec[i], ]
        
        # get product type
        intv = x$pos3-x$pos2
        
        if (counter["cis"] < Ns["cis"] & intv > 0 & x$pos3 > x$pos2) {
          
          pspAll$product[i] = x$product
          pspAll$positions[i] = paste(x$pos1, x$pos2, x$pos3, x$pos4, sep = "_")
          pspAll$type[i] = "cis"
          counter["cis"] = counter["cis"]+1
          
        } else if (counter["revCis"] < Ns["revCis"] & intv <= 0 & x$pos4 <= x$pos1) {
          
          pspAll$product[i] = x$product
          pspAll$positions[i] = paste(x$pos1, x$pos2, x$pos3, x$pos4, sep = "_")
          pspAll$type[i] = "revCis"
          counter["revCis"] = counter["revCis"]+1
          
          
        } else if (counter["trans"] < Ns["trans"] & intv <= 0 & x$pos4 > x$pos1) {
          
          pspAll$product[i] = x$product
          pspAll$positions[i] = paste(x$pos1, x$pos2, x$pos3, x$pos4, sep = "_")
          pspAll$type[i] = "trans"
          counter["trans"] = counter["trans"]+1
          
        }
        
        # update counter
        i = i+1
        
        # print the progress
        if (i %% 1000 == 0) {
          print(i)
          print(counter)
        }
        
      }
      
      # minimum length cutoff
      pcpALL = pcpALL[-which(nchar(pcpALL$product) < minLength), ]
      # sample PCPs
      if (nrow(pcpALL) > Ns["PCP"]) {
        pcpALL = pcpALL[sample(nrow(pcpALL), Ns["PCP"], replace = F), ]
      }
      
      ALL = rbind(pcpALL %>% select(product,positions,type),
                  pspAll %>% select(product,positions,type)) %>%
        as.data.frame() %>%
        na.omit
      names(ALL) = c("pepSeq","positions","spliceType")
      
      ALL$digestTime = tps[t]
      ALL$substrateID = S[s]
      ALL$substrateSeq = subSeq
      ALL$productType = "PSP"
      ALL$productType[ALL$spliceType == "PCP"] = "PCP"
      
      randomDB.allTPs[[t]] = ALL
      
    }
  }
  
  randomDB.allSubs[[s]] = plyr::ldply(randomDB.allTPs)
}

d = plyr::ldply(randomDB.allSubs)
# get multi-mappers
rndDB = mapping(d) %>%
  uniquePeptides()

# make sure no PCPs are assigned as PSP
k = which(str_detect(rndDB$positions, "^[:digit:]+_[:digit:]+$"))
rndDB$spliceType[k] = "PCP"
rndDB$productType[k] = "PCP"

# ----- 4) quality control (heatmaps) -----
# colour should be <= 2e-04

# pdf("results/RANDOM_quant.pdf", height = 4, width = 8)
pdf("results/RANDOM.pdf", height = 4, width = 8)
plotHeatmap(rndDB)

cnt = disentangleMultimappers.Type(rndDB)
print((table(cnt$spliceType) / sum(table(cnt$spliceType)) *100) %>% round(2))
print((table(DB$spliceType) / sum(table(DB$spliceType)) *100) %>% round(2))
dev.off()


### OUTPUT ###
# save(rndDB, file = "data/randomDB_aSPIre.RData")
save(rndDB, file = "data/randomDB.RData")

