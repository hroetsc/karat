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
DB = Kinetics %>% uniquePeptides()
# DB = ProteasomeDB
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
      L = nchar(subSeq)
      pcpALL = genPCP(subSeq)
      
      pcpALL$pos1 = str_split_fixed(pcpALL$positions, coll("_"), Inf)[,1] %>% as.numeric()
      pcpALL$pos2 = str_split_fixed(pcpALL$positions, coll("_"), Inf)[,2] %>% as.numeric()
      
      # maximum length cutoff
      if (any(nchar(pcpALL$product) > maxLength)) {
        pcpALL = pcpALL[-which(nchar(pcpALL$product) > maxLength), ]
      }
      
      # get number of products
      Ns = c(
        "cis" = length(which(cntDB$spliceType == "cis")) * rndSize,
        "revCis" = length(which(cntDB$spliceType == "revCis")) * rndSize,
        "trans" = length(which(cntDB$spliceType == "trans")) * rndSize,
        "PCP" = length(which(cntDB$spliceType == "PCP")) * rndSize
      )
      
      counter = c( "cis" = 0, "revCis" = 0, "trans" = 0,
                   "PCP" = length(which(cntDB$spliceType == "PCP")) * rndSize)
      
      pspAll = matrix(NA, ncol = 3, nrow = 2*(Ns["cis"]+Ns["revCis"]+Ns["trans"])) %>% as.data.frame()
      colnames(pspAll) = c("product","positions","type")
      
      
      i = 1
      while (!all(counter == Ns)) {
        
        # sample four positions and check the product type
        pos = sample(1:L, size = 4, replace = T)
        N = pos[2]-pos[1]+pos[4]-pos[3]+2
        
        if (N >= minLength & N <= maxLength & pos[1]<=pos[2] & pos[3]<=pos[4]) {
          
          # get product type
          intv = pos[3]-pos[2]
          
          if (counter["cis"] < Ns["cis"] & intv > 0 & pos[3] > pos[2]) {
            
            pspAll$product[i] = paste(substr(subSeq, pos[1], pos[2]), substr(subSeq, pos[3], pos[4]), sep = "")
            pspAll$positions[i] = paste(pos[1], pos[2], pos[3], pos[4], sep = "_")
            pspAll$type[i] = "cis"
            counter["cis"] = counter["cis"]+1
            
            i = i+1
            
          } else if (counter["revCis"] < Ns["revCis"] & intv <= 0 & pos[4] < pos[1]) {
            
            pspAll$product[i] = paste(substr(subSeq, pos[1], pos[2]), substr(subSeq, pos[3], pos[4]), sep = "")
            pspAll$positions[i] = paste(pos[1], pos[2], pos[3], pos[4], sep = "_")
            pspAll$type[i] = "revCis"
            counter["revCis"] = counter["revCis"]+1
            
            i = i+1
            
          } else if (counter["trans"] < Ns["trans"] & intv <= 0 & pos[4] >= pos[1]) {
            
            pspAll$product[i] = paste(substr(subSeq, pos[1], pos[2]), substr(subSeq, pos[3], pos[4]), sep = "")
            pspAll$positions[i] = paste(pos[1], pos[2], pos[3], pos[4], sep = "_")
            pspAll$type[i] = "trans"
            counter["trans"] = counter["trans"]+1
            
            i = i+1
          }
          
        }
        
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

pdf("results/RANDOM_quant.pdf", height = 4, width = 8)
# pdf("results/RANDOM.pdf", height = 4, width = 8)
plotHeatmap(rndDB)

cnt = disentangleMultimappers.Type(rndDB)
print((table(cnt$spliceType) / sum(table(cnt$spliceType)) *100) %>% round(2))
print((table(DB$spliceType) / sum(table(DB$spliceType)) *100) %>% round(2))
dev.off()


### OUTPUT ###
save(rndDB, file = "data/randomDB_aSPIre.RData")
# save(rndDB, file = "data/randomDB_smart.RData")

