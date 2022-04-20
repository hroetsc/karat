### PCPS brainstorming ###
# description:  positions of fwd/rev cis spliced peptides of length N generated from substrate of length L
# input:        substrate length L, product length N
# output:       positions
# author:       HPR

library(dplyr)


# ----- forward cis -----
# L: substrate length
# N: peptide length
# Lext: minimal splice-reactant length
getCis = function(L, N, Lext=1) {
  
  if (L-1 >= N) {
    
    allCis = sapply(seq(1,L-N), function(i){
      
      sapply(seq(i+Lext-1,N-Lext+i-1), function(j){
        
        sr2 = N-j+i-1
        sapply(seq(j+2, L-sr2+1), function(k){
          
          n = k+sr2-1
          return(paste(i,j,k,n,sep = "_"))
        })
      })
    })
    
    CIS = unlist(allCis)
    
  } else {
    CIS = NA
  }
  
  return(CIS)
}


numCis = function(L, N, Lext=1){
  return((L-N)*(N-2*Lext+1)*(L-N+1)*0.5)
}




# ----- reverse cis -----
getRevCis = function(L, N, Lext=1) {
  
  if (L >= N) {
    
    allRevCis = sapply(seq(1,L-N+1), function(k){
      
      sapply(seq(k+Lext-1,N-Lext+k-1), function(n){
        
        sr1 = N-n+k-1
        sapply(seq(n+1, L-sr1+1), function(i){
          
          j = i+sr1-1
          return(paste(i,j,k,n,sep = "_"))
        })
      })
    })
    
    REVCIS = unlist(allRevCis)
    
    
  } else {
    REVCIS = NA
  }
  
  
  return(REVCIS)
}

numRevCis = function(L, N, Lext=1){
  return((L-N+1)*(N-2*Lext+1)*(L-N+2)*0.5)
}



# ----- trans -----

getTrans = function(L, N, Lext=1) {
  
  allTrans = sapply(seq(1,L-Lext), function(i){
    
    sapply(seq(i+Lext-1,N-Lext+i-1), function(j){
      
      sr2 = N-j+i-1
      kmin = if(i-sr2+1 > 0) i-sr2+1 else 1
      sapply(seq(kmin, j), function(k){
        
        n = k+sr2-1
        return(paste(i,j,k,n,sep = "_"))
      })
    })
  })
  
  TRANS = unlist(allTrans)
  
  return(TRANS)
}

