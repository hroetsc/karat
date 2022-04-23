### PCPS brainstorming ###
# description:  positions of fwd/rev cis spliced peptides of length N generated from substrate of length L
#               number of spliced/non-spliced peptides of length N generated from substrate L
# input:        substrate length L, product length N
# output:       positions, number of peptides
# author:       HPR

library(dplyr)


# ----- forward cis -----
# L: substrate length
# N: peptide length
# Lext: minimal splice-reactant length

# peptide positions
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

# number of peptides
numCis = function(L, N, Lext=1){
  if (L-1 >= N) {
    return((L-N)*(N-2*Lext+1)*(L-N+1)*0.5)
  } else {
    return(0)
  }
}

numCis_1term = function(L, N, Lext=1){
  if (L-1 >= N) {
    return((N-2*Lext+1)*(L-N))
  } else {
    return(0)
  }
}

numCis_bothterm = function(L, N, Lext=1){
  if (L-1 >= N) {
    return((N-2*Lext+1))
  } else {
    return(0)
  }
}


# number of splice sites
numCis_sites = function(L,N,Lext=1) {
  if (L-1 >= N) {
    return((N-2*Lext+1)*(L-N) + 0.5*(L-N-1)*(L-N))
  } else {
    return(0)
  }
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

# number of peptides
numRevCis = function(L, N, Lext=1){
  if (L >= N) {
    return((L-N+1)*(N-2*Lext+1)*(L-N+2)*0.5)
  } else {
    return(0)
  }
}

numRevCis_1term = function(L, N, Lext=1){
  if (L >= N) {
    return((N-2*Lext+1)*(L-N+1))
  } else {
    return(0)
  }
}

numRevCis_bothterm = function(L, N, Lext=1){
  if (L >= N) {
    return((N-2*Lext+1))
  } else {
    return(0)
  }
}

# number of splice sites
numRevCis_sites = function(L,N,Lext=1) {
  if (L-1 >= N) {
    return(0.5*(L-N+1)*(L-N+2))
  } else {
    return(0)
  }
}


# ----- trans -----

getTrans = function(L, N, Lext=1) {
  
  allTrans = sapply(seq(1,L-Lext), function(i){
    
    sapply(seq(i+Lext-1,N-Lext+i-1), function(j){
      
      sr2 = N-j+i-1
      kmin = if(i-sr2+1 > 0) i-sr2+1 else 1
      print(i-sr2+1)
      sapply(seq(kmin, j), function(k){
        
        n = k+sr2-1
        return(paste(i,j,k,n,sep = "_"))
      })
    })
  })
  
  TRANS = unlist(allTrans)
  
  return(TRANS)
}


sapply(1:length(i), function(x){
  any(c(i[x]:j[x]) %in% c(k[x]:n[x]))
}) %>%
  all()

k2 = which(pos[,2]==pos[,3])


numTrans = function(L, N, Lext=1){
  
  x = L^2+L*(1-2*N)+(1+Lext)*(2*N-Lext-2)
  part1 = (1/8)*x^2 + (1/4)*x
  part2 = 0.5*(N-1)*(-1*Lext+L-1)*(-1*Lext+L)
  
  return(part1+part2)
}


L=15
N=5
trans = getTrans(L,N)
length(trans)

numTrans(L,N)


# ----- trans-spliced from two different proteins -----

getTransProt = function(L1, L2, N, Lext=1) {
  
  transProt = sapply(seq(Lext, N-Lext), function(sr1){
    sr2 = N-sr1
    
    # SR1 from protein 1, SR2 from protein 2
    I1 = seq(1,L1-sr1+1)
    J1 = sapply(I1, function(i1){
      j1 = sr1+i1-1
    })
    I2 = seq(1,L2-sr2+1)
    J2 = sapply(I2, function(i2){
      j2 = sr2+i2-1
    })
    
    PROT1 = paste(I1,J1,sep = "_")
    PROT2 = paste(I2,J2,sep = "_")
    P1 = outer(PROT1,PROT2,paste,sep="+") %>% as.vector()
    
    return(P1)
    
  })
  
  TRANSPROT = unlist(transProt)
  
  return(TRANSPROT)
}


# L1: length of protein 1
# L2: length of protein 2
numTransProt = function(L1, L2, N, Lext=1) {
  
  # multiply (L1-sr1+1)*(L2-N+sr1+1) over [Lext, N-Lext] out and separate terms
  # not dependent on sr1
  not_dependent = (N-2*Lext+1)*(L1*L2+L2-L1*N+L1-N+1)
  # linear dependency on sr1
  linear_dependent = 0.5*N*(-L2+L1+N)*(N-2*Lext+1)
  # square dependency on sr1
  square_part = (1/6)*Lext*(Lext+1)*(2*Lext+1) - (1/6)*(N-Lext-1)*((N-Lext-1)+1)*(2*(N-Lext-1)+1)
  
  return(not_dependent + linear_dependent - square_part)
}



# ----- non-spliced -----

getPCP = function(L,N) {
  
  if (L > N) {
    
    allPCP = sapply(seq(1,L-N+1), function(i){
      j = N+i-1
      return(paste(i,j,sep = "_"))
    })
    
    PCP = unlist(allPCP)
    
  } else {
    PCP = NA
  }
  
  return(PCP)
}


numPCP = function(L,N) {
  if (L > N) {
    return(L-N+1)
  } else {
    return(0)
  }
}

numCleavageSites = function(L,N) {
  if (L > N) {
    return(L-N+1)
  } else {
    return(0)
  }
}





