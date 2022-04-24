### PCPS brainstorming ###
# description:  positions of fwd/rev cis spliced peptides of length N generated from substrate of length L
#               number of spliced/non-spliced peptides of length N generated from substrate L
# input:        substrate length L, product length N
# output:       positions, number of peptides, number of splice sites
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

# PSP containing either N or C term of substrate
numCis_1term = function(L, N, Lext=1){
  if (L-1 >= N) {
    return((N-2*Lext+1)*(L-N))
  } else {
    return(0)
  }
}

# PSP containing both N and C term of substrate
numCis_bothterm = function(L, N, Lext=1){
  if (L-1 >= N) {
    return((N-2*Lext+1))
  } else {
    return(0)
  }
}

# number of sequences with restricted intervening sequence length
# Imax: maximum intervening sequence length
numCis_maxInterv = function(L, N, Imax, Lext=1){
  if (L-1 >= N) {
    # return((L+3-N-Imax)*(N-2*Lext+1)*(Imax-2) + 0.5*(N-2*Lext+1)*(Imax-3)*(Imax-2))
    return(-0.5*(-2+Imax)*(1-2*Lext+N)*(-3+Imax-2*L+2*N))
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
# peptide positions
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

# PSP containing either N or C term of substrate
numRevCis_1term = function(L, N, Lext=1){
  if (L >= N) {
    return((N-2*Lext+1)*(L-N+1))
  } else {
    return(0)
  }
}

# PSP containing both N and C term of substrate
numRevCis_bothterm = function(L, N, Lext=1){
  if (L >= N) {
    return((N-2*Lext+1))
  } else {
    return(0)
  }
}

# number of sequences with restricted intervening sequence length
# Imax: maximum intervening sequence length
numRevCis_maxInterv = function(L, N, Imax, Lext=1){
  if (L-1 >= N & Imax >= N) {
    # return((N-2*Lext+1)*((L-Imax+1)*(Imax-N+1) + 0.5*(Imax-N)*(Imax-N+1)))
    return((1+Imax-N)*(1-0.5*Imax+L-0.5*N)*(1-2*Lext+N))
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
# peptide positions
getTrans = function(L, N, Lext=1) {
  
  if (L+Lext >= N) {
    
    allTrans = sapply(seq(Lext,N-Lext), function(sr1){
      
      # SR1
      Is = seq(1, L-sr1+1)
      Js = Is+sr1-1
      
      # SR2
      Ks = seq(1,L-N+sr1+1)
      Ns = Ks+N-sr1-1
      
      SR1 = paste(Is,Js, sep = "_")
      SR2 = paste(Ks,Ns, sep = "_")
      # all combinations
      PEP = outer(SR1, SR2, paste, sep="_") %>% as.vector()
      
      # get overlap
      idx = sapply(PEP, function(x){
        y = strsplit(x, "_") %>% unlist() %>% as.numeric()
        z = (max(y)-min(y)+1 <= N-Lext) & (max(y)-min(y)+1 > N/2) | any(y[1]:y[2] %in% y[3]:y[4])
        return(z)
      }) %>% which()

      return(PEP[idx])
      
    })
    
    TRANS = unlist(allTrans) %>% as.vector()
    
  } else {
    TRANS = NA
  }
  
  return(TRANS)
}


# number of peptides
numTrans = function(L, N, Lext=1){
  
  if (L+Lext >= N) {
    # # multiply (L-sr1+1)*(L-N+sr1+1) over [Lext, N-Lext] out and separate terms
    # # not dependent on sr1
    # not_dependent = (N-2*Lext+1)*(L^2+2*L-L*N-N+1)
    # # linear dependency on sr1
    # linear_dependent = 0.5*N^2*(N-2*Lext+1)
    # # square dependency on sr1
    # square_part = (1/6)*Lext*(Lext+1)*(2*Lext+1) - (1/6)*(N-Lext-1)*((N-Lext-1)+1)*(2*(N-Lext-1)+1)
    # print(square_part)
    # all_sliding = not_dependent + linear_dependent - square_part
    # print(all_sliding)
    
    # all_sliding = -(1/6)*Lext*(1+Lext)*(1+2*Lext) - (1/6)*(1+2*Lext-2*N)*(Lext-N)*(1+Lext-N) + (1+L)*(1+L-N)*(1-2*Lext+N) + 0.5*N^2*(1-2*Lext+N)
    # print(all_sliding)
    
    all_sliding = (1/6)*(-1+2*Lext-N)*(-6-12*L-6*L^2-2*Lext+2*Lext^2+7*N+6*L*N-2*Lext*N-N^2)
    # print(all_sliding)
    # get number of cis PSPs and PCPs
    spliced = numCis(L,N) + numRevCis(L,N)
    # print(numCis(L,N))
    # print(numRevCis(L,N))
    # print(spliced)
    nonspliced = if (numPCP(L,N) > 0) numPCP(L,N) else 1
    # print((N-2*Lext+1)*nonspliced)
    
    # subtract from all sliding windows
    x = all_sliding - spliced - (N-2*Lext+1)*nonspliced
    return(x)
     
  } else {
    return(0)
  }
}

L=50
N=20
trans = getTrans(L,N)
length(trans)
numTrans(L,N)


# ----- trans-spliced from two different proteins -----

getTransProt = function(L1, L2, N, Lext=1) {
  
  if (min(L1,L2)+Lext >= N) {
    
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
    
  } else {
    
    TRANSPROT = NA
    
  }
  
  return(TRANSPROT)
}


# L1: length of protein 1
# L2: length of protein 2
numTransProt = function(L1, L2, N, Lext=1) {
  if (min(L1,L2)+Lext >= N) {
    
    # # multiply (L1-sr1+1)*(L2-N+sr1+1) over [Lext, N-Lext] out and separate terms
    # # not dependent on sr1
    # not_dependent = (N-2*Lext+1)*(L1*L2+L2-L1*N+L1-N+1)
    # # linear dependency on sr1
    # linear_dependent = 0.5*N*(-L2+L1+N)*(N-2*Lext+1)
    # # square dependency on sr1
    # square_part = (1/6)*Lext*(Lext+1)*(2*Lext+1) - (1/6)*(N-Lext-1)*((N-Lext-1)+1)*(2*(N-Lext-1)+1)
    # # not_dependent + linear_dependent - square_part
    
    x = (-1/6)*Lext*(1+Lext)*(1+2*Lext) - (1/6)*(1+2*Lext-2*N)*(Lext-N)*(1+Lext-N) + (1+L1)*(1+L2-N)*(1-2*Lext+N) + 0.5*N*(L1-L2+N)*(1-2*Lext+N)
    
    return(x)
    
  } else {
    return(0)
  }
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
    return(L-N)
  } else {
    return(0)
  }
}





