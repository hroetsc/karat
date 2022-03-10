

quantileNorm = function(M) {
  
  # order matrix column-wise
  ranks = apply(M,2,base::order) %>% as.matrix()
  Mo = apply(M,2,function(x){x[base::order(x)]})
  
  # get row means
  mu = rowMeans(Mo)
  newRanks = mu[order(mu)]
  
  # replace ordered matrix with row means
  Mo_rpl = apply(Mo,2,function(x){mu})
  
  # bring in old order
  Mrec = apply(ranks,2,function(x){
    Mo_rpl[x]
  })
  
  return(Mrec)
}


medianPolish = function(M) {
  
  M2 = as.matrix(M)
  
  # get row medians
  rMd = apply(M2,1,median)
  # overall effect is median of row means
  overall = median(rMd)
  # subtract row meadians
  M2 = M2 - rMd
  # subtract overall effect from each row median
  rMd = rMd-overall
  
  # do the same for column
  cMd = apply(M,2,median)
  overall = overall + median(cMd)
  M2 = M2-as.matrix(cMd)[rep(1,nrow(M2)), ]
  cMd = cMd-overall
  
  M3 = medpolish(M2)$residuals
  
}
