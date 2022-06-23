
# ----- ROC curve ------

getROCcurve = function(ttrue_mean, tsims_mean, substrate, threshPerc = 0.01, retValsOnly=F) {
  
  PRECISION = function(TP = "", FP = "") {
    return(as.numeric(TP) / (as.numeric(TP) + as.numeric(FP)))
  }
  
  RECALL = function(TP = "", P = "") {
    return(as.numeric(TP) / as.numeric(P))
  }
  
  SENSITIVITY = function(TP = "", P = "") {
    return(as.numeric(TP) / as.numeric(P))
  }
  
  SPECIFICITY = function(TN = "", N = "") {
    return(as.numeric(TN) / as.numeric(N))
  }
  
  k = order(ttrue_mean, decreasing = F)
  ttrue_mean = ttrue_mean[k]
  
  predThresh = log(threshPerc+pseudo)
  trueBinary = ifelse(ttrue_mean > predThresh, 1, 0)
  
  tsims_mean = tsims_mean[k]
  
  th_range = c(-Inf, seq(min(tsims_mean), max(tsims_mean), length.out = 300), Inf)
  
  sens = rep(NA, length(th_range))
  spec = rep(NA, length(th_range))
  
  prec = rep(NA, length(th_range))
  rec = rep(NA, length(th_range))
  
  for (tr in seq_along(th_range)) {
    
    cntDat = trueBinary
    cntPred = ifelse(tsims_mean > th_range[tr], 1, 0)
    
    # P, N, TP, TN, FP
    P = which(cntDat == 1) %>% length()
    N = which(cntDat == 0) %>% length()
    
    TP = which(cntPred == 1 & cntDat == 1) %>% length()
    TN = which(cntPred == 0 & cntDat == 0) %>% length()
    
    FP = which(cntPred == 1 & cntDat == 0) %>% length()
    FN = which(cntPred == 0 & cntDat == 1) %>% length()
    
    sens[tr] = SENSITIVITY(TP, P)
    spec[tr] = SPECIFICITY(TN, N)
    
    prec[tr] = PRECISION(TP, FP)
    rec[tr] = RECALL(TP, P)
  }
  
  curve = data.frame(score = th_range,
                     precision = prec,
                     recall = rec,
                     sensitivity = sens,
                     specificity = spec)
  
  # AUC
  pr.na = which(! is.na(curve$precision | curve$recall))
  pr.auc = AUC(curve$recall[pr.na],
               curve$precision[pr.na])
  
  
  roc.na = which(! is.na(curve$sensitivity | curve$specificity))
  roc.auc = AUC(curve$specificity[roc.na],
                curve$sensitivity[roc.na])
  
  if (retValsOnly) {
    
    curve$pr_auc = pr.auc
    curve$roc_auc = roc.auc
    
    return(curve)
    
  } else {
    
    roc.curve = curve %>%
      ggplot() +
      geom_path(aes(1 - specificity, sensitivity)) + 
      geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
      xlim(c(0,1)) +
      ylim(c(0,1)) +
      ggtitle(substrate,
              subtitle = paste0("ROC curve, AUC: ", roc.auc %>% round(4)))
    
    pr.curve = curve %>%
      ggplot() +
      geom_path(aes(recall, precision)) +
      geom_abline(intercept = length(which(trueBinary == 1))/length(trueBinary),
                  slope = 0, linetype = "dotted") +
      xlim(c(0,1)) +
      ylim(c(0,1)) +
      ggtitle(substrate,
              subtitle = paste0("PC curve, AUC: ", pr.auc %>% round(4)))
    
    
    return(list(roc.curve, pr.curve))
    
  }
  
}


