### MAKE PREDICTION FROM OPTIMIZED PARAMETERS #####

predict_msm <- function(par_opt, newdata, newdata_p, t = 1) {
  
  form_all <- formula(paste("~ ", paste(colnames(newdata), collapse= "+")))
  form_p <- formula(paste("~ ", paste(colnames(newdata_p), collapse= "+")))
  
  mm <- model.matrix(form_all, newdata)
  mm_p <-  model.matrix(form_p, newdata_p)
  names_params <- names(par_opt)
  n_obs <- ifelse(is.null(nrow(mm)), 1, nrow(mm))
  
  # Model 
  BM <- exp(mm %*% par_opt[startsWith(names_params, "bm")])
  BP <- exp(mm_p %*% par_opt[startsWith(names_params, "bp")])
  BB <- -base::rowSums(cbind(BP, BM))
  
  MB <- exp(mm %*% par_opt[startsWith(names_params, "mb")])
  MP <- exp(mm_p %*% par_opt[startsWith(names_params, "mp")])
  MT <- exp(mm %*% par_opt[startsWith(names_params, "mt")])
  MM <- -base::rowSums(cbind(MB, MP, MT))
  
  PB <- exp(mm %*% par_opt[startsWith(names_params, "pb")])
  PM <- exp(mm %*% par_opt[startsWith(names_params, "pm")])
  PT <- exp(mm %*% par_opt[startsWith(names_params, "pt")])
  PP <- -base::rowSums(cbind(PB, PM, PT))
  
  TM <- exp(mm %*% par_opt[startsWith(names_params, "tm")])
  TP <- exp(mm_p %*% par_opt[startsWith(names_params, "tp")])
  TT <- -base::rowSums(cbind(TM, TP))
  
  BT <- TB <- rep(0, n_obs) # impossible transition
  
  # Q matrix
  Q_df <- as.data.frame(t(cbind.data.frame(BB, BM, BP, BT, 
                                           MB, MM, MP, MT, 
                                           PB, PM, PP, PT, 
                                           TB, TM, TP, TT)))
  
  P <- list()
  # Compute P matrix using matrix exponential 
  for(i in t) { 
    Q_tmp <- rbind.data.frame(Q_df, t=i)
    P_tmp <- lapply(Q_tmp, FUN = pmat)
    P_tmp <- lapply(P_tmp, "dimnames<-", list(c("B", "M", "P", "T"), c("B", "M", "P", "T")))
    P[[i]] <- P_tmp
  }
  
  
  
  return(P)
}

