### Function to compute the lolikelihood of my multi state model ####
multistate_LL = function(params, dat, mm, mm_p = mm, names_params) {
  
  # params = intial parameters
  # dat = df containing state from (From), state to (To), time interval (time_interv), covariates
  # mm = model matrix with covariates for all the state transition
  # mm_p = model matrix with covariates for transition to pioneer
  # names_params = parameter names 
  
  st0 <- dat$From
  st1 <- dat$To
  time_interv <- dat$time_interv
  
  n_obs <- length(st0)
  
  # State num
  st0_num <- as.factor(st0)
  levels(st0_num) <- c("1","2","3","4")
  st0_num <- as.numeric(st0_num)
  
  st1_num <- as.factor(st1)
  levels(st1_num) <- c("1","2","3","4")
  st1_num <- as.numeric(st1_num)

  
  # Q matrix:
  
  # Model
  BM <- exp(mm %*% params[startsWith(names_params, "bm")])
  BP <- exp(mm_p %*% params[startsWith(names_params, "bp")])
  BB <- -base::rowSums(cbind(BP,BM))
  
  MB <- exp(mm %*% params[startsWith(names_params, "mb")])
  MP <- exp(mm_p %*% params[startsWith(names_params, "mp")])
  MT <- exp(mm %*% params[startsWith(names_params, "mt")])
  MM <- -base::rowSums(cbind(MB, MP, MT))
  
  PB <- exp(mm %*% params[startsWith(names_params, "pb")])
  PM <- exp(mm %*% params[startsWith(names_params, "pm")])
  PT <- exp(mm %*% params[startsWith(names_params, "pt")])
  PP <- -base::rowSums(cbind(PB, PM, PT))
  
  TM <- exp(mm %*% params[startsWith(names_params, "tm")])
  TP <- exp(mm_p %*% params[startsWith(names_params, "tp")])
  TT <- -base::rowSums(cbind(TM, TP))
  
  BT <- TB <- rep(0, n_obs) # impossible transition
  
  # Combine all qrs in a Q matrix and add time
  
  Q_df <- as.data.frame(do.call(rbind, list(BB, BM[,1], BP[,1], BT, 
                              MB[,1], MM, MP[,1], MT[,1], 
                              PB[,1], PM[,1], PP, PT[,1], 
                              TB, TM[,1], TP[,1], TT, 
                              time_interv)))
  
  
  # Compute P matrix using matrix exponential 
    # initiate the cluster
  P <- parLapply(cl, X= Q_df, fun = pmat)
   
  
  #P <- lapply(X= Q_df, FUN = pmat)
  
  loglik <- sum(log(mapply(get_q, P=P, st0=st0_num, st1=st1_num)))
  
  if(is.infinite(loglik)) loglik = -.Machine$double.xmax
  if(is.nan(loglik)) loglik = -.Machine$double.xmax
  if(is.na(loglik)) print("loglik: na values!")
  
  return(-loglik)
}

# Function to compute matrix exponential on lines of a df ####
pmat <- function(qt) {
  qmat <- matrix(qt[-17], nrow = 4, ncol = 4, byrow = TRUE)
  P <- expm(qt[17] * qmat, method = "Higham08.b")
}

get_q <- function(P, st0, st1) { P[st0, st1] }

# test
# np <- detectCores(logical = FALSE)
# cl <- makeForkCluster(np)
# system.time({multistate_LL(params = par_init, dat = states_trans, mm = mm_list[[1]], mm_p = mm_list[[2]], names_params = names(par_init))})
#stopCluster(cl)


