### MARKOV CHAIN ###

# HILL 2002

# 1. create transition matrix

trans_mat <- function(from, to) {
  tab <- table(from, to)
  trans_mat <- tab/apply(tab,1,sum)
  list(trans = trans_mat, raw = tab)
}

# tr_ls <- trans_mat(states_trans2$From, states_trans2$To)
# tr <- tr_ls$trans
# 
# freq_mat <- tr_ls$raw



# 2. compute stationary state

steady_state <- function(pmat = NULL, qmat = NULL, plot = FALSE){
  
  if(!is.null(pmat)) {
    eig <- eigen(t(pmat))
    
    pi <- t(eig$vectors[,which(abs(eig$values-1) < 1e-12)])
  } 
  
  if(!is.null(qmat)) {
    eig <- eigen(t(qmat))
    
    pi <- t(eig$vectors[,which(abs(eig$values) < 1e-12)])
  }
  
  
  steady <- pi/sum(pi)
  colnames(steady) <- colnames(pmat)
  
  if(plot) barplot(steady)
  
  return(steady)
  
}

# steady <- steady_state(tr, plot = TRUE, freq = freq_mat)


# 3. Successional Transitions
# probability of colonization after disturbances (pioneer)
# P(colonization from pioneer by species j) = psj

# P(replacement of species j) = 1 - pjj

diag2 <- function(x) {
  if(is.matrix(x)) { diag(x)}
  if(is.array(x)) { apply(x,3,diag) }
}

eigen_left <- function(x) {
  if(is.matrix(x)) { eigen(t(x)) }
  if(is.array(x)) { eigen(t(x[,,"estimate"])) }
}

succession <- function(pmat, jump, ci = T) {
  
  if(ci) {
    states <- colnames(pmat$estimates)
    pmat <- pmat[1:4,1:4]
    
    dimnames(pmat)[[1]] <- states
    dimnames(pmat)[[2]] <- states
    
    # colonization from pioneer to 
    p_col <- pmat[,"Pioneer",1:3]
    
    # replacement by a species
    pmatmp <- pmat
    diag(pmatmp[,,1]) <- 0
    diag(pmatmp[,,2]) <- 0
    diag(pmatmp[,,3]) <- 0
    p_replaceby <- colSums(pmatmp)/(length(states)-1)
    
  } else {
    states <- colnames(pmat)
    class(pmat) <- "matrix"
    
    # colonization from pioneer to 
    p_col <- pmat[,"Pioneer"]
    
    # replacement by a species
    pmatmp <- pmat
    diag(pmatmp) <- 0
    p_replaceby <- colSums(pmatmp)/(length(states)-1)
  }
  
  
  # replacement of a species
  p_replaceof <- 1 - diag2(pmat)
  
  # peristence
  p_persist <- diag2(pmat)
  
  # turnover
  turn_time <- 1/p_replaceof
  
  # convergence
  eig <- eigen_left(pmat)
  convergence <- eig$values[1]/eig$values[2]
  
  # Damping ratio (*10 because pmatrix on 10 years) = half life to steady state
  damping = log(2)/log(convergence)*10
  
  # Recurrence
  pi <- t(eig$vectors[,1])
  
  steady <- pi/sum(pi)
  colnames(steady) <- states
  
  if(is.matrix(pmat)) {
    recur <- (1 - steady)/(steady*(p_replaceof))
  } else {
    recur <- (1 - steady)/(steady*(p_replaceof[,"estimate"]))
  }
  
  # Entropy
  jump <- jump[1:4,]
  
  log_jump <- log(jump)
  log_jump[which(log_jump=="-Inf")] <- 0
  
  entropy_st <- -apply(jump*log_jump, 3, rowSums)
  
  entropy_all <- sum(steady*entropy_st[,"estimate"])
  
  entropy_rel <- entropy_all/log(length(states)-1)
  
  # result
  list(p_colonization = p_col, 
       p_replaceof = p_replaceof, 
       p_replaceby = p_replaceby, 
       p_persist = p_persist,
       turn_time = turn_time, 
       steady = steady,
       convergence = convergence, 
       damping = damping,
       recurrence_time = recur,
       entropy_st = entropy_st,
       entropy_all = entropy_all, entropy_rel = entropy_rel)
}
# succession(tr, disturb = "Pioneer")


##log linear model

# x=table(From=states_trans2$From, To=states_trans2$To, 
#         nat_disturb=states_trans2$nat_disturb, 
#         ecoreg=states_trans2$Subzone)
# sat.model = loglm(~ From * To * nat_disturb * ecoreg, data=x)
# 
# ti = as.data.frame(x)
# glm.model = glm(Freq ~ From * To * nat_disturb * ecoreg, data=ti, family=poisson)
# anova(glm.model, test="Chisq")
# 
# 
# anova(update(glm.model,.~.-(To:nat_disturb:ecoreg)),test="Chisq")
# 
# step.model=step(glm.model, direction="backward") 

