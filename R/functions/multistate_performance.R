### Performance Measures for Multi-Class Problems ####

# from https://www.datascienceblog.net/post/machine-learning/performance-measures-multi-class-problems/

### Prediction ###
# pmat_fun <- function(mod, t, covariates, from) {
#   pmat <- pmatrix.msm(x = mod, t = t, covariates = as.list(covariates[i,]))
#   st_from <- which(states == from[i])
#   pvec <- pmat[st_from,]
#   p_ls[[i]] <- pvec
# }
# apply(dat, 1, function(y) pmat_fun(mod, dat['t'], dat['cov_names'], dat['from']))


### A. Hard classifier method (use predicted states, not probability) ####

### 1. Accuracy ####

calculate.accuracy <- function(predictions, ref.labels) {
  return(length(which(predictions == ref.labels)) / length(ref.labels))
}

### 2. Micro and macro averages of the F1-score from confusion matrix ####


get.conf.stats <- function(cm) {
  out <- vector("list", length(cm))
  for (i in seq_along(cm)) {
    x <- cm[[i]]
    tp <- x$table[x$positive, x$positive] 
    fp <- sum(x$table[x$positive, colnames(x$table) != x$positive])
    fn <- sum(x$table[colnames(x$table) != x$positie, x$positive])
    # TNs are not well-defined for one-vs-all approach
    elem <- c(tp = tp, fp = fp, fn = fn)
    out[[i]] <- elem
  }
  df <- do.call(rbind, out)
  rownames(df) <- unlist(lapply(cm, function(x) x$positive))
  return(as.data.frame(df))
}
get.micro.f1 <- function(cm) {
  cm.summary <- get.conf.stats(cm)
  tp <- sum(cm.summary$tp)
  fn <- sum(cm.summary$fn)
  fp <- sum(cm.summary$fp)
  pr <- tp / (tp + fp)
  re <- tp / (tp + fn)
  f1 <- 2 * ((pr * re) / (pr + re))
  return(f1)
}


get.macro.f1 <- function(cm) {
  c <- cm[[1]]$byClass # a single matrix is sufficient
  re <- sum(c[, "Recall"]) / nrow(c)
  pr <- sum(c[, "Precision"]) / nrow(c)
  f1 <- 2 * ((re * pr) / (re + pr))
  return(f1)
}


### B. Soft classifier ####

### 3. AUC generalization from Hand and Till ####

compute.A.conditional <- function(pred.matrix, i, j, ref.outcome) {
  # computes A(i|j), the probability that a randomly 
  # chosen member of class j has a lower estimated probability (or score) 
  # of belonging to class i than a randomly chosen member of class i
  
  # select predictions of class members
  i.idx <- which(ref.outcome == i)
  j.idx <- which(ref.outcome == j)
  pred.i <- pred.matrix[i.idx, i] # p(G = i) assigned to class i observations
  pred.j <- pred.matrix[j.idx, i] # p(G = i) assigned to class j observations
  all.preds <- c(pred.i, pred.j)
  classes <- c(rep(i, length(pred.i)), rep(j, length(pred.j)))
  o <- order(all.preds)
  classes.o <- classes[o]
  # Si: sum of ranks from class i observations
  Si <- sum(which(classes.o == i))
  ni <- length(i.idx)
  nj <- length(j.idx)
  # calculate A(i|j)
  A <- (Si - ((ni * (ni + 1))/2)) / (ni * nj)
  return(A)
}

multiclass.auc <- function(pred.matrix, ref.outcome) {
  labels <- colnames(pred.matrix)
  A.ij.cond <- utils::combn(labels, 2, function(x, pred.matrix, ref.outcome) {x
    i <- x[1]
    j <- x[2]
    A.ij <- compute.A.conditional(pred.matrix, i, j, ref.outcome)
    A.ji <- compute.A.conditional(pred.matrix, j, i, ref.outcome)
    pair <- paste0(i, "/", j)
    return(c(A.ij, A.ji))
  }, simplify = FALSE, pred.matrix = pred.matrix, ref.outcome = ref.outcome)
  c <- length(labels)
  pairs <- unlist(lapply(combn(labels, 2, simplify = FALSE), function(x) paste(x, collapse = "/")))
  A.mean <- unlist(lapply(A.ij.cond, mean))
  names(A.mean) <- pairs
  A.ij.joint <- sum(unlist(A.mean))
  M <- 2 / (c * (c-1)) * A.ij.joint 
  attr(M, "pair_AUCs") <- A.mean
  return(M)
}


### 4. Chi2?? ####

chi_contrib <- function(observ, expect) {
  rowSums(((observ - expect)^2)/expect)
}

### PSEUDO R2 ####

PseudoR2 <- function (mod0, mod, to = NULL, dat = NULL) {
  
  
  D.full <- mod$minus2loglik
  L.full <- D.full/-2
  D.base <- mod0$minus2loglik
  L.base <- D.base/-2
  
  G2 <- -2 * (L.base - L.full)
  n <- nrow(mod$data$mf)
  
  edf <- mod$paramdata$npars
  
  McFadden <- 1 - (L.full/L.base)
  McFaddenAdj <- 1 - ((L.full - edf)/L.base)
  Nagelkerke <- (1 - exp((D.full - D.base)/n))/(1 - exp(-D.base/n))
  CoxSnell <- 1 - exp(-G2/n)
  
  BIC <- D.full + log(n)*edf
  
  if(!is.null(dat)) {
    to <- DescTools::Dummy(to, "full")
    y.hat.resp <- msm_pred(mod, t = dat[,1], from = dat[,2], covariates = dat[,-c(1,2)])
    Efron <- (1 - (sum(rowSums((to - y.hat.resp)^2)))/(sum(rowSums((to - apply(to,2,mean))^2))))
    biserial_cor <- c()
    for(s in (colnames(to))) {
      biserial_cor[s] <- cor(to[,s], y.hat.resp[,s])^2
    }
    
  } else {
    Efron = NA
    biserial_cor = NA
  }
  res <- c(McFadden = McFadden, McFaddenAdj = McFaddenAdj,
           CoxSnell = CoxSnell, Nagelkerke = Nagelkerke, 
           Efron = Efron, biserial_cor = biserial_cor,
           AIC = AIC(mod), logLik = L.full,
           logLik0 = L.base, G2 = G2)
  
  return(res)
}

### WALD TEST OF COVARIATE SIGNIFICANCE ####
wald <- function(coeff, cov, index = NULL, h0 = NULL)
{
  # check input is still missing
  
  if (is.null(index)) { index <- 1:length(coeff) }
  w <- length(index)
  
  if (is.null(h0)) { h0 <- rep(0, w) } else { h0 <- rep(h0, w) }
  
  l <- matrix(rep(0, length(coeff) * w), ncol = length(coeff))
  for(i in 1:w) { l[i, index[i]] <- 1 }
  
  f <- (l %*% coeff) - h0
  mat <- qr.solve(l %*% cov %*% t(l))
  stat <- t(f) %*% mat %*% f
  p <- 1 - pchisq(stat, df = w)
  
  output <- c(stat, w, p)
  names(output) <- c("chi2", "df", "p")
  output
}

