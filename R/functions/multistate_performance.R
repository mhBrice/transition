### Performance Measures for Multi-Class Problems ####

# from https://www.datascienceblog.net/post/machine-learning/performance-measures-multi-class-problems/

### Prediction ###

msm_pred <- function(mod, covariates = NULL, t, from, ci = "none",
                      states = c("Boreal", "Mixed", "Pioneer", "Temperate")) {
  
  n_trans <- length(from)
  
  # Predicted probability
  p_ls <- list()  

  for(i in 1:n_trans) {
    pmat <- pmatrix.msm(x = mod, t = t[i], covariates = as.list(covariates[i,]), ci = ci)
    st_from <- which(states == from[i])
    pvec <- pmat[st_from,]
    p_ls[[i]] <- pvec
  }
  
  if(ci == "none") {
    p_df <- do.call(rbind, p_ls)
  } else {
    p_df <- list()
    p_df[["estimate"]] <- do.call(rbind, lapply(p_ls, function(x) x[,"estimate"]))
    p_df[["lower"]] <- do.call(rbind, lapply(p_ls, function(x) x[,"lower"]))
    p_df[["upper"]] <- do.call(rbind, lapply(p_ls, function(x) x[,"upper"]))
    
    p_df <- lapply(p_df, "colnames <-", states)
  }
  

}  

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