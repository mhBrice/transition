### Performance Measures for Multi-Class Problems ####

# from https://www.datascienceblog.net/post/machine-learning/performance-measures-multi-class-problems/


### 1. AUC generalization from Hand and Till ####

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
  A.ij.cond <- combn(labels, 2, function(x, pred.matrix, ref.outcome) {x
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

### 2. LOGARITHMIC SCORE ####

logscore_state <- function(to, pred, states = c("Boreal", "Mixed", "Pioneer", "Temperate")) {
  tmp <- c()
  for(s in states) {
    tmp <- cbind(tmp, logscore(to[,s] ~ pred[,s]))
  }
  colnames(tmp) <- states
  tmp
}

### PLOT PERFORMANCE MEASURES ####

plot_score <- function(score_ls, ylab = "Score", 
                       mod_names = c("Baseline", "Climate", "Soil", "Disturbances", "Full"),
                       states = c("Boreal", "Mixed", "Pioneer", "Temperate"),
                       col_mod = c("grey55", "#FDAE61", "#68340e", "#D7191C",  "#780a56"),
                       srt = 0, stats = NULL, ylim = NULL, lwd = 1.3,
                       text.width = 1.9,...) 
  {
  
  ord <- order(colMeans(score_ls[["cv_msm_glb"]]), decreasing = T)
  ord_mod <- order(unlist(lapply(score_ls, mean)), decreasing = T)
  
  nst <- length(states)
  if(is.null(ylim)) ylim = range(score_ls)
  
  plot0(xlim = c(1, nst), ylim = ylim, grid.col = "grey85", ...)
  axis(2, tick = F, line = -1, las = 1, cex.axis = .8)
  mtext(ylab, 2, line = 1.5, cex = .95)
  for(cv in 1:length(score_ls)) {
    avg <- colMeans(score_ls[[cv]])[ord]
    sdev <- apply(score_ls[[cv]], 2, sd)[ord]
    
    points(avg, type = "b", pch = 21, 
           col = col_mod[cv], bg = alpha(col_mod[cv],.5), 
           cex = .9, lwd = lwd)
    arrows(1:nst, avg-sdev, 1:nst, avg+sdev, length=0, col = col_mod[cv], lwd = lwd)
  }
  myusr <- par("usr")
  if(srt!=0) { adj = 1 ; pos = NULL } else { adj = NULL ; pos = 1}
  text(1:nst, myusr[3], states[ord], cex = 0.9, xpd = NA, adj = adj, pos = pos, srt = srt)
  
  lgd=legend("topright", legend = mod_names[ord_mod], 
             pch = 21, col = col_mod[ord_mod], pt.bg = alpha(col_mod[ord_mod],.5), pt.cex = .9,
             lwd = lwd, cex = .75, bg = alpha("white", .8), box.lwd = 0, box.col="white",
             text.width = text.width)

if(!is.null(stats)) {
  text(x = lgd$text$x+.6*lgd$rect$w, y = lgd$text$y, 
         labels = format(stats[ord_mod],3), cex = .75, bty = "n", xpd = NA)
}
  
}
