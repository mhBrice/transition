plot_SS_change <- function(mod, tp_grad,
                           d_grad,
                           other_covar,
                           bg_pts = c("white", "grey55", "black"), 
                           st_col = c("#158282", "#D43650"),
                           xlab = "Temperature gradient", 
                           ylab = "State proportion at equilibrium", 
                           main = NULL, unscale = NULL,
                           axes = c(1,2), tp_ecoreg = NULL,
                           col_reg = c("#D53E4F", "#FC8D59", "#FEE08B", "#99D594", "#3288BD")) {
  
 
  # Scale axis
  if(!is.null(unscale)) {
    x <- tp_grad * sc_sTP[2] + sc_sTP[1]
    tp_ecoreg <- tp_ecoreg * sc_sTP[2] + sc_sTP[1]
  } else {
    x <- tp_grad
  }
  
  # Empty plot
  plot0(ylim = c(0,1), xlim = range(x), xaxs = "i", frame.plot = TRUE)
  xy <- par('usr')


  axis(1, labels = F, tcl = -.3)
  axis(2, labels = F, tcl = -.3)
  if(any(axes == 1)) axis(1, cex.axis = 1, line = -.5, tick = F)
  if(any(axes == 2)) axis(2, cex.axis = 1, line = -.5, tick = F, las = 1)
  mtext(xlab, 1, line = 1.8, cex = .9)
  mtext(ylab, 2, line = 1.8, cex = .9)
  
  mtext(main, 3, line = -1.5, font = 2, cex = .9, adj = 0.2)
  
  # Steady state along temperature gradient
  int <- list()

  dlevel <- length(d_grad)>1
  dlen <- ifelse(dlevel, length(d_grad), length(d_grad[[1]])) 
  
  for(d in 1:dlen) {
    SS <- c()
    if(dlevel) dd <- d_grad[[d]] else dd <- lapply(d_grad, function(x) x[d])
  
    for(tp in 1:length(tp_grad)) {
      covar_ls <- c(other_covar, 
                    sTP = tp_grad[tp], 
                    dd)
      
      qtmp <- qmatrix.msm(mod, 
                          covariates = covar_ls, 
                          ci = 'none')
      
      SS <- rbind.data.frame(SS, steady_state(qmat = qtmp))
    }
    bb <- SS[,1]
    tt <- SS[,2]+SS[,4]
    
    # Intersect between boreal and mixed+temperate SS curves
    int[[d]] <- curve_intersect(curve1 = cbind.data.frame(x, bb), 
                                curve2 = cbind.data.frame(x, tt))
    
    # Plot lines
    lty <- 1
    if(dlevel) {
      lty <- d ; col <- st_col
    } else {
      col <- alpha(st_col, d_grad[[1]][d])
    }

    lines(bb ~ x, col = col[1], lty = lty, lwd = 1.5)
    lines(tt ~ x, col = col[2], lty = lty, lwd = 1.5)
    
  }
  # Ecoreg
  if(!is.null(tp_ecoreg)) {
    for(r in 1:nrow(tp_ecoreg)) {
      lines(x = c(tp_ecoreg[r,1], tp_ecoreg[r,2]), 
            y = rep(xy[4], 2), 
            col = col_reg[r], lwd = 3, lend = 2, xpd = NA)
    }
  }
  
  # Arrows of intersect change
  int <- do.call(rbind.data.frame, int)
  
  if(dlevel) {
    points(int[,1], rep(xy[4], dlen), xpd = NA, 
           pch = 21, col = "black", bg = bg_pts, cex = 2)
    
    arrows(x0 = int[1,1], y0 = xy[4]+.035, x1 = int[2,1], y1 = xy[4]+.035, 
           length = 0.06, col =  bg_pts[2], lwd = 1.5, xpd = NA)
    arrows(x0 = int[1,1], y0 = xy[4]+.06, x1 = int[3,1], y1 = xy[4]+.06, 
           length = 0.06, col =  bg_pts[3], lwd = 1.5, xpd = NA)
  } else {
    points(int[,1], rep(xy[4], dlen), xpd = NA, 
           pch = 19, col = alpha("black", d_grad[[1]]), cex = 1)
  }
  

  
  invisible(int[,1])
}

### MARKOV METRICS ALONG TEMPERATURE GRADIENT ####

index_gradient <- function(mod, covar_d, tp_grad, 
                           states = c("Boreal", "Mixed", "Pioneer", "Temperate")) {
  
  entropy_df <- halflife_df <- nbvisit_df <- soj_df <-  c()
 
  soj_st_ls <- entropy_st_ls <- nbvisit_ls <- list()
  
  for(d in 1:length(covar_d)) {
    
    entropy_commrel <- halflife <- entropy_st_contrib <- soj_st_contrib <- c()
    
    for(tp in 1:length(tp_grad)) {
      
      # Q matrix
      qmat <- qmatrix.msm(mod,
                          covariates = c(covar_d[[d]], sTP = tp_grad[tp]), 
                          ci = "none")
      
      ### Steady state from qmat
      eig <- eigen(t(qmat)) 
      pi <- t(eig$vectors[,which(abs(eig$values) < 1e-12)])
      steady <- pi/sum(pi)
      colnames(steady) <- states
      
      ### Convergence to steady state - Damping ratio
      lambda <- sort(eig$values, decreasing = TRUE)[2]
      damping <- exp(abs(lambda))
      
      # eigenval_Pmat = exp(eigenval_Qmat)
      # damping = eigenval_Pmat_1 / eigenval_Pmat_2 = 1 / eigenval_Pmat_2
      # damping = exp(eigenval_Qmat_4) / exp(eigenval_Qmat_3) = exp(eigenval_Qmat_4 - eigenval_Qmat_3) = exp(abs(eigenval_Qmat_3))
      
      
      ### Half life to steady state
      halflife <- c(halflife, log(2)/log(damping))
      
      # Jump matrix
      jump <- pnext.msm(mod, 
                        covariates = c(covar_d[[d]], sTP = tp_grad[tp]), 
                        ci = "none")[[1]]
      log_jump <- log(jump)
      log_jump[which(log_jump=="-Inf")] <- 0
      
      ### Entropy
      
      entropy_st <- -rowSums(jump*log_jump)
      
      entropy_st_contrib <- rbind(entropy_st_contrib, entropy_st*steady)
      
      entropy_comm <- sum(entropy_st*steady)
      
      entropy_commrel <- c(entropy_commrel, entropy_comm/log(length(states)-1))


      ### Sojourn time = turnover time from Hill
      soj_tmp <- sojourn.msm(mod, covariates = c(covar_d[[d]], sTP = tp_grad[tp]),
                          ci = "none")[[1]]
      
      soj_st_contrib <- rbind(soj_st_contrib, steady*soj_tmp)
 
    }

    # Mean community sojourn time
    soj_df <- cbind(soj_df, as.numeric(rowSums(soj_st_contrib)))
    
    # Community entropy
    entropy_df <- cbind(entropy_df, as.numeric(entropy_commrel))
    
    # Half life to steady state
    halflife_df <- cbind(halflife_df, as.numeric(halflife))

    # State contribution to entropy
    entropy_st_ls[[d]] <- entropy_st_contrib
    
    # State contribution to sojourn time
    soj_st_ls[[d]] <- soj_st_contrib

  }
  res <- list(entropy_df = entropy_df, 
              halflife_df = halflife_df, 
              soj_df = soj_df,
              entropy_st_ls = entropy_st_ls,
              soj_st_ls = soj_st_ls
              )
}


plot_index_gradient <- function(index, tp_grad, 
                           bg_pts = c("white", "grey55", "black"), 
                           xlab = NULL, ylab = NULL, main = NULL,
                           unscale = NULL,
                           axes = c(1,2), ylim = NULL,
                           tp_ecoreg = NULL, 
                           col_reg = c("#D53E4F", "#FC8D59", "#FEE08B", "#99D594", "#3288BD")) {

  # Scale axis
  if(!is.null(unscale)) {
    tp_grad <- tp_grad * sc_sTP[2] + sc_sTP[1]
    tp_ecoreg <- tp_ecoreg * sc_sTP[2] + sc_sTP[1]
  } 
  # Empty plot
  if(is.null(ylim)) ylim <- range(index)
  plot0(ylim = ylim, xlim = range(tp_grad), xaxs = "i", yaxs = "i", frame.plot = TRUE)
  xy <- par('usr')
  
  # Axis
  axis(1, labels = F, tcl = -.3)
  axis(2, labels = F, tcl = -.3)
  if(any(axes == 1)) axis(1, cex.axis = .9, line = -.5, tick = F)
  if(any(axes == 2)) axis(2, cex.axis = .9, line = -.5, tick = F, las = 1)
  mtext(xlab, 1, line = 1.7, cex = .8, xpd = NA)
  mtext(ylab, 2, line = 1.7, cex = .8, xpd = NA)
  
  if(!is.null(main))
    legend("topleft", legend = main, col = "transparent", bty = "n", 
           cex = 1.1, text.font = 2, inset = c(-.05, 0))
  
  # Ecoreg
  if(!is.null(tp_ecoreg)) {
    for(r in 1:nrow(tp_ecoreg)) {
      lines(x = c(tp_ecoreg[r,1], tp_ecoreg[r,2]), 
            y = rep(xy[4], 2), 
            col = col_reg[r], lwd = 3, lend = 2, xpd = NA)
    }
  }


  for(i in 1:ncol(index)) {
    lines(index[,i] ~ tp_grad, col = "black", lty = i, lwd = 1.5)
  }
}

########

plot_sojourn <- function(mod, covar_d, tp_grad, 
                           #bg_pts = c("white", "grey55", "black"), 
                           st_col = c("#158282", "#A1BD93","#FEAC19", "#D43650"),
                           xlab = "Temperature gradient", 
                           ylab = "Sojourn time", 
                           main = NULL, unscale = NULL,
                           axes = c(1,2)) {
  
  # Empty plot
  plot0(ylim = c(1,22), xlim = range(tp_grad), xaxs = "i", frame.plot = TRUE)
  xy <- par('usr')
  
  # Axis
  if(!is.null(unscale)) {
    x <- pretty(tp_grad)
    x <- round(x * sc_sTP[2] + sc_sTP[1], 1)
  } else {
    x = NULL
  }
  
  axis(1, labels = F, tcl = -.3)
  axis(2, labels = F, tcl = -.3)
  if(any(axes == 1)) axis(1, labels = x, at = pretty(tp_grad), 
                          cex.axis = .8, line = -.5, tick = F)
  if(any(axes == 2)) axis(2, cex.axis = .8, line = -.5, tick = F, las = 1)
  mtext(xlab, 1, line = 1.7)
  mtext(ylab, 2, line = 1.8)
  
  mtext(main, 3, line = -1, adj = 0.05)
  
  # Steady state along temperature gradient
  int <- list()
  
  for(d in 1:length(covar_d)) {
    soj <- c()
    for(tp in 1:length(tp_grad)) {
      sojtmp <- envisits.msm(mod, start = 2, tot = 1000,
                          covariates = c(covar_d[[d]], sTP = tp_grad[tp]), 
                          ci = "none")
      
      soj <- rbind.data.frame(soj, sojtmp)
    }
    
    # Plot lines
    for(s in 1:4) {
      lines(soj[,s] ~ tp_grad, col = st_col[s], lty = d, lwd = 1.5)
    }

  }

}

### FIND INTERSECT BETWEEN 2 CURVES ####
curve_intersect <- function(curve1, curve2) {

  colnames(curve1) <- colnames(curve2) <- c("x", "y")
  # Approximate the functional form of both curves
  curve1_f <- approxfun(curve1$x, curve1$y, rule = 2)
  curve2_f <- approxfun(curve2$x, curve2$y, rule = 2)
  
  # Calculate the intersection of curve 1 and curve 2 along the x-axis
  point_x <- uniroot(function(x) curve1_f(x) - curve2_f(x),
                     c(min(curve1$x), max(curve1$x)))$root
  
  # Find where point_x is in curve 2
  point_y <- curve2_f(point_x)
  
  return(list(x = point_x, y = point_y))
}
