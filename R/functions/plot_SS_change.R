plot_SS_change <- function(mod, df, tp_mixed = NULL, dist = "logging",
                           xlab = "Temperature gradient", 
                           ylab = "State proportion at equilibrium", 
                           main = NULL, unscale = NULL,
                           axes = c(1,2)) {
  
  
  # Scale axis
  if(!is.null(unscale)) {
    x <- tp_grad * unscale[2] + unscale[1]
    tp_mixed <- tp_mixed * unscale[2] + unscale[1]
  } else {
    x <- tp_grad
  }
  
  # Compute Q matrices
  qmats <- apply(df, 1, function(x)
    qmatrix.msm(mod, covariates = as.list(x), ci = "none"))
  
  SS <- apply(qmats, 2, function(x) steady_state(qmat = matrix(x, 4, 4)))
  
  bb <- SS[1,]
  tt <- SS[4,]
  
  # Colors
  col_bb <- "#158282"
  col_tt <- "#D43650"
  col_pts <- c("grey75", "grey45", "grey15")
  
  # Empty plot
  plot0(ylim = c(0,1), xlim = range(x), xaxs = "i", frame.plot = TRUE, yaxs = "i")
  xy <- par('usr')
  
  # Polygon of ecotone
  polygon(x = c(tp_mixed, rev(tp_mixed)), y = c(0,0,1,1),
          col = alpha("grey", .15), border = NA)
  
  # Axis and labels
  axis(1, labels = F, tcl = -.3)
  axis(2, labels = F, tcl = -.3)
  if(any(axes == 1)) axis(1, cex.axis = 1, tick = F, line = -.3)
  if(any(axes == 2)) axis(2, cex.axis = 1, tick = F, las = 1, line = -.3)
  mtext(xlab, 1, line = 2.1, cex = .85)
  mtext(ylab, 2, line = 2.1, cex = .85)
  
  mtext(main, 3, line = .2, font = 2, cex = .85, adj = 0)
  mtext("Boreal", 3, line = -1.3, col = col_bb, cex = .75, adj = 0.02)
  mtext("Temperate", 3, line = -1.3, col = col_tt, cex = .75, adj = 0.98)


  # Lines
  dlev <- unique(df[,dist])
  d3 <- length(dlev)<=3
  
  intl <- list()
  
  for(i in 1:length(dlev)) {
    ll <- which(df[,dist] == dlev[i])
    
    lty <- 1
    if(d3) {
      lty <- i; a <- 1; coli <- col_pts[i]
    } else {
      a <- logging[i]; coli <- "black"
    }
    
    lines(bb[ll] ~ x, col = alpha(col_bb, a), lwd = 1.5, lty = lty)
    lines(tt[ll] ~ x, col = alpha(col_tt, a), lwd = 1.5, lty = lty)
    
    # Intersect between boreal and temperate SS curves
    int <- curve_intersect(curve1 = cbind.data.frame(x, bb[ll]), 
                           curve2 = cbind.data.frame(x, tt[ll]))
    intl[[i]] <- int
    points(int$x, 1, xpd = NA, 
           pch = 21, col = "white", bg = alpha(coli, a), cex = 1.6, lwd = .7)
    if(i==1 & d3) int0 <- int 
    if(i!=1 & d3) {
      arrows(x0 = int0$x, y0 = 1+i*.02, x1 = int$x, y1 = 1+i*.02, 
             length = 0.06, col =  col_pts[i], lwd = 1.5, xpd = NA)
    }
    
  }


  invisible(intl)
}

### BARPLOT STEADY STATE & TRANSIENT ####


barplot_index <- function(index = NULL, ylim = NULL,
                          ylab = NULL,
                          lgd = c("Observed proportion",
                                  "Minor", 
                                  "Moderate natural", 
                                  "Moderate logging", 
                                  "Major natural",
                                  "Major logging"), 
                          colss = c("grey95","grey75", 
                                    "grey45", "grey45",
                                    "grey15","grey15"),
                          border = c("black", "grey75", 
                                     "grey45", "grey45",
                                     "grey15","grey15"),
                          density = c(0, 0,
                                      0, 25,
                                      0, 25)) {
  
  if(is.null(ylim)) ylim <- range(index)
  bp <- barplot(as.matrix(index), beside = T, 
                space = c(0.18,1), plot = FALSE)
  
  barplot(as.matrix(index), beside = T, 
                space = c(0.18,1), xlim = range(bp), ylim = ylim, yaxs = "i",
                col = colss, border = NA, 
                axes = FALSE, axisnames = FALSE)
  par(lwd = .8)
  barplot(as.matrix(index), beside = T, 
          space = c(0.18,1), ylim = ylim, yaxs = "i",
          col = "white", border = border, add = TRUE,
          axes = FALSE, axisnames = FALSE, density = density)
  par(lwd = 1)
  box2(1:4)
  
  axis(2, labels = F, tcl = -.3)
  axis(2, cex.axis = 1, tick = F, las = 1, line = -.3)
  
  mtext(states, 1, at = colMeans(bp), line = .2, xpd = NA, 
        cex = .75)
  mtext(ylab, 2, line = 3, font = 2, cex = 1.8)
  
  lg <-legend("topleft", legend = lgd, cex = 1.1, 
         fill = colss, border = border, pt.cex = 1.1, pt.lwd = .9,
         bty = "n")
  legend("topleft", legend = lgd, cex = 1.1, text.col = "transparent",
         fill = "white", border = border, density = density,
         pt.cex = 1.1, pt.lwd = .9,
         bty = "n")
}




### MARKOV METRICS ALONG TEMPERATURE GRADIENT ####


index_gradient <- function(mod, df, 
                           states = c("Boreal", "Mixed", "Pioneer", "Temperate")) {
  
  # Compute Q matrices
  qmats <- apply(df, 1, function(x)
    qmatrix.msm(mod, covariates = as.list(x), ci = "none"))
  
  # Steady
  steady <- apply(qmats, 2, function(x) steady_state(qmat = matrix(x, 4, 4)))
  
  ### Convergence to steady state - Damping ratio
  eig <- apply(qmats, 2, function(x) eigen(t(matrix(x, 4, 4))))
  lambda <- lapply(eig, function(x) sort(x$values, decreasing = TRUE)[2])
  damping <- lapply(lambda, function(x) exp(abs(x)))
  
  ### Half life to steady state
  halflife <- lapply(damping, function(x) log(2)/log(x)) 
  halflife <- unlist(halflife)
  
  ### Sojourn time = turnover time from Hill
  
  sojs <- apply(df, 1, function(x)
    sojourn.msm(mod, covariates = as.list(x), ci = "none")[[1]])
  
  soj_st_contrib <- steady*sojs
  row.names(soj_st_contrib) <- states
  
  soj_comm <- colSums(soj_st_contrib)
  
  ### Entropy using Jump matrix
  jump <- apply(df, 1, function(x)
    pnext.msm(msm_glb, covariates = as.list(x), ci = "none")[[1]])
  
  log_jump <- log(jump)
  log_jump[which(log_jump=="-Inf")] <- 0
  
  entropy_st <- jump*log_jump
  
  entropy_st <- apply(entropy_st, 2, function(x) -rowSums(matrix(x, 4, 4)))
  
  entropy_st_contrib <- entropy_st*steady
  row.names(entropy_st_contrib) <- states
  
  entropy_comm <- colSums(entropy_st_contrib)
  
  entropy_commrel <- entropy_comm/log(length(states)-1)
  
  
  res <- list(halflife = halflife, 
              entropy_comm = entropy_commrel, 
              soj_comm = soj_comm,
              entropy_st_contrib = entropy_st_contrib,
              soj_st_contrib = soj_st_contrib)
}


plot_transients <- function(index, tp_grad, tp_mixed = NULL, 
                            xlab = NULL, ylab = NULL, main = NULL,
                            unscale = NULL,
                            axes = c(1,2), ylim = NULL) {
  
  
  # Scale axis
  if(!is.null(unscale)) {
    tp_grad <- tp_grad * unscale[2] + unscale[1]
    tp_mixed <- tp_mixed * unscale[2] + unscale[1]
  } 
  # Empty plot
  if(is.null(ylim)) ylim <- range(index)
  plot0(ylim = ylim, xlim = range(tp_grad), xaxs = "i", yaxs = "i", 
        frame.plot = TRUE, xpd = NA)
  xy <- par('usr')
  
  # Polygon of ecotone
  polygon(x = c(tp_mixed, rev(tp_mixed)), y = c(0,0,max(ylim),max(ylim)),
          col = alpha("grey", .15), border = NA)
  
  # Axis
  axis(1, labels = F, tcl = -.3, xpd = NA)
  axis(2, labels = F, tcl = -.3, xpd = NA)
  if(any(axes == 1)) axis(1, cex.axis = .9, line = -.5, tick = F, xpd = NA)
  if(any(axes == 2)) axis(2, cex.axis = .9, line = -.5, tick = F, las = 1, xpd = NA)
  mtext(xlab, 1, line = 1.8, cex = .9, xpd = NA)
  mtext(ylab, 2, line = 2, cex = .9, xpd = NA)
  
  if(!is.null(main))
    legend("topleft", legend = main, col = "transparent", bty = "n", 
           cex = 1.1, text.font = 2, inset = c(-.05, 0))
  
  
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
