# function to plot risk ratio ####

plot_risk <- function(mod, mod0=NULL, varnames=NULL, 
                      states=c("Boreal", "Mixed", "Pioneer", "Temperate"), logy=T) {
  
  # Model fit - McFadden R2
  # rsq <- round((1 - mod$minus2loglik/mod0$minus2loglik)*100, 2)
  
  
  # Get estimates, CI and p-value
  HR <- hazard.msm(mod)
  
  if(is.null(varnames)) varnames <- names(HR)

  n_var <- length(varnames)
  
  trans <- rownames(HR[[1]])
  
  trans_list <- list()
  
  for(i in 1:length(trans)) {
    coef <- lapply(HR, function(x) x[i,])
    coef <- do.call(rbind, coef)
    
    signif <- ifelse(coef[,2] <= 1 & coef[,3] <= 1 | coef[,2] >= 1 & coef[,3] >= 1, 1, 0)
     
    
    # Color
    col_pt <- rep(NA, n_var)
    
    col_pt[which(signif==1 & coef[,1]>1)] <- "dodgerblue4"
    col_pt[which(signif==0 & coef[,1]>1)] <- "#C3D3E2"
    col_pt[which(signif==1 & coef[,1]<1)] <- "#B41414"
    col_pt[which(signif==0 & coef[,1]<1)] <- "#ECC4C4"
    
    if(!(logy)) coef <- log2(coef)
    trans_list[[trans[i]]] <- cbind.data.frame(coef, signif, col_pt=col_pt)
  }

  # Layout
  mat <- matrix(c(0,1,2,0,
                  3,0,4,5,
                  6,7,0,8,
                  0,9,10,0),
                4, 4, byrow = T)
  diag(mat) <- 11:14
  mat <- rbind(mat, 15:18)
  
  layout(mat, heights = c(1,1,1,1,.5))  
  
  # Plot
  par(mar = c(.3,2,.3,.5), oma = c(0,2,1.5,0))
  for(i in trans) {
    
    tmp <- trans_list[[i]]
    
    if(logy) {
      log <- "y"
      #ylim <- c(0.05,100)
      h <- 1
      labx <- "Hazard ratio"
    } else {
      log <- ""
      #ylim <- c(-6,7)
      h <- 0
      labx <- "Hazard ratio (log base 2)"
    }
    
    ylim <- range(lapply(trans_list, function(x) range(x[,1:3])))
    
    plot(tmp$HR, log = log, ylim = ylim, col = "transparent", 
         ann=F, xaxt="n", yaxt="n",  bty = "l", frame.plot = T)
    
    abline(h = h, col = "grey65")
    
    # bars for confidence interval
    arrows(x0 = 1:n_var, 
           y0 = tmp$L,
           y1 = tmp$U, 
           angle = 90, code = 1, 
           length = 0, col = "grey65", lwd = 1.5, xpd = NA)
    
    # points
    points(tmp$HR, pch = 19, 
           col = as.character(tmp$col_pt), cex = 1.4)
    
    # axis
    if(logy) {
      axis(2, at=seq(0,100,1), tcl= -0.2, labels=F, col = "grey35")
      axis(2, at=c(.1,1,10,100), labels = c(.1,1,10,100), las=1, cex.axis = .8)
    } else {
      axis(2, tcl= -0.2, labels=F, col = "grey35")
      axis(2, las=1, cex.axis = .8)
    }
  }
  
  # Labels
  
  par(mar=c(0.5,2,0.5,0.5))
  for(st in states) {
    plot0(text=st,fill="grey95",font=2, cex=1.5)
  }

 
  par(mar = c(.5,2,.1,.5))
  for(i in 1:4) {
    plot0(x=1:n_var, y=rep(0,n_var))
    text(x=1:n_var, y=rep(1,n_var), labels = varnames, font=2, cex = .8,
         srt = 90, xpd = NA, adj = 1) 
  }

  mtext(labx, 2,  line = 0.5, cex=.8, font = 2, las=0, outer=T)
  
  mtext("From", 2, adj = .96, line = -1.5, cex=.85, font = 2, las=0, outer=T)
  
  mtext("To", 3, adj = 0.06, line = 0, cex=.85, font = 2, outer=T)
}


### plot risk ratio in column ####

plot_risk2 <- function(mod, mod0=NULL, varnames=NULL, 
                      states=c("Boreal", "Mixed", "Pioneer", "Temperate"),
                      logy = T) {
  
  # Model fit - McFadden R2
  # rsq <- round((1 - mod$minus2loglik/mod0$minus2loglik)*100, 2)
  
  
  # Get estimates, CI and p-value
  HR <- hazard.msm(mod)
  
  if(is.null(varnames)) varnames <- names(HR)
  
  n_var <- length(varnames)
  
  trans <- rownames(HR[[1]])
  
  trans_list <- list()
  
  for(i in 1:length(trans)) {
    coef <- lapply(HR, function(x) x[i,])
    coef <- do.call(rbind, coef)
    
    signif <- ifelse(coef[,2] <= 1 & coef[,3] <= 1 | coef[,2] >= 1 & coef[,3] >= 1, 1, 0)
    
    
    # Color
    col_pt <- rep(NA, n_var)
    
    col_pt[which(signif==1 & coef[,1]>1)] <- "dodgerblue4"
    col_pt[which(signif==0 & coef[,1]>1)] <- "#C3D3E2"
    col_pt[which(signif==1 & coef[,1]<1)] <- "#B41414"
    col_pt[which(signif==0 & coef[,1]<1)] <- "#ECC4C4"
    
    if(!(logy)) coef <- log2(coef)
    trans_list[[trans[i]]] <- cbind.data.frame(coef, signif, col_pt=col_pt)
  }
  
  # Layout
  mat <- matrix(c(1:20),
                10, 2, byrow = T)
  mat <- rbind(mat, c(21,0))
  
  layout(mat,widths = c(1,.5))  
  
  # Plot
  par(oma = c(0,1.5,0,0))
  for(i in trans) {
    
    tmp <- trans_list[[i]]
    
    if(logy) {
      log <- "y"
      #ylim <- c(0.05,100)
      h <- 1
      labx <- "Hazard ratio"
    } else {
      log <- ""
      #ylim <- c(-6,7)
      h <- 0
      labx <- "Hazard ratio (log base 2)"
    }
    ylim <- range(lapply(trans_list, function(x) range(x[,1:3])))
    
    par(mar = c(.5,2,.5,.5))
    plot(tmp$HR, log = log, ylim = ylim, col = "transparent", 
         ann=F, xaxt="n", yaxt="n",  bty = "l", frame.plot = T)

    abline(h=h, col = "grey65")

    
    # bars for confidence interval
    arrows(x0 = 1:n_var, 
           y0 = tmp$L,
           y1 = tmp$U, 
           angle = 90, code = 1, 
           length = 0, col = "grey65", lwd = 1.5, xpd = NA)
    
    # points
    points(tmp$HR, pch = 21, bg = as.character(tmp$col_pt), col = as.character(tmp$col_pt))
    
    # axis
    
    if(logy) {
      axis(2, at=seq(0,100,1), tcl= -0.2, labels=F, col = "grey35")
      axis(2, at=c(.1,1,10,100), labels = c(.1,1,10,100), las=1, cex.axis = .8)
    } else {
      axis(2, tcl= -0.2, labels=F, col = "grey35")
      axis(2, las=1, cex.axis = .8)
    }
    
    par(mar = c(.5,.1,.5,.1))
    plot0(text = i, fill = "grey95", cex = 1.1, font = 2)
  }
  

  par(mar = c(.5,2,.1,.5))
  plot0(x=1:n_var, y=rep(0,n_var))
  text(x=1:n_var, y=rep(1,n_var), labels = varnames, font=2, cex = .8,
       srt = 90, xpd = NA, adj = 1) 

  mtext(labx, 2,  line = 0.5, cex=.8, font = 2, las=0,outer=T)
}

### Function to plot transition matrix ####

plot_trans <- function(pmat, ci = NULL, cols = c("#FDF7F7", "red3", "#060000"), 
                       states_lab = NULL, labels = FALSE, main = NULL) {
  
  col <- colorRampPalette(cols)(200)
  
  if(!is.null(ci)) {
    ci <- matrix(paste0("(", round(ci[,,1],2), ", ", round(ci[,,2],2), ")"), 4)
  }
  
  if(is.null(states_lab)) states_lab = colnames(tr)
  
  n <- nrow(pmat)
  
  # Plot matrix
  image2(pmat, col = col, border = "white", lwd = 2)
  
  # Axis labels
  coordx <- seq(0, 1, len = n)
  coordy <- rev(coordx)
  axis(3, at = coordx, labels = states_lab, font = 2, 
       tick = FALSE, cex.axis = 1, line = -1, col.axis = "grey15")
  axis(2, at = coordy, labels = states_lab, font = 2, 
       tick = FALSE, cex.axis = 1, las = 1, line = -.8, col.axis = "grey15")
  if(labels) {
    mtext("From", 2, font = 3, at = 1.1, las = 1, line = 1.5, cex = .8)
    mtext("To", 3, font = 3, at = -.2, line = 1.3, cex = .8)
  }
  
  # Main
  mtext(main, 3, line = 1.7, font = 1, cex = .85)

  # Probabilities
  for(i in 1:n) {
    if(is.null(ci)) {
    text(x = coordx, y = coordy[i], labels = round(pmat[i,], 2), 
         col = ifelse(pmat[i,]<.5, "black", "white"), cex = 1, xpd = NA)
    } else {
      text(x = coordx, y = coordy[i]+.1, labels = round(pmat[i,], 2), 
           col = ifelse(pmat[i,]<.5, "black", "white"), cex = 1, xpd = NA)
      text(x = coordx, y = coordy[i]-.2, labels = ci[i,], 
           col = ifelse(pmat[i,]<.5, "black", "white"), cex = .8, xpd = NA)
    }
  }
  
  
}

### plot transition probabilities through time ####

plot_pmatrix <- function(mod, t = 1:40, covar = "mean", ci = "none", 
                         st_col = c("#158282", "#A1BD93","#FEAC19", "#D43650"),
                         states = c("Boreal", "Mixed", "Pioneer", "Temperate"), 
                         main = F, yaxis = T) {
  
  
  
  if(is.list(covar)) {
    p <- list()
    for(i in 1:length(covar)) {
      p_tmp <- pmatrix.msm(mod, t = t, covariates = covar[[i]], ci = ci)
      class(p_tmp) <- "array"
      p[[i]] <- p_tmp
    }
  } else {
    p <- pmatrix.msm(mod, t = t, covariates = covar, ci = ci)
    class(p) <- "array"
  }

  
  for(st_from in 1:4) {
    plot0(x = t, xlim = c(0, max(t)), ylim = c(0,1), xaxs = "i", yaxs = "i")
    
    axis(1, labels = F, tcl= -0.5, col = "grey35")
    if(st_from==4) axis(1, tick = F, cex = .9)
    
    axis(2, labels = F, tcl= -0.5, col = "grey35")
    if(yaxis) axis(2, las = 1, tick = F, cex = .9)
    #box2()
    
    if(main) mtext(paste("From", states[st_from]), 3)
    
    if(is.list(p)) {
      for(i in 1:length(p)){
        from <- p[[i]][st_from, , ]
        
        for(st_to in 1:4) {
          lty = 1:3
          lines(t, from[st_to,], col = st_col[st_to], lwd = 1.3, lty = lty[i])
        }  
      }
    } else {
      from <- p[st_from, , ]
      
      for(st_to in 1:4) {
        lines(t, from[st_to,], col = st_col[st_to], lwd = 1.3, lty = 1)
      }  
    }
    
  }

}