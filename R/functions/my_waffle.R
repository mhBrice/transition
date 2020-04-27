my_waffle <- function(x, nrows, ncols = 50, 
                      pal = "Dark2", cols = NULL, main = NULL, lgd = TRUE, ...) {
  nd <- length(x)
  if(is.null(cols)){
    if(nd != 1) {
      cols <- brewer.pal(ncol(x), name = pal)
    } else {
      cols <- brewer.pal(3, name = pal)[1]
    }
  } 
  
  xx <- rep(cols, times = x)
  lx <- length(xx)
  
  m <- matrix(data = c(xx, rep(NA,(nrows*ncols-length(xx)))), 
              nrow = nrows, ncol = ncols, byrow = T)
  
  o <- cbind(c(row(m)), c(col(m))) + 1
  plot0(xlim = c(0, max(o[, 2]) + 1), ylim = c(0, max(o[, 1]) + 1),
        asp = 1, xaxs = 'i', yaxs = 'i')
  
  mtext(main, 3, line = .85, cex = .8, font = 2)
  #mtext(paste("n =", sum(x)), 3, line = 0, cex = .65)
  usr <- par("usr")
  rect(xleft=o[, 2], 
       ybottom=o[, 1], 
       xright=o[, 2] + .85, 
       ytop=o[, 1] + .85, 
       col = c(m), border = NA)
  
  if(lgd) {
    lgd_pos <- which(is.na(m), arr.ind = T)[1,2]
    legend(lgd_pos+1, mean(usr[3:4]), yjust = 0.5, cex = .8,
           legend = names(x), #paste0(names(x), " ", round(x/sum(x)*100, 2), "%"), 
           fill = cols, border = NA,
           bty = "n", x.intersp = .7, y.intersp = 1,
           inset = c(0, -1), xpd = NA)
  }
  
}


sort_col <- function(x) {
  if(is.null(dim(x))) {
    x
  }else{
    x[,order(colSums(x), decreasing = TRUE)]} 
}


