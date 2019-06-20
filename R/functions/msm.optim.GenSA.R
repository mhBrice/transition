
msm.optim.GenSA <- function(p, gr, hessian, msmdata, qmodel, qcmodel, cmodel, hmodel, ...) {
  
  
  optim.args <- list(...)
  
  if (is.null(optim.args$control)) optim.args$control <- list()
  
  optim.args <- c(optim.args, list(par=p$inits, fn=lik.msm,
                                   lower = rep(-10, length(p$inits)), 
                                   upper = rep(10, length(p$inits)),
                                   msmdata=msmdata, qmodel=qmodel, qcmodel=qcmodel,
                                   cmodel=cmodel, hmodel=hmodel, paramdata=p))
  
  
  if (requireNamespace("GenSA", quietly = TRUE)){
    opt <- do.call("GenSA", optim.args)
  } else stop("\"GenSA\" package not available")
   
  if (hessian)
    opt$hessian <- optimHess(par=opt$par, fn=lik.msm, gr=gr,
                             msmdata=msmdata, qmodel=qmodel, qcmodel=qcmodel,
                             cmodel=cmodel, hmodel=hmodel, paramdata=p)

  p$lik <- opt$value
  p$params[p$optpars] <- opt$par
  p$opt <- opt
  p$
  p
}

# msm0 <- msm(states_num ~ year_measured, subject = plot_id, data = filter(states_ba, plot_id<500),
#             qmatrix = Q, 
#             gen.inits = TRUE,
#             obstype = 1,
#             opt.method = "GenSA")
