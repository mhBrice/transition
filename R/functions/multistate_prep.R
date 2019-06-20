### Function to create model matrix ####

my_mm <- function(covar = list(), data) {
  res <- list()
  for(i in 1:length(covar)) {
    form <- paste0("~ ", paste(covar[[i]], collapse = "+"))
    mm <- model.matrix(as.formula(form), data)
    res[[form]] <- mm
  }
  return(res)
}


### Function to initialize parameters using nnet::multinom ####

init_param <- function(covar, covar_p=covar, dat_trans) {
  require(nnet)
  
  full_formula <- paste0("To ~ ", paste(covar, collapse = "+"))
  
  dat_trans$To <- relevel(dat_trans$To, ref = "Boreal")
  B_multi <- multinom(full_formula, data = dat_trans, 
                      subset = From=="Boreal", trace = F)
  
  dat_trans$To <- relevel(dat_trans$To, ref = "Mixed")
  M_multi <- multinom(full_formula, data = dat_trans, 
                      subset = From=="Mixed", trace = F)
  
  dat_trans$To <- relevel(dat_trans$To, ref = "Pioneer")
  P_multi <- multinom(full_formula, data = dat_trans, 
                      subset = From=="Pioneer", trace = F)
  
  dat_trans$To <- relevel(dat_trans$To, ref = "Temperate")
  T_multi <- multinom(full_formula, data = dat_trans, 
                      subset = From=="Temperate", trace = F)
  
  # Coefficients
  coef_b <- coefficients(B_multi)
  coef_m <- coefficients(M_multi)
  coef_t <- coefficients(T_multi)
  coef_p <- coefficients(P_multi)
  
  params <- list(bm = coef_b["Mixed",],
                 bp = coef_b["Pioneer", c("(Intercept)", covar_p)],
                 
                 mb = coef_m["Boreal",],
                 mp = coef_m["Pioneer", c("(Intercept)", covar_p)],
                 mt = coef_m["Temperate",],
                 
                 pb = coef_p["Boreal",],
                 pm = coef_p["Mixed",],
                 pt = coef_p["Temperate",],
                 
                 tm = coef_t["Mixed",],
                 tp = coef_t["Pioneer", c("(Intercept)", covar_p)])
  
  params <- unlist(params)
  names(params) <- gsub("\\.", "_", names(params))
  params
}

