### CROSS-VALIDATION FUNCTIONS #####

### MAKE FORMULAS ####
make_forms <- function(covar, covar_p = NULL) {
  
  if(is.null(covar_p)) covar_p <- covar
  
  form_all <- as.formula(paste0("~ ", paste(covar, collapse = "+")))
  form_p <- as.formula(paste0("~ ", paste(covar_p, collapse = "+")))
  
  covariates =  list(
    # From boreal
    "1-2" =  form_all,
    "1-3" =  form_p,
    # From mixed
    "2-1" = form_all,
    "2-3" =  form_p,
    "2-4" = form_all,
    # From pioneer
    "3-1" = form_all,
    "3-2" = form_all,
    "3-4" = form_all,
    # From temperate
    "4-2" = form_all,
    "4-3" = form_p
  )
  return(covariates)
}

### FAST TRANSITION DATA ####
to_trans <- function(data, covar_names = NULL) {

  states_trans <- data %>% 
    rename(year1 = year_measured, from = states) %>%
    select(plot_id, year1, from, covar_names) %>% data.table()
  
  cols = c("year1","from")
  anscols = c("year2", "to")
  states_trans[, (anscols) := shift(.SD, 1, NA, "lead"), .SDcols=cols, by=plot_id]
  
  states_trans <- states_trans %>% 
    mutate(t = year2 - year1, trans = paste0(from,"-",to)) %>% 
    filter(!is.na(to)) %>% select(plot_id, t, from, to, trans, covar_names)
}

### CREATE FOLD  ####

kfold <- function(strata, k = 10, id) {
  if(is.null(dim(strata))) {
    fold <- createFolds(y = strata, k = k)
  } else {
    strata <- interaction(strata)
    fold <- createFolds(y = strata, k = k)
  }
  fold <- lapply(fold, function(x) id[x])
}

### CHECK COVARIATE DISTRIBUTION IN FOLD ####
check_fold <- function(fold, data, data_trans) {
  covar_names <- c("sTP", "CMI", "TYPEHUMUS",  "natural", "logging")
  valid <- data %>% filter(plot_id %in% fold)
  var_summ <- summary(valid[,covar_names])
  
  train_trans <- data_trans %>% filter(!(plot_id %in% fold))
  
  var_table <- list()
  for(v in covar_names) {
    tmp <- tapply(train_trans[[v]], train_trans$trans, summary)
    var_table[[v]] <- do.call(rbind, tmp)
  }
  res <- list(var_summ, var_table)
}


### CROSS VALIDATION OF MSM ####
cv_msm <- function(data, fold, Q, covar_form, covar_names, covinits = NULL) {
  cv_res <- list()
  for(i in names(fold)){
    print(paste("Cross validation of", i))
    train <- data %>% filter(!(plot_id %in% fold[[i]])) #training set
    validation <- data %>% filter(plot_id %in% fold[[i]]) #validation set
    
    #fit model on the train data
    msm_train <- msm(states_num ~ year_measured, 
                     subject = plot_id, 
                     data = train,
                     qmatrix = Q, 
                     gen.inits = TRUE,
                     obstype = 1, 
                     control = list(trace=1, 
                                    maxit=5000, 
                                    fnscale = nrow(train)-5000),
                     opt.method = "optim",
                     covariates = covar_form, 
                     covinits = covinits) 
    
    if(is.null(msm_train$covmat)) {
      print(paste(i, "Optimisation has probably not converged to the maximum likelihood - Hessian is not positive definite."))
      
      cv_res[[i]] <- list(msm_train = msm_train,
                          pred_valid = NULL)
    } else {
      valid_trans <- to_trans(validation, covar_names = covar_names)
      
      if(is.null(covar_names)) { 
        covariates <- NULL
      } else {
        covariates <- valid_trans[covar_names]
      }
      
      p_tmp <- msm_pred(mod = msm_train, 
                        covariates = covariates, 
                        t = valid_trans$t, from = valid_trans$from)
      
      
      cv_res[[i]] <- list(msm_train = msm_train,
                          pred_valid = cbind.data.frame(plot_id = valid_trans$plot_id,
                                                        from = valid_trans$from, 
                                                        to = valid_trans$to, p_tmp))
    }
    
  }
  
  return(cv_res)
}


### PREDICTION ####

msm_pred <- function(mod, covariates = NULL, t, from, 
                     states = c("Boreal", "Mixed", "Pioneer", "Temperate")) {
  
  n_trans <- length(from)
  
  # Predicted probability
  p_ls <- list()  
  
  for(i in 1:n_trans) {
    pmat <- pmatrix.msm(x = mod, t = t[i], covariates = as.list(covariates[i,]))
    st_from <- which(states == from[i])
    pvec <- pmat[st_from,]
    p_ls[[i]] <- pvec
  }
  
  # if(ci == "none") {
  p_df <- do.call(rbind, p_ls)
  # } else {
  #   p_df <- list()
  #   p_df[["estimate"]] <- do.call(rbind, lapply(p_ls, function(x) x[,"estimate"]))
  #   p_df[["lower"]] <- do.call(rbind, lapply(p_ls, function(x) x[,"lower"]))
  #   p_df[["upper"]] <- do.call(rbind, lapply(p_ls, function(x) x[,"upper"]))
  #   
  #   p_df <- lapply(p_df, "colnames <-", states)
  # }
  
  
}  
