# st0 = "b"
# st1 = "c"
# tr_dat = grTrans[,c("From", "To", "transition")]
# clim_dat = delta_clim[,-1]
# disturb_dat = env_df[,3:5]

RegTreeTransition <- function(st0 , st1, tr_dat, clim_dat = NULL, disturb_dat = NULL, method= "class")
{
  # print(st0)
  # print("->")
  # print(st1)
  
  var_clim = names(clim_dat)
  var_disturb = names(disturb_dat)
  
  
  if(is.null(disturb_dat)) {
    expl_var = clim_dat
  } else{
    expl_var = cbind.data.frame(clim_dat, disturb_dat)
  }
  
  
  # Select by transition
  l_st0 <- which(tr_dat[,1] %in% st0)
  
  datst0 = tr_dat[l_st0, ]
  datst0$transition =  ifelse(datst0[,2] %in% st1, 1, 0) #factor(datst0$transition, levels = unique(datst0$transition)) #
  
  envst0 = expl_var[l_st0, ]
  
  datst0 = cbind.data.frame("transition" = datst0$transition, envst0)
  
  mod <- rpart(transition ~ ., data = datst0, method = method) 

  res <- summary(mod)

  name <- paste(st0, st1, sep = "->")
  
  return(list(mod = mod, res = mod, name = name, st0 = st0, st1 = st1))
}

#test
# regTree_ab <- RegTreeTransition(st0 = "b", st1 = "c",
#                                 tr_dat = grTrans[,c("From", "To", "transition")],
#                                 clim_dat = delta_clim[,-1],
#                                 disturb_dat = env_df[,3:5])

#### Figure of regression tree ####

fig_regTree <- function(mod) {
  st0 <- mod$st0
  st1 <- mod$st1
  mod <- mod$mod
  #if(mod$method = "class") extra = 106
  pfit <- prune(mod, cp = mod$cptable[which.min(mod$cptable[,"xerror"]), "CP"])
  rpart.plot(mod, main = bquote("Transition" ~ .(st0) %->% ~ .(st1)))
  #rpart.plot(pfit, type=4,                   
             #  box.palette="GnBu",
             # branch.lty = 3, shadow.col="gray", nn=TRUE,
             # main = bquote("Transition" ~ .(st0) %->% ~ .(st1)), 
             # cex.main=1)
}



fig_all_regTree <- function(tr_dat, clim_dat = NULL, disturb_dat = NULL) {
  #models
  mod_ab <- RegTreeTransition(st0 = "a", st1 = "b", 
                            tr_dat = tr_dat, 
                            clim_dat = clim_dat, 
                            disturb_dat = disturb_dat)
  mod_ba <- RegTreeTransition(st0 = "b", st1 = "a", 
                              tr_dat = tr_dat, 
                              clim_dat = clim_dat, 
                              disturb_dat = disturb_dat)
  mod_bc <- RegTreeTransition(st0 = "b", st1 = "c", 
                              tr_dat = tr_dat, 
                              clim_dat = clim_dat, 
                              disturb_dat = disturb_dat)
  mod_cb <- RegTreeTransition(st0 = "c", st1 = "b", 
                              tr_dat = tr_dat, 
                              clim_dat = clim_dat, 
                              disturb_dat = disturb_dat)
  mod_bd <- RegTreeTransition(st0 = "b", st1 = "d", 
                              tr_dat = tr_dat, 
                              clim_dat = clim_dat, 
                              disturb_dat = disturb_dat)
  mod_db <- RegTreeTransition(st0 = "d", st1 = "b", 
                              tr_dat = tr_dat, 
                              clim_dat = clim_dat, 
                              disturb_dat = disturb_dat)
  
  quartz(width = 10, height = 6)
  layout(matrix(c(1:6), ncol = 3, byrow = F))
  par(mar = c(0,0,0,0), oma= c(0,0,0,0))
  
  fig_regTree(mod_ab)
  fig_regTree(mod_ba)
  fig_regTree(mod_bc)
  fig_regTree(mod_cb)
  fig_regTree(mod_bd)
  fig_regTree(mod_db)

}


#test

# fig_all_regTree(tr_dat = grTrans[,c("From", "To", "transition")],
#                 clim_dat = delta_clim[,-1],
#                 disturb_dat = env_df[,2:5])
