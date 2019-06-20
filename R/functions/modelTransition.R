### function for GLM of state transition 

# tr.dat = data frame containing st0 in the first column and st1 in the second column
# clim.dat = explanatory climate variables
# disturb.dat = explanatory disturbance related variables
# st0 = state at time 0
# st1 = state at time 1
# degree = 1 or 2 for clim.dat only

 # st0 = "Temperate"
 # st1 = "Mixed"
# tr_dat = grTrans[,c("From", "To")] 
# clim_dat = delta_clim[,-1]
# disturb_dat = env_df[,3:5]

modelTransition <- function(st0 , st1, tr_dat, clim_dat = NULL, disturb_dat = NULL, degree=2, name = NULL)
{
  # print(st0)
  # print("->")
  # print(st1)
  
  var_clim = names(clim_dat)
  var_disturb = names(disturb_dat)
  
  # Scaling climate variables (disturb.dat binary for now)
  clim_dat = scale(clim_dat)
  
  if(is.null(disturb_dat)) {
    expl_var = clim_dat
  } else{
    expl_var = cbind.data.frame(clim_dat, disturb_dat)
  }
   
  
  # Select by transition
  l_st0 <- which(tr_dat[,1] %in% st0)
  
  datst0 = tr_dat[l_st0, ]
  datst0$transition = ifelse(datst0[,2] %in% st1, 1, 0)
  
  envst0 = expl_var[l_st0, ]

  datst0 = cbind.data.frame("transition" = datst0$transition, envst0)
 
  if(degree == 2) { 
    # if polynomial
    f <- paste("transition ~", 
               paste0("poly(", var_clim,", degree = 2, raw=T)", collapse=" + "))
    if(!is.null(disturb_dat))  f <- paste0(f, " + ", paste0(var_disturb, collapse=" + "))
    mod <- glm(as.formula(f), 
              family = "binomial", data = datst0)
    name_coeff2 <- names(mod$coefficients)[-1]
    name_coeff2 <- cbind(name_coeff2, name_coeff2)
    name_coeff2[startsWith(name_coeff2[,2], "poly"), 2] <- paste0(rep(var_clim, each = 2), "_", 1:2)
  } else { 
    # if not
    mod <- glm(transition ~ ., 
              family = "binomial", data = datst0)  
    }
    
  
  stepMod <- stepAIC(mod, trace = 0)
  
  #print(summary(stepMod))
  
  pred <- predict(stepMod, new = datst0,"response")
  # overall performance
  R2 <- NagelkerkeR2(stepMod)$R2
  #discrimination
  perf <- performance(prediction(pred, datst0$transition), "auc")
  AUC <- perf@y.values[[1]]
  
  ## selected vars
  coeff = summary(stepMod)$coefficients
  vars = rownames(coeff)[-1]
  if(degree == 2) vars = name_coeff2[which(name_coeff2[,1] %in% vars),2]
  effect = coeff[-1,1]
  pval = coeff[-1,4]
  ranges = apply(datst0[unlist(lapply(1:ncol(datst0), 
                                      function(x)is.numeric(datst0[,x])))], 2, range)
  # print(pval)
  
  if(is.null(name)) name = paste(st0, st1, sep = "->")
  
  return(list(mod = stepMod, vars = vars, effect = effect, pval = pval, R2 = R2, AUC = AUC, 
              ranges = ranges, 
              name = name, st0 = st0, st1 = st1))
}


###################################################################
#####    figures                            #######
###################################################################

pal <- colorRampPalette(c("lightblue", "yellow", "orange"), space = "rgb")

# figure for climate models
fig_glm <- function(mod)
{
  
  temp = seq(as.numeric(mod$ranges[1, "TP"]), as.numeric(mod$ranges[2, "TP"]), length.out = 50)
  pp = seq(as.numeric(mod$ranges[1, "PP"]), as.numeric(mod$ranges[2, "PP"]), length.out = 50)
  prob = predict(mod$mod, newdata = data.frame(expand.grid(TP = temp, PP = pp) ), type = "response")
  st0 <- mod$st0
  st1 <- mod$st1
  image(x = temp, y = pp, z = matrix(prob, ncol = length(pp), nrow = length(temp)), 
        xlab = " ", ylab = " ", col = pal(12), las=1,
        main = bquote(~ .(st0) %->% ~ .(st1)))
  contour(x=temp, y=pp, z = matrix(prob, ncol = length(pp), nrow = length(temp)), add=TRUE)
  mtext("Precipitation", side=2, outer=TRUE, cex = 1, line =2, las=0)
  mtext("Temperature", side = 1, outer = TRUE, cex = 1, line = 2)
}

# coefficient plot for all variables model
fig_coef <- function(mod)
{
  st0 <- mod$st0
  st1 <- mod$st1
  coefplot2(mod$mod, varnames = mod$vars, main = bquote( ~ .(st0) %->% ~ .(st1)), 
            cex.pts=1.2, col.pts = "#257CDE")
}


fig_all_glm <- function(tr_dat, clim_dat = NULL, disturb_dat = NULL, degree=2, gr.name, graph= "image")
{
  #tr_tbl <- paste0(trDat[,])
  #models
  mod_ab <- modelTransition(st0 = gr.name[1], st1 = gr.name[2], 
                            tr_dat = tr_dat, 
                            clim_dat = clim_dat, 
                            disturb_dat = disturb_dat,
                            degree = degree)
  mod_ba <- modelTransition(st0 = gr.name[2], st1 = gr.name[1], 
                            tr_dat = tr_dat, 
                            clim_dat = clim_dat, 
                            disturb_dat = disturb_dat,
                            degree = degree)
  mod_bc <- modelTransition(st0 = gr.name[2], st1 = gr.name[3], 
                            tr_dat = tr_dat, 
                            clim_dat = clim_dat, 
                            disturb_dat = disturb_dat,
                            degree = degree)
  mod_cb <- modelTransition(st0 = gr.name[3], st1 = gr.name[2], 
                            tr_dat = tr_dat, 
                            clim_dat = clim_dat, 
                            disturb_dat = disturb_dat,
                            degree = degree)
  mod_bd <- modelTransition(st0 = gr.name[2], st1 = gr.name[4], 
                            tr_dat = tr_dat, 
                            clim_dat = clim_dat, 
                            disturb_dat = disturb_dat,
                            degree = degree)
  mod_db <- modelTransition(st0 = gr.name[4], st1 = gr.name[2], 
                            tr_dat = tr_dat, 
                            clim_dat = clim_dat, 
                            disturb_dat = disturb_dat,
                            degree = degree)
  list(mod_ab, mod_ba, mod_bc, mod_cb, mod_db, mod_bd)
 # pdf(paste("../figures/glm_transitions", name, ".pdf",sep=""), width = 15, height = 8)
  #quartz(width = 7, height = 8)
  layout(matrix(c(1:6), ncol = 2, byrow = TRUE))
  if(graph == "image") {
    fig_glm(mod_ab)
    fig_glm(mod_ba)
    fig_glm(mod_bc)
    fig_glm(mod_cb)
    fig_glm(mod_bd)
    fig_glm(mod_db)
  } 
  if(graph == "coeff") {
    fig_coef(mod_ab)
    fig_coef(mod_ba)
    fig_coef(mod_bc)
    fig_coef(mod_cb)
    fig_coef(mod_bd)
    fig_coef(mod_db)
  }
  
 # dev.off()
  
}
