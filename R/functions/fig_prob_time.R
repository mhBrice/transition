fig_prob_time <- function(state_from, state_to, harvest, axes = T){
  # predict on scaled data
  tp <- seq(-2.5, 2, length.out = 200)
  time_interv <- seq(-1.5, 2.5, length.out = 200)
  newdata <- data.frame(expand.grid(TP = tp, 
                                    PP = 0,
                                    delta_TPy = 0,
                                    delta_PPy = 0,
                                    stress_mecha = 0,
                                    harvest01 = as.factor(harvest),
                                    harvested = harvest,
                                    time2harvest = harvest,
                                    age_mean = 0.5408876, # mean(Trans.df$age_mean),
                                    insect = as.factor(0), 
                                    time2insect = 0,
                                    fire = as.factor(0),
                                    time2fire = 0,
                                    plantation = as.factor(0), 
                                    time_interv = time_interv))
  prob <- as.data.frame(predict(stepmod[[state_from]], newdata = newdata, type = "probs"))
  
  # unscaled data
  tp <- tp * attr(scale_tp,"scaled:scale") + attr(scale_tp,"scaled:center")
  time_interv <- time_interv * attr(scale_time,"scaled:scale") + attr(scale_time,"scaled:center")
  
  if(axes) { 
    xlab = "Time (years)" 
    ylab = "Mean temperature" 
  } else { 
    xlab = "" 
    ylab = ""
  }
  image(x=time_interv, y=tp, 
        z = t(matrix(prob[,which(names(prob)==state_to)], ncol = length(time_interv), nrow = length(tp))),
        xlab = xlab, ylab = ylab, 
        col = pal.image(100), zlim=c(0,1), ylim = rev(range(tp)), axes = axes)
  contour(x=time_interv, y=tp, 
          z = t(matrix(prob[,which(names(prob)==state_to)], ncol = length(time_interv), nrow = length(tp))),
          cex = 0.5, add = TRUE)
}
