fig_prob_ima <- function(state_from, state_to, harvest, axes = T){
  # predict on scaled data
  tp <- seq(-2.5, 2, length.out = 200)
  delta_tp <- seq(-2, 2, length.out = 200)
  newdata <- data.frame(expand.grid(TP = tp, 
                                    PP = 0,
                                    delta_TP = delta_tp,
                                    delta_PP = 0,
                                    stress_mecha = 0,
                                    harvest01 = as.factor(harvest),
                                    harvested = harvest,
                                    age_mean = mean(Trans.df$age_mean),
                                    insect = as.factor(0), 
                                    fire = as.factor(0),
                                    plantation = as.factor(0), 
                                    time = 10))
  prob <- as.data.frame(predict(stepmod[[state_from]], newdata = newdata, type = "probs"))
  
  # unscaled data
  tp <- tp * attr(scale_tp,"scaled:scale") + attr(scale_tp,"scaled:center")
  delta_tp <- delta_tp*attr(scale_delta_tp,"scaled:scale") + attr(scale_delta_tp,"scaled:center")
  
  if(axes) { 
    xlab = "Mean temperature" 
    ylab = "Temperature change (°C)" 
  } else { 
    xlab = "" 
    ylab = ""
  }
  image(x=tp, y=delta_tp, 
        z = (matrix(prob[,which(names(prob)==state_to)], ncol = length(tp), nrow = length(delta_tp))),
        xlab = xlab, ylab = ylab, 
        col = pal.image(100), zlim=c(0,1), axes = axes)
  contour(x=tp, y=delta_tp, 
          z = (matrix(prob[,which(names(prob)==state_to)], ncol = length(tp), nrow = length(delta_tp))),
          cex = 0.5, add = TRUE)
}
