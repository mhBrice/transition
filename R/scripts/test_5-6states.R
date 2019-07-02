### 5 states ####

Q5  <-  rbind(c(0.7, 0.1, 0, 0.2, 0),
             c(0.1, 0.6, 0.1, 0.2, 0),
             c(0, 0.3, 0.5, 0.1, 0.1),
             c(0.2, 0.1, 0.1, 0.5, 0.1),
             c(0, 0, .1, .1, .8))
Q.crude  <- crudeinits.msm(states_num5 ~ year_measured, plot_id, data=states_ba, qmatrix=Q5)
rownames(Q5) <- colnames(Q5) <- levels(states_ba$states_ba5)

covar <- c("sTP", "CMI", "natural", "logging")
covar_p <- c("natural", "logging")

form_all <- as.formula(paste0("~ ", paste(covar, collapse = "+")))
form_p <- as.formula(paste0("~ ", paste(covar_p, collapse = "+")))

covariates5 =  list(
  # From boreal
  "1-2" =  form_all,
  "1-4" =  form_p,
  # From mixed_b
  "2-1" = form_all,
  "2-3" = form_all,
  "2-4" = form_p,
  # From mixed_t
  "3-2" = form_all,
  "3-4" = form_p,
  "3-5" = form_all,
  # From pioneer
  "4-1" = form_all,
  "4-2" = form_all,
  "4-3" = form_all,
  "4-5" = form_all,
  # From temperate
  "5-3" = form_all,
  "5-4" = form_p
)

msm0.5 <- msm(states_num5 ~ year_measured, subject = plot_id, data = states_ba,
              qmatrix = Q5, 
              gen.inits=TRUE,
              obstype = 1, 
              control = list(trace=1, maxit=5000, fnscale=29990),
              opt.method = "optim")

msm5.5 <- msm(states_num5 ~ year_measured, subject = plot_id, data = states_ba,
              qmatrix = Q5, 
              gen.inits=TRUE,
              obstype = 1, 
              control = list(trace=1, maxit=5000, fnscale=30990),
              opt.method = "optim", 
              covariates = covariates5)

1 - msm5.5$minus2loglik/msm0.5$minus2loglik

### 6 states ####

Q6  <-  rbind(c(0.7, 0.1, 0, 0.2, 0, 0),
              c(0.1, 0.6, 0.1, 0.2, 0, 0),
              c(0, 0.3, 0.5, 0.1, 0.1, 0),
              c(0.1, 0.1, 0.1, 0.5, 0.1, .1),
              c(0, 0, .1, .1, .8, .1),
              c(0, 0, 0, .1, .8, .1))
Q.crude  <- crudeinits.msm(states_num6 ~ year_measured, plot_id, data=states_ba, qmatrix=Q6)
rownames(Q6) <- colnames(Q6) <- levels(states_ba$states_ba6)


covariates6 =  list(
  # From boreal
  "1-2" =  form_all,
  "1-4" =  form_p,
  # From boreal_m
  "2-1" = form_all,
  "2-3" = form_all,
  "2-4" = form_p,
  # From mixed
  "3-2" = form_all,
  "3-4" = form_p,
  "3-5" = form_p,
  # From pioneer
  "4-1" = form_all,
  "4-2" = form_all,
  "4-3" = form_all,
  "4-5" = form_all,
  "4-6" = form_all,
  # From temperate_m
  "5-3" = form_all,
  "5-4" = form_p,
  "5-6" = form_all,
  # From temperate
  "6-5" = form_all,
  "6-4" = form_p
)

msm5.6 <- msm(states_num6 ~ year_measured, subject = plot_id, data = states_ba,
              qmatrix = Q6, 
              gen.inits=TRUE,
              obstype = 1, 
              control = list(trace=1, maxit=5000, fnscale=35990),
              opt.method = "optim", 
              covariates = covariates6)



temperature.trend <- function(climate_var, year_measured, plot_id) {
  n = length(unique(plot_id))
  res = matrix(NA, n, 4)
  # rownames(res) = paste('Site', unique(plot_id), sep='.')
  colnames(res) = c("slope", "Rsquare", "adjRsquare", "p-value")
  for (i in 1:n) {
    curr.plot <- unique(plot_id)[i]
    sample.year <- year_measured[plot_id == curr.plot]
    sample.clim <- climate_var[plot_id == curr.plot]
    if(length(sample.year)>1) {re <- lm(sample.clim ~ sample.year)
    res[i, 1] = summary(re)$coefficients[2, 1]
    res[i, 2] = summary(re)$r.squared
    res[i, 3] = summary(re)$adj.r.squared
    f.vector = summary(re)$fstatistic
    res[i, 4] = pf(f.vector[1], f.vector[2], f.vector[3], lower.tail = FALSE)}
    
  }
  res <- cbind.data.frame(plot_id = unique(plot_id), res)
  
  
  
}



d <- data.frame(TP=states_envba$sTP,
                Boreal=states_envba$Boreal/states_envba$tot_TB, 
                year=states_envba$year_measured, 
                id=states_envba$plot_id,
                st = states_envba$states_ba) 
d2 <- d %>% group_by(id) %>% mutate(From = first(st), To = last(st)) %>%
  slice(n()) %>% mutate(col = case_when(To == "Temperate" & From != "Temperate" ~ "red", 
                                                                   st != "Temperate" ~ "grey"))
d1 = d %>% group_by(id) %>% slice(1)

d <- na.omit(d)

ggplot(d1, aes(x=TP, y=Boreal)) +
  geom_point(colour = d2$col, alpha = .5)  + theme_bw()

arrow_data <- layer_data(p1, 2) 

directions <- diff(arrow_data$y)[seq(1,nrow(arrow_data),2)]
colors <- case_when(trend$slope < 0 ~ "red",
                    trend$slope == 0 ~ "grey50",
                    trend$slope > 0  ~ "green")


d = d %>% group_by(id) %>% slice(1)
ggplot(d, aes(x=TP, y=Boreal)) +
  geom_point(alpha = .5)  + theme_bw()
