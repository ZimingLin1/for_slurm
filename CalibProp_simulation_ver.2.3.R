#####################################################################
## A data generating process with changeable overlapping constant
#####################################################################

true_prop_score <- function(w_1, w_2, w_3, overlap){
  prob <- plogis(overlap * (sin(4*w_1) + w_2 + 2*w_3^2))
  return(prob)
}

# we can modify the value of our overlapping constant by changing overlap in our setup
# here we just let it to be 1

true_cond_expectation <- function(a, w_1, w_2, w_3, overlap){
  pscore <- true_prop_score(w_1, w_2, w_3, overlap)
  cond <- plogis(qlogis(pscore)+ cos(3*w_2) +  a * (1 + sin(3*w_1) + w_2) )
  return(cond)
}

DGP <- function(n_sample, dim_w=3, overlap){
  sample_w_1 <- runif(n_sample, min=-2, max=2)
  sample_w_2 <- runif(n_sample, min=-2, max=2)
  sample_w_3 <- runif(n_sample, min=-2, max=2)
  
  prob <- true_prop_score(sample_w_1, sample_w_2, sample_w_3, overlap)
  sample_a <- rbinom(n_sample, 1, prob)
  sample_y <- rbinom(n_sample, 1, true_cond_expectation(sample_a, sample_w_1, sample_w_2, sample_w_3, overlap)) # binary y
  
  # sample_a <- as.factor(sample_a)
  
  dat <- data.frame(w_1=sample_w_1, w_2=sample_w_2, w_3=sample_w_3, a=sample_a, y=sample_y)
  return(dat)
}


# Adaptive propensity score truncation

truncate_pscore_adaptive <- function(A, pi, min_trunc_level = 1e-8) {
  risk_function <- function(cutoff, level=1) {
    pi <- pmax(pi, cutoff) # parallel max 
    pi <- pmin(pi, 1 - cutoff) # parallel min, so pi is truncated
    alpha <- A/pi - (1-A)/(1-pi) # Riesz-representor
    alpha1 <- 1/pi
    alpha0 <- - 1/(1-pi)
    mean(alpha^2 - 2*(alpha1 - alpha0)) # return this value
  }
  cutoff <- optim(par = 0.25, fn = risk_function, method = "Brent", lower = min_trunc_level, upper = 0.5)$par 
  return(cutoff)
}

truncate_pi <- function(pi_pred, cutoff){
  pi_pred <- pmin(pi_pred, 1 - cutoff)
  pi_pred <- pmax(pi_pred, cutoff)
  return(pi_pred)
}

# Isotonic calibration for inverse propensity score

iso_reg <- function(dat, pi_est){
  
  pi1 <- pi_est
  pi0 <- 1-pi1
  A <- dat[,c("a")]
  
  calibrator_pi1 <- as.stepfun(isoreg(pi1, A))
  calibrator_pi0 <- as.stepfun(isoreg(pi0, 1-A))
  return(list(calibrator_pi1 = calibrator_pi1, calibrator_pi0 = calibrator_pi0))  
  
}

#####################################################################
## A data generating process with changeable overlapping constant
#####################################################################

get_mu_logit <- function(dat){
  train_logit_mu <- function(dat){
    fit <- glm(y~a+w_1+w_2+w_3, data = dat, family = "binomial")
    return(fit)
  }
  fit_logit_mu <- train_logit_mu(dat)
  return(fit_logit_mu)
}

get_mu_GAM <- function(dat, task_mu){
  
  cols <- paste0("w_", 1:3)
  formula_a_lin <- paste0("y~", paste0("s(", cols, ", k = 20, bs='bs',m=c(1,0))", collapse = " + "))
  
  stack_GAM_a_lin <-  Lrnr_gam$new(family = "binomial", formula = formula_a_lin)
  sl <- Pipeline$new(Lrnr_cv$new(stack_GAM_a_lin), Lrnr_cv_selector$new(loss_squared_error))
  
  # set.seed(1917)
  sl_fit <- sl$train(task = task_mu)
  return(sl_fit)
}


get_mu_earth <- function(dat, task_mu){
  
  lrn_earth <- Lrnr_earth$new(degree=1)
  sl <- Pipeline$new(Lrnr_cv$new(lrn_earth), Lrnr_cv_selector$new(loss_squared_error))
  
  # set.seed(1917)
  sl_fit <- sl$train(task = task_mu)
  return(sl_fit)
}


get_mu_RF <- function(dat, task_mu){
  
  lrn_ranger <- Lrnr_ranger$new()
  sl <- Pipeline$new(Lrnr_cv$new(lrn_ranger), Lrnr_cv_selector$new(loss_squared_error))
  
  # set.seed(1917)
  sl_fit <- sl$train(task = task_mu)
  return(sl_fit)
}


get_mu_xgb <- function(dat, task_mu){
  
  lrn_xgb <- Lrnr_xgboost$new(max_depth = 4, min_child_weight = 15)
  sl <- Pipeline$new(Lrnr_cv$new(lrn_xgb), Lrnr_cv_selector$new(loss_squared_error))
  
  # set.seed(1917)
  sl_fit <- sl$train(task = task_mu)
  return(sl_fit)
}

# convert the predicted probabilities to categories
prob_to_cat <- function(prob){
  ifelse(prob>=0.5, 1, 0)
}


# train a single model mu for outcome regression

# get_mu <- function(dat, outcome_reg_type){
#   
#   dat_sl <- dat
#   covars <- c("w_1", "w_2", "w_3", "a")
#   outcome <- "y"
#   task_mu <- make_sl3_Task(data = dat_sl, covariates = covars,
#                            outcome = outcome, outcome_type="binomial")
#   
#   x_1 <- data.frame(a=1, w_1=dat$w_1, w_2=dat$w_2, w_3=dat$w_3)
#   x_0 <- data.frame(a=0, w_1=dat$w_1, w_2=dat$w_2, w_3=dat$w_3)
#   x_a <- data.frame(a=dat$a, w_1=dat$w_1, w_2=dat$w_2, w_3=dat$w_3)
#   
#   dat_sl1 <- dat_sl
#   dat_sl0 <- dat_sl
#   
#   dat_sl1$a <- 1
#   dat_sl0$a <- 0
#     
#   task_mu1 <- make_sl3_Task(data = dat_sl1, covariates = covars,
#                            outcome = outcome, outcome_type="binomial")
#   task_mu0 <- make_sl3_Task(data = dat_sl0, covariates = covars,
#                            outcome = outcome, outcome_type="binomial")
#   
#   if(outcome_reg_type == "logit"){
#     mu <- get_mu_logit(dat)
#     diff <- mu %>% predict(x_1) %>% prob_to_cat() - mu %>% predict(x_0) %>% prob_to_cat()
#     mu_xa <- mu %>% predict(x_a) %>% prob_to_cat()
#   }
#   
#   else if(outcome_reg_type == "GAM"){
#     mu <- get_mu_GAM(dat, task_mu)
#     diff <- mu$predict(task_mu1) %>% prob_to_cat() - mu$predict(task_mu0) %>% prob_to_cat()
#     mu_xa <- mu$predict(task_mu) %>% prob_to_cat()
#   }
#   
#   else if(outcome_reg_type == "earth"){
#     mu <- get_mu_earth(dat, task_mu)
#     diff <- mu$predict(task_mu1) %>% prob_to_cat() - mu$predict(task_mu0) %>% prob_to_cat()
#     mu_xa <- mu$predict(task_mu) %>% prob_to_cat()
#   }
#   
#   else if(outcome_reg_type == "random forest"){
#     mu <- get_mu_RF(dat, task_mu)
#     diff <- mu$predict(task_mu1) - mu$predict(task_mu0)
#     mu_xa <- mu$predict(task_mu)
#   }
#   
#   else if(outcome_reg_type == "xgboost"){
#     mu <- get_mu_xgb(dat, task_mu)
#     diff <- mu$predict(task_mu1) - mu$predict(task_mu0)
#     mu_xa <- mu$predict(task_mu)
#   }
#   else if(outcome_reg_type == "oracle"){
#     diff <- true_cond_expectation(a=1, w_1=dat$w_1, w_2=dat$w_2, w_3=dat$w_3, overlap=1) - true_cond_expectation(a=0, w_1=dat$w_1, w_2=dat$w_2, w_3=dat$w_3, overlap=1)
#     mu_xa <- true_cond_expectation(a=dat$a, w_1=dat$w_1, w_2=dat$w_2, w_3=dat$w_3, overlap=1)
#   }
#   
#   return(list(diff=diff, mu_xa=mu_xa))
#   
# }

# train two separate models for outcome regression

get_mu <- function(dat, outcome_reg_type){
  
  dat_sl <- dat
  dat_sl1 <- dat[dat$a==1,]
  dat_sl0 <- dat[dat$a==0,]
  covars <- c("w_1", "w_2", "w_3")
  outcome <- "y"
  
  # used for training mu1
  task_mu1 <- make_sl3_Task(data = dat_sl1, covariates = covars,
                           outcome = outcome, outcome_type="binomial")
  # used for training mu0
  task_mu0 <- make_sl3_Task(data = dat_sl0, covariates = covars,
                            outcome = outcome, outcome_type="binomial")
  
  # used for making predictions
  task_mu <- make_sl3_Task(data = dat_sl, covariates = covars,
                           outcome = outcome, outcome_type="binomial")
  

  if(outcome_reg_type == "logit"){
    
    mu1 <- glm(y~w_1+w_2+w_3, data = dat_sl1, family = "binomial")
    mu0 <- glm(y~w_1+w_2+w_3, data = dat_sl0, family = "binomial")
    
    mu1_preds <- mu1 %>% predict(task_mu) %>% prob_to_cat()
    mu0_preds <- mu0 %>% predict(task_mu) %>% prob_to_cat()
    
    diff <- mu1_preds - mu0_preds
    mu_xa <- ifelse(dat$a, mu1_preds, mu0_preds)
    
  }
  
  else if(outcome_reg_type == "GAM"){
    
    mu1 <- get_mu_GAM(dat, task_mu1)
    mu0 <- get_mu_GAM(dat, task_mu0)
    
    mu1_preds <- mu1$predict(task_mu)%>% prob_to_cat()
    mu0_preds <- mu0$predict(task_mu)%>% prob_to_cat()
    
    diff <- mu1_preds - mu0_preds
    mu_xa <- ifelse(dat$a, mu1_preds, mu0_preds)
    
  }
  
  else if(outcome_reg_type == "earth"){
    
    mu1 <- get_mu_earth(dat, task_mu1)
    mu0 <- get_mu_earth(dat, task_mu0)
    
    mu1_preds <- mu1$predict(task_mu)%>% prob_to_cat()
    mu0_preds <- mu0$predict(task_mu)%>% prob_to_cat()
    
    diff <- mu1_preds - mu0_preds
    mu_xa <- ifelse(dat$a, mu1_preds, mu0_preds)
  }
  
  else if(outcome_reg_type == "random forest"){
    
    mu1 <- get_mu_RF(dat, task_mu1)
    mu0 <- get_mu_RF(dat, task_mu0)
    
    mu1_preds <- mu1$predict(task_mu)%>% prob_to_cat()
    mu0_preds <- mu0$predict(task_mu)%>% prob_to_cat()
    
    diff <- mu1_preds - mu0_preds
    mu_xa <- ifelse(dat$a, mu1_preds, mu0_preds)
  }
  
  else if(outcome_reg_type == "xgboost"){
    
    mu1 <- get_mu_xgb(dat, task_mu1)
    mu0 <- get_mu_xgb(dat, task_mu0)
    
    mu1_preds <- mu1$predict(task_mu)%>% prob_to_cat()
    mu0_preds <- mu0$predict(task_mu)%>% prob_to_cat()
    
    diff <- mu1_preds - mu0_preds
    mu_xa <- ifelse(dat$a, mu1_preds, mu0_preds)
  }
  
  else if(outcome_reg_type == "oracle"){
    diff <- true_cond_expectation(a=1, w_1=dat$w_1, w_2=dat$w_2, w_3=dat$w_3, overlap=1) - true_cond_expectation(a=0, w_1=dat$w_1, w_2=dat$w_2, w_3=dat$w_3, overlap=1)
    mu_xa <- true_cond_expectation(a=dat$a, w_1=dat$w_1, w_2=dat$w_2, w_3=dat$w_3, overlap=1)
  }
  
  return(list(diff=diff, mu_xa=mu_xa))
  
}


get_pi_pred_logit <- function(dat){
  
  train_logit_pi <- function(dat){
    fit <- glm(a~w_1+w_2+w_3, data = dat, family = "binomial")
    return(fit)
  }
  fit_logit_pi <- train_logit_pi(dat)
  pi_logit_pred <- fit_logit_pi %>% predict(dat, type = "response")
  
  return(pi_logit_pred)
}


get_pi_pred_GAM <- function(dat, task_pi){
  
  cols <- paste0("w_", 1:3)
  
  formula_a_lin <- paste0("a~", paste0("s(", cols, ", k = 20, bs='bs',m=c(1,0))", collapse = " + "))
  stack_GAM_a_lin <-  Lrnr_gam$new(family = "binomial", formula = formula_a_lin)
  sl <- Pipeline$new(Lrnr_cv$new(stack_GAM_a_lin), Lrnr_cv_selector$new(loss_squared_error))
  
  # set.seed(1917)
  sl_fit <- sl$train(task = task_pi)
  pi_GAM_lin_cv <- sl_fit$predict(task = task_pi)
  return(pi_GAM_lin_cv)
}


get_pi_pred_earth <- function(dat, task_pi){
  
  lrn_earth <- Lrnr_earth$new(degree=1)
  sl <- Pipeline$new(Lrnr_cv$new(lrn_earth), Lrnr_cv_selector$new(loss_squared_error))
  
  # set.seed(1917)
  sl_fit <- sl$train(task = task_pi)
  pi_earth_cv <- sl_fit$predict(task = task_pi)
  
  return(pi_earth_cv)
}


get_pi_pred_RF <- function(dat, task_pi){
  
  lrn_ranger <- Lrnr_ranger$new()
  sl <- Pipeline$new(Lrnr_cv$new(lrn_ranger), Lrnr_cv_selector$new(loss_squared_error))
  
  # set.seed(1917)
  sl_fit <- sl$train(task = task_pi)
  pi_RF_cv <- sl_fit$predict(task = task_pi)
  
  return(pi_RF_cv)
}


get_pi_pred_xgb <- function(dat, task_pi){
  
  lrn_xgb <- Lrnr_xgboost$new(max_depth = 4, min_child_weight = 15)
  sl <- Pipeline$new(Lrnr_cv$new(lrn_xgb), Lrnr_cv_selector$new(loss_squared_error))
  
  # set.seed(1917)
  sl_fit <- sl$train(task = task_pi)
  pi_xgb_cv <- sl_fit$predict(task = task_pi)
  
  return(pi_xgb_cv)
}


get_pi_est <- function(dat, ps_est_type){
  
  dat_sl <- dat
  covars <- c("w_1", "w_2", "w_3")
  outcome <- "a"
  task_pi <- make_sl3_Task(data = dat_sl, covariates = covars,
                           outcome = outcome, outcome_type="binomial")
  
  if(ps_est_type == "logit"){
    pi_est <- get_pi_pred_logit(dat)
  }
  
  else if(ps_est_type == "GAM"){
    pi_est <- get_pi_pred_GAM(dat, task_pi)
  }
  
  else if(ps_est_type == "earth"){
    pi_est <- get_pi_pred_earth(dat, task_pi)
  }
  
  else if(ps_est_type == "random forest"){
    pi_est <- get_pi_pred_RF(dat, task_pi)
  }
  
  else if(ps_est_type == "xgboost"){
    pi_est <- get_pi_pred_xgb(dat, task_pi)
  }
  else if(ps_est_type == "oracle"){
    pi_est <- true_prop_score(w_1 = dat$w_1, w_2 = dat$w_2, w_3 = dat$w_3, overlap = 1)
  }
  
  return(pi_est)
}

get_alpha_star <- function(dat, pi_est, calibration_method){
  
  if(calibration_method == "deterministic truncation"){
    cutoff <- 25/sqrt(nrow(dat)) / log(nrow(dat))
    pi1_star <- truncate_pi(pi_est, cutoff)
    pi0_star <- 1-pi1_star
  }
  
  else if(calibration_method == "adaptive truncation"){
    cutoff <- truncate_pscore_adaptive(dat$a, pi_est)
    pi1_star <- truncate_pi(pi_est, cutoff)
    pi0_star <- 1-pi1_star
  }
  
  else if(calibration_method == "isotonic regression"){
    
    calibrator_pi1 <- iso_reg(dat, pi_est)$calibrator_pi1
    calibrator_pi0 <- iso_reg(dat, pi_est)$calibrator_pi0
    
    pi1_star <- calibrator_pi1(pi_est)
    pi0_star <- calibrator_pi0(1-pi_est)
  }
  
  else if(calibration_method == "none"){
    pi1_star <- pi
    pi0_star <- 1-pi
  }
  
  A <- dat[,c("a")]
  alpha_star  <- ifelse(A==1, 1/pi1_star, - 1/pi0_star)
  return(alpha_star)
}

# for those AIPW estimators first estimate pi and then convert pi to alpha

get_ATE_AIPW_with_pi <- function(dat, outcome_reg_type, ps_est_type, calibration_method){
  
  mu <- get_mu(dat, outcome_reg_type)
  pi <- get_pi_est(dat, ps_est_type)
  alpha <- get_alpha_star(dat, pi, calibration_method)
  reg_dif <- mu$diff
  riesz_prod <- alpha * (dat$y - mu$mu_xa)
  
  est <- mean(reg_dif + riesz_prod)
  std <- sd(reg_dif + riesz_prod) 
  
  return(list(est=est, std=std))
}

################################################################################

loss_inv_xg <- function (pred, observed) {
  A <- observed
  err <- A * pred^2  - 2 * pred
  return(err)
}

# Custom objective function

objective_invprop <- function(preds, dtrain) {
  A <- xgboost::getinfo(dtrain, "label")
  grad <- 2 * (A* preds - 1)
  hess <-  2 * A
  return(list(grad = grad, hess = hess))
}

# Custom Metric
evalerror_invprop <- function(preds, dtrain) {
  A <- xgboost::getinfo(dtrain, "label")
  err <- A * preds^2  - 2 * preds
  return(list(metric = "MyError", value = mean(err)))
}

make_ranger_inv_prop <- function(max_depth) {
  
  objective <- objective_invprop
  eval_metric <- evalerror_invprop
  
  
  return(Lrnr_xgboost$new(
    nrounds=1,
    objective = objective,
    eval_metric = eval_metric,
    eta = 1,
    num_parallel_tree = 500,
    subsample = 0.63,
    colsample_bynode = 1,
    lambda = 0,
    max_depth = max_depth,
    min_child_weight = 15))
}

make_xgboost_inv_prop <- function(nrounds, max_depth, eta = 0.2, gamma = 0, min_child_weight = 15 ) {
  
  objective <- objective_invprop
  eval_metric <- evalerror_invprop
  
  return(Lrnr_xgboost$new(nrounds=nrounds,
                          max_depth = max_depth,
                          objective = objective,
                          eval_metric = eval_metric,
                          gamma = gamma,
                          min_child_weight = min_child_weight,
                          eta = eta))
}

get_alpha_pred_RF <- function(dat){
  
  dat_sl1 <- dat
  dat_sl0 <- dat
  dat_sl0$a <- 1-dat$a
  covars <- c("w_1", "w_2", "w_3")
  outcome <- "a"
  
  task_alpha1 <- make_sl3_Task(data = dat_sl1, covariates = covars,
                           outcome = outcome, outcome_type="binomial")
  task_alpha0 <- make_sl3_Task(data = dat_sl0, covariates = covars,
                               outcome = outcome, outcome_type="binomial")
  
  lrn_rf_alpha <- make_ranger_inv_prop(max_depth = 8)
  
  sl <- Pipeline$new(Lrnr_cv$new(lrn_rf_alpha), Lrnr_cv_selector$new(loss_inv_xg))
  
  # set.seed(1917)
  sl_fit_1 <- sl$train(task = task_alpha1)
  sl_fit_0 <- sl$train(task = task_alpha0)
  
  alpha1 <- sl_fit_1$predict(task = task_alpha1)
  alpha0 <- sl_fit_0$predict(task = task_alpha0)
  alpha <- dat$a*alpha1 + (1-dat$a)*alpha0
  
  return(alpha)
}


get_alpha_pred_xgb <- function(dat){
  
  dat_sl1 <- dat
  dat_sl0 <- dat
  dat_sl0$a <- 1-dat$a
  covars <- c("w_1", "w_2", "w_3")
  outcome <- "a"
  
  task_alpha1 <- make_sl3_Task(data = dat_sl1, covariates = covars,
                               outcome = outcome, outcome_type="binomial")
  task_alpha0 <- make_sl3_Task(data = dat_sl0, covariates = covars,
                               outcome = outcome, outcome_type="binomial")
  
  lrn_xgb_alpha <- make_xgboost_inv_prop(nrounds = 20, max_depth = 4, eta=0.2, min_child_weight = 15)
  
  sl <- Pipeline$new(Lrnr_cv$new(lrn_xgb_alpha), Lrnr_cv_selector$new(loss_inv_xg))
  
  # set.seed(1917)
  sl_fit_1 <- sl$train(task = task_alpha1)
  sl_fit_0 <- sl$train(task = task_alpha0)
  
  alpha1 <- sl_fit_1$predict(task = task_alpha1)
  alpha0 <- sl_fit_0$predict(task = task_alpha0)
  alpha <- dat$a*alpha1 + (1-dat$a)*alpha0
  
  return(alpha)
}

get_alpha_est <- function(dat, alpha_est_type){
  
  if(alpha_est_type == "random forest"){
    alpha_est <- get_alpha_pred_RF(dat)
  }
  else if(alpha_est_type == "xgboost"){
    alpha_est <- get_alpha_pred_xgb(dat)
  }
  
  return(alpha_est)
}

# for those AIPW estimators estimate alpha directly

get_ATE_AIPW_with_alpha <- function(dat, outcome_reg_type, alpha_est_type){
  
  mu <- get_mu(dat, outcome_reg_type)
  alpha <- get_alpha_est(dat, alpha_est_type)
  
  reg_dif <- mu$diff
  riesz_prod <- alpha * (dat$y - mu$mu_xa)
  
  est <- mean(reg_dif + riesz_prod)
  std <- sd(reg_dif + riesz_prod)
  
  return(list(est=est, std=std))
}

################################################################################
# Monte Carlo simulation
################################################################################

get_true_ATE <- function(overlap=1){
  
  dat <- DGP(1000000, overlap=1)
  # dat1 <- dat[dat$a == 1, ]
  # dat0 <- dat[dat$a == 0, ]
  
  # not
  # cond1 <- mean(dat1$y)
  # cond0 <- mean(dat0$y)
  
  # also not
  # cond1 <- mean(true_cond_expectation(dat1$a, dat1$w_1, dat1$w_2, dat1$w_3, overlap))
  # cond0 <- mean(true_cond_expectation(dat0$a, dat0$w_1, dat0$w_2, dat0$w_3, overlap))
  
  cond1 <- mean(true_cond_expectation(1, dat$w_1, dat$w_2, dat$w_3, overlap))
  cond0 <- mean(true_cond_expectation(0, dat$w_1, dat$w_2, dat$w_3, overlap))
  
  return(cond1-cond0)
}

# these two function is for model diagnostics

sim_one_ATE_est_with_pi <- function(K_times, n_samples, overlap=1, outcome_reg_type, ps_est_type, calibration_method){
  
  AIPW_est <- numeric(K_times)
  AIPW_std <- numeric(K_times)
  
  for(k in 1:K_times){
    dat <- DGP(n_sample=n_samples, overlap=overlap)
    AIPW <- get_ATE_AIPW_with_pi(dat, outcome_reg_type, ps_est_type, calibration_method)
    AIPW_est[k] <- AIPW$est
    AIPW_std[k] <- AIPW$std
  }
  
  true_ATE <- get_true_ATE()
  
  bias <- abs(mean(AIPW_est - true_ATE)) 
  std <- sd(AIPW_est)
  mse <- bias^2+std^2
  Multiplier_CI <-  qnorm(1-{1-0.95}/2)
  se <- AIPW_std/sqrt(n_samples) 

  lower_conf_bd <- AIPW_est - Multiplier_CI * se
  upper_conf_bd <- AIPW_est + Multiplier_CI * se
  cvrg_prob <-  mean((lower_conf_bd < true_ATE) & (upper_conf_bd > true_ATE))
  
  return(list(bias=bias, std=std, mse=mse, cvrg_prob=cvrg_prob,
              AIPW_est=AIPW_est, AIPW_std=AIPW_std, se=se,
              lower_conf_bd=lower_conf_bd, upper_conf_bd=upper_conf_bd,
              true_ATE=true_ATE))  
}

sim_one_ATE_est_with_alpha <- function(K_times, n_samples, overlap=1, outcome_reg_type, alpha_est_type){
  
  AIPW_est <- numeric(K_times)
  AIPW_std <- numeric(K_times)
  
  for(k in 1:K_times){
    dat <- DGP(n_sample=n_samples, overlap=overlap)
    AIPW <- get_ATE_AIPW_with_alpha(dat, outcome_reg_type, alpha_est_type)
    AIPW_est[k] <- AIPW$est
    AIPW_std[k] <- AIPW$std
  }
  
  true_ATE <- get_true_ATE()
  bias <- abs(mean(AIPW_est - true_ATE)) 
  std <- sd(AIPW_est)
  mse <- bias^2+std^2
  Multiplier_CI <-  qnorm(1-{1-0.95}/2)
  se <- AIPW_std/sqrt(n_samples) 
  
  lower_conf_bd <- AIPW_est - Multiplier_CI * se
  upper_conf_bd <- AIPW_est + Multiplier_CI * se
  cvrg_prob <-  mean((lower_conf_bd < true_ATE) & (upper_conf_bd > true_ATE))
  
  return(list(bias=bias, std=std, mse=mse, cvrg_prob=cvrg_prob,
              AIPW_est=AIPW_est, AIPW_std=AIPW_std, se=se,
              lower_conf_bd=lower_conf_bd, upper_conf_bd=upper_conf_bd,
              true_ATE=true_ATE))  
}

# go through all the calibration methods for specified outcome regression and ps estimation methods with a particular sample size

sim_some_ATE_est_with_pi <- function(K_times, n_samples, overlap=1, outcome_reg_type, ps_est_type){
  
  
  calibration_types <- c("none", "deterministic truncation", "adaptive truncation", "isotonic regression")
  df <- data.frame(matrix(ncol = length(calibration_types), nrow = K_times))
  colnames(df) <- calibration_types
  
  df_est <- df # a dataframe used to store estimated ATE
  df_std <- df # a dataframe used to store estimated std of ATE, for computing coverage probability
  
  true_ATE <- get_true_ATE()
  Multiplier_CI <-  qnorm(1-{1-0.95}/2)
  
  for(k in 1:K_times){
    
    dat <- DGP(n_sample=n_samples, overlap=overlap)
    mu <- get_mu(dat, outcome_reg_type)
    pi <- get_pi_est(dat, ps_est_type)
    
    for(calibration_method in calibration_types){
      
      alpha <- get_alpha_star(dat, pi, calibration_method)
      reg_dif <- mu$diff
      riesz_prod <- alpha * (dat$y - mu$mu_xa)
      
      df_est[k, calibration_method] <- mean(reg_dif + riesz_prod)
      df_std[k, calibration_method] <- sd(reg_dif + riesz_prod)
    }
  }
  
  evaluation_metrics <- c("bias", "standard deviation", "root MSE", "coverage probability")
  evaluation <- data.frame(matrix(ncol = length(calibration_types), nrow = length(evaluation_metrics)))
  colnames(evaluation) <- calibration_types
  rownames(evaluation) <- evaluation_metrics
  
  for(calibration_method in calibration_types){
    
    evaluation["bias", calibration_method] <- abs(mean(df_est[,calibration_method] - true_ATE)) 
    evaluation["standard deviation", calibration_method] <- sd(df_est[,calibration_method])
    evaluation["root MSE", calibration_method] <- sqrt(evaluation["bias", calibration_method]^2 + evaluation["standard deviation", calibration_method]^2)
    
    se_ATE_est <- df_std[,calibration_method]/sqrt(n_samples) 
    lower_conf_bd <- df_est[,calibration_method] - Multiplier_CI * se_ATE_est
    upper_conf_bd <- df_est[,calibration_method] + Multiplier_CI * se_ATE_est
    cvrg_prob <-  mean((lower_conf_bd < true_ATE) & (upper_conf_bd > true_ATE))
    evaluation["coverage probability", calibration_method] <- cvrg_prob
    
  }
  
  return(evaluation)
}

# for MC simulation

sim_ATE_est_with_pi_n_varying <- function(K_times, n_samples_range, overlap=1, outcome_reg_type, ps_est_type){
  
  mylist = vector("list", length = length(n_samples_range))
  cnt <- 1
  
  for(n_samples in n_samples_range){
    
    evaluation_for_n_samples <- sim_some_ATE_est_with_pi(K_times, n_samples, overlap, outcome_reg_type, ps_est_type)
    evaluation_for_n_samples <- tibble::rownames_to_column(evaluation_for_n_samples, "evaluation metric")
    evaluation_for_n_samples['n samples'] <- rep(n_samples, 4)
    
    mylist[[cnt]] <- evaluation_for_n_samples
    cnt <- cnt+1
  }
  
  evaluation <- do.call(rbind, mylist)
  
  return(evaluation)
}

# for visualization

make_plots_for_results <- function(df){
  
  df_bias <- df[df$`evaluation metric`=="bias",] %>% subset(select = -c(`evaluation metric`))
  df_std <- df[df$`evaluation metric`=="standard deviation",] %>% subset(select = -c(`evaluation metric`))
  df_rMSE <- df[df$`evaluation metric`=="root MSE",] %>% subset(select = -c(`evaluation metric`))
  df_cvrg <- df[df$`evaluation metric`=="coverage probability",] %>% subset(select = -c(`evaluation metric`))
  
  bias_plot <- df_bias %>% gather(`calibration method`, y, -`n samples`) %>% ggplot(aes(x = `n samples`, y=y, color=`calibration method`)) +
    geom_line(linewidth = 0.8) + theme_bw() +
    theme(text = element_text(size=12), axis.text.x = element_text(size = 10 , hjust = 1, vjust = 0.5))  + 
    labs(x = "Sample Size (n)", y = "bias") + scale_y_continuous()
  
  std_plot <- df_std %>% gather(`calibration method`, y, -`n samples`) %>% ggplot(aes(x = `n samples`, y=y, color=`calibration method`)) +
    geom_line(linewidth = 0.8) + theme_bw() +
    theme(text = element_text(size=12), axis.text.x = element_text(size = 10 , hjust = 1, vjust = 0.5))  + 
    labs(x = "Sample Size (n)", y = "std") + scale_y_continuous()
  
  rMSE_plot <- df_rMSE %>% gather(`calibration method`, y, -`n samples`) %>% ggplot(aes(x = `n samples`, y=y, color=`calibration method`)) +
    geom_line(linewidth = 0.8) + theme_bw() +
    theme(text = element_text(size=12), axis.text.x = element_text(size = 10 , hjust = 1, vjust = 0.5))  + 
    labs(x = "Sample Size (n)", y = "root MSE") + scale_y_continuous()
  
  cvrg_plot <- df_cvrg %>% gather(`calibration method`, y, -`n samples`) %>% ggplot(aes(x = `n samples`, y=y, color=`calibration method`)) +
    geom_line(linewidth = 0.8) + theme_bw() +
    theme(text = element_text(size=12), axis.text.x = element_text(size = 10 , hjust = 1, vjust = 0.5))  + 
    labs(x = "Sample Size (n)", y = "coverage probability") + scale_y_continuous()
  
  return(list(bias_plot=bias_plot, std_plot=std_plot, rMSE_plot=rMSE_plot, cvrg_plot=cvrg_plot))
}

# simulation for xgboost and random forest, which includes the additional types of directly estimating alpha

sim_some_ATE_est_additional <- function(K_times, n_samples, overlap=1, outcome_reg_type, ps_est_type){
  
  
  calibration_types <- c("none", "deterministic truncation", "adaptive truncation", "isotonic regression")
  df <- data.frame(matrix(ncol = length(calibration_types)+1, nrow = K_times))
  colnames(df) <- append(calibration_types, "direct")
  
  df_est <- df # a dataframe used to store estimated ATE
  df_std <- df # a dataframe used to store estimated std of ATE, for computing coverage probability
  
  true_ATE <- get_true_ATE()
  Multiplier_CI <-  qnorm(1-{1-0.95}/2)
  
  for(k in 1:K_times){
    
    dat <- DGP(n_sample=n_samples, overlap=overlap)
    mu <- get_mu(dat, outcome_reg_type)
    pi <- get_pi_est(dat, ps_est_type)
    
    for(calibration_method in calibration_types){
      
      alpha <- get_alpha_star(dat, pi, calibration_method)
      reg_dif <- mu$diff
      riesz_prod <- alpha * (dat$y - mu$mu_xa)
      
      df_est[k, calibration_method] <- mean(reg_dif + riesz_prod)
      df_std[k, calibration_method] <- sd(reg_dif + riesz_prod)
    }
    
    # new codes
    alpha_direct <- get_alpha_est(dat, ps_est_type)
    reg_dif <- mu$diff
    riesz_prod <- alpha_direct * (dat$y - mu$mu_xa)
    
    df_est[k, "direct"] <- mean(reg_dif + riesz_prod)
    df_std[k, "direct"] <- sd(reg_dif + riesz_prod)
  }
  
  evaluation_metrics <- c("bias", "standard deviation", "root MSE", "coverage probability")
  evaluation <- data.frame(matrix(ncol = length(colnames(df)), nrow = length(evaluation_metrics)))
  colnames(evaluation) <- colnames(df)
  rownames(evaluation) <- evaluation_metrics
  
  for(calibration_method in colnames(df)){
    
    evaluation["bias", calibration_method] <- abs(mean(df_est[,calibration_method] - true_ATE)) 
    evaluation["standard deviation", calibration_method] <- sd(df_est[,calibration_method])
    evaluation["root MSE", calibration_method] <- sqrt(evaluation["bias", calibration_method]^2 + evaluation["standard deviation", calibration_method]^2)
    
    se_ATE_est <- df_std[,calibration_method]/sqrt(n_samples) 
    lower_conf_bd <- df_est[,calibration_method] - Multiplier_CI * se_ATE_est
    upper_conf_bd <- df_est[,calibration_method] + Multiplier_CI * se_ATE_est
    cvrg_prob <-  mean((lower_conf_bd < true_ATE) & (upper_conf_bd > true_ATE))
    evaluation["coverage probability", calibration_method] <- cvrg_prob
    
  }
  
  return(evaluation)
}

# for MC simulation of xgboost and random forest as ps_est_type, which includes the additional types of directly estimating alpha

sim_ATE_est_additional_n_varying <- function(K_times, n_samples_range, overlap=1, outcome_reg_type, ps_est_type){
  
  mylist = vector("list", length = length(n_samples_range))
  cnt <- 1
  
  for(n_samples in n_samples_range){
    
    evaluation_for_n_samples <- sim_some_ATE_est_additional(K_times, n_samples, overlap, outcome_reg_type, ps_est_type)
    evaluation_for_n_samples <- tibble::rownames_to_column(evaluation_for_n_samples, "evaluation metric")
    evaluation_for_n_samples['n samples'] <- rep(n_samples, 4)
    
    mylist[[cnt]] <- evaluation_for_n_samples
    cnt <- cnt+1
  }
  
  evaluation <- do.call(rbind, mylist)
  
  return(evaluation)
}

