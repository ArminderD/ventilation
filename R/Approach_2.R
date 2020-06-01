rm(list=ls())

# Install these packages first and then load with the below code

library(boot)
library(deSolve)
library(FME)
library(ggplot2)

# Load the saved dataframe
data_df <- read.csv("./data/ClinicA.csv", stringsAsFactors = FALSE)
data_df$Room <- NULL #remove first column which only contains the room name in first cell
# format times (to give seconds)
data_df$Time <- vapply(data_df$Time, function(x) as.numeric(strptime(x, format = "%Y-%m-%dT%H:%M:%SZ")), numeric(1))
data_df$Time <- data_df$Time - data_df$Time[1]
colnames(data_df) <- c("time", "C_out", "C_in1", "C_in2", "C_in3", "I", "Adults", "Children", "Infants")
# convert ppm to litres of CO2 per litres of air (by dividing by 1x10^6)
data_df[,c("C_out", "C_in1", "C_in2", "C_in3")] <- data_df[,c("C_out", "C_in1", "C_in2", "C_in3")] / 1e6
# take mean of three C_in measurements
data_df$C_in <- (data_df$C_in1 + data_df$C_in2 + data_df$C_in3) / 3


odefunc <- function(dat, G_adult, G_child, G_infant, V) {
  # The bootstrap code needs these defined within the function
  G_adult <- 0.00377
  G_child <-0.001975
  G_infant <- 0.00105
  V <- 135.363 * 1000
  
  C_out_func <- approxfun(dat[,c("time", "C_out")], method = "linear", rule = 2) # linear interpolation
  t_vec_fine <- seq(0, max(dat[,1]), by = 60) # calculate C_out every 60 seconds
  C_out <- C_out_func(t_vec_fine) # calculate C_out every 60 seconds
  C_out <- mean(C_out) 
  
  dat$IG_adults <- dat$Adults * G_adult 
  dat$IG_children <- dat$Children * G_child 
  dat$IG_infants <- dat$Infants * G_infant 
  dat$IG  <- (dat$IG_adults + dat$IG_children + dat$IG_infants) # the total generation rate of CO2 by the people in the room
  dat$diff_C <- dat$C_in - C_out 
  
  lm1 <- lm(formula = IG ~ diff_C + 0, data = dat) # constrained linear regression with intercept = 0
  Q_guess <- unname(lm1$coefficients)
  
  ACH <- (Q_guess*3.6)/(V/1000)
  IG_func <- approxfun(dat[,c("time", "IG")], method = "linear", rule = 2) 
  
  calc_C_in <- function(Q, V, IG_func, C_out, C_in_0, t_vec) {
    model <- function(t, y, parms) {  
      with(as.list(c(y, parms)), {
        dy <- ((C_out - y) * Q + IG_func(t)) / V
        return(list(dy))
      })
    }
    
    parms <- c(Q = Q, C_out = C_out, V = V)
    
    out <- as.data.frame(ode(y = c(C_in = C_in_0), times = t_vec, func = model, parms = parms))
    return(out)
  }
  
  pars <- c(Q = Q_guess, C_in_0 = dat$C_in[1]) # initial guess for parameter values.  Note that
  # the initial inside concentration is a fitted parameter because we don't assume the first observation is perfect
  lower <- c(Q = Q_guess / 10, C_in_0 = dat$C_in[1] / 2) # lower bound for parameter search
  upper <- c(Q = Q_guess * 10, C_in_0 = dat$C_in[1] * 2) # upper bound for parameter search
  
  calc_cost_func_closure <- function(IG_func, C_out, V, dat){
    calc_cost <- function(p){
      C_in <- calc_C_in(Q = p[["Q"]], IG_func = IG_func, V = V, C_out = C_out, C_in_0 = p[["C_in_0"]], t_vec = dat$time)
      cost <- modCost(obs = dat[,c("time", "C_in")], model = C_in)
      return(cost)
    }
    calc_cost
  }
  
  calc_cost <- calc_cost_func_closure(IG_func, C_out, V, dat)
  
  fit <- modFit(p = pars, f = calc_cost, lower = lower, upper = upper,
                method = "Pseudo", control = c(npop = 50, numiter = 1000, varleft = 1e-05))
  
  # non-steady state:
  opt_par <- c(fit$par)
  opt_model <- calc_C_in(Q = opt_par[["Q"]], IG_func = IG_func, V = V, C_out = C_out,
                         C_in_0 = opt_par[["C_in_0"]], t_vec = dat$time)
  
  
  Qfit <- opt_par[["Q"]] # Non-steady state, from running fit (around line 94)
  ACHfit <- (opt_par[["Q"]]*3.6)/(V/1000)
  
  Q_ACH_res_model <- list(Q = Qfit, ACH = ACHfit, 
                          res = fit$res, model_pred = opt_model)
  return(Q_ACH_res_model)
}

# solve for the actual data
Q_ACH_res_model <- odefunc(data_df,G_adult, G_child, G_infant, V)
# get residuals
res <- Q_ACH_res_model$res
# get moel predictions
model_pred <- data_df
model_pred$C_in <- Q_ACH_res_model$model_pred$C_in
set.seed(1)
bootstrap_func <- function() {
  resampled_df <- model_pred
  # resample the residuals and add them to the model predictions
  # to get the bootstrapped data
  resampled_df$C_in <- resampled_df$C_in + sample(res, length(res), replace = TRUE)
  Q_ACH_res_model <- odefunc(resampled_df,G_adult, G_child, G_infant, V) # estimate Q and ACh from resampled data
  c(Q = Q_ACH_res_model$Q,
    ACH = Q_ACH_res_model$ACH)
}

# do bootstrap replicates
n_replicates <- 10
res <- replicate(n_replicates, bootstrap_func())

Qs <- sort(res[1,])
ACHs <- sort(res[2,])

# to get the Q and ACH from the model: Q_ACH_res_model

#to get the 95% CIs:
quantile(Qs, c(.025, .975)) 
quantile(ACHs, c(.025, .975)) 


