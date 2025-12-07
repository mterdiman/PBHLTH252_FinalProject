
library(tidyverse)
library(readr)
library(zoo)
#redoing so this file can stand alone
hiv_inc <- read.csv("/Users/miraterdiman/Desktop/My Folder Folder/MPH Year 1/PBHLTH 252/FINAL PROJECT_HIV/month_hiv.csv")

hiv_inc_ts <- ts(hiv_inc$Cases, start = c(2016,1), frequency = 12)
decomp <- stl(hiv_inc_ts, s.window = "periodic")

trend <- decomp$time.series[, "trend"]
hiv_inc$trend <- as.numeric(trend)
hiv_inc$month_num <- 1:nrow(hiv_inc)

time_vals <- time(hiv_inc_ts)
hiv_inc$year <- floor(time_vals)
hiv_inc$month_in_year <- round((time_vals - hiv_inc$year) * 12 + 1)
hiv_inc$year_month <- as.Date(as.yearmon(time_vals))

hiv_inc_parse <- hiv_inc %>%
  select(Months..num., trend)

hiv_inc_parse %>%
  ggplot(aes(x = Months..num., y = trend)) +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed") +
  geom_vline(xintercept = 40, color = "blue", linetype = "dashed") +
  geom_vline(xintercept = 65, color = "red", linetype = "dashed") +
  geom_vline(xintercept = 75, color = "red", linetype = "dashed") +
  geom_rect(aes(xmin = 0, xmax = 40, ymin = -Inf, ymax = Inf), 
            fill = "lightblue", alpha = 0.4) +
  geom_rect(aes(xmin = 65, xmax = 75, ymin = -Inf, ymax = Inf), 
            fill = "pink", alpha = 0.4) + 
  geom_point() +
  labs(x = "Month", y = "Cases") +
  ggtitle("Log HIV Incidence Data 2016-2022") +
  theme_minimal()

hiv_inc_pre <- hiv_inc_parse %>%
  filter(Months..num. < 41) %>%
  rename(Month = Months..num.,
         Cases = trend)

hiv_inc_post <- hiv_inc_parse %>%
  filter(Months..num. < 76 & Months..num. > 64) %>%
  rename(Month = Months..num.,
         Cases = trend)

#installing necessary libraries
library(deSolve)
library(ggplot2)
library(coda)

#current variable names are hard to use
time_months <- hiv_inc_pre$Month

#rounding trend cases so they aren't in decimal form (not realistic)
hiv_inc_pre$Cases <- round(hiv_inc_pre$Cases, digits = 0)
IncData <- hiv_inc_pre$Cases
numPoints <- length(time_months)

#converting to years, since our other parameter estimations are in years
time_years <- time_months / 12


## Coding the SPITU ODE (same ODE as before)
HIV_ode <- function(time, state, theta){
  
  #states
  S <- state["S"]
  P <- state["P"]
  I <- state["I"]
  T <- state["T"]
  U <- state["U"]
  
  N <- S + P + I + T + U
  
  #force of infection
  lambda <- theta["c"] * (theta["beta_I"] * I + theta["beta_T"] * T) / N
  
  #ODEs --> in year units
  dS <- theta["a"] * N - lambda * S - theta["m"] * S - theta["p"] * S + theta["b"] * P
  dP <- theta["p"] * S - theta["b"] * P - theta["m"] * P
  dI <- lambda * S - theta["gamma"] * I - theta["mu"] * I - theta["m"] * I
  dT <- theta["gamma"] * I  - theta["s"] * T - theta["m"] * T
  dU <- theta["s"] * T - theta["m"] * U
  
  list(c(dS, dP, dI, dT, dU))
}


## Coding likelihood functions
#only coding Prior for estimated parameters
logPrior <- function(theta_MH) {
  beta_I <- theta_MH[["beta_I"]]
  
  logPriorbeta_I <- dlnorm(beta_I, meanlog = log(0.5), sdlog = 1, log = TRUE)
  
  return(logPriorbeta_I)
}

# likelihood function for single data point
pointLogLike <- function(i, IncData, IncModel) {
  # Incidence is observed through a Poisson process.
  poissonLike <- dpois(x=IncData[i], lambda=IncModel[i], log=TRUE)
  if (is.na(poissonLike)) {
    return(-Inf)
  } else {
    return(poissonLike)
  }
}

# Likelihood function for all data points
trajLogLike <- function(time_years, IncData, theta, initState) {
  # Solve ODE at the observation times (in YEARS)
  traj <- data.frame(ode(
    y = initState,
    times = time_years,
    func  = HIV_ode,
    parms = theta,
    method = "ode45"
  ))
  
  # Compute modelled incidence for each observation month:
  # lambda(t) * S(t) * (delta t), with  delta t = 1/12 year
  N        <- traj$S + traj$P + traj$I + traj$T + traj$U
  lambda_t <- theta["c"] * (theta["beta_I"] * traj$I + theta["beta_T"] * traj$T) / N
  IncModel <- lambda_t * traj$S * (1/12)   # expected new cases per *month*
  
  logLike <- 0
  for (i in seq_along(IncData)) {
    logLike <- logLike + pointLogLike(i, IncData, IncModel)
  }
  logLike
}

# Posterior function

logPosterior <- function(time_years, IncData, theta_MH, theta_fixed, initState) {
  theta <- c(theta_MH, theta_fixed)
  lp    <- logPrior(theta_MH)
  ll    <- trajLogLike(time_years, IncData, theta, initState)
  return(lp + ll)
}


## Setting fixed parameters
theta_fixed <- c(
  # Infectiousness on treatment: beta_T = rho * beta_I, rho fixed
  beta_T = NA_real_,      # will be filled in inside likelihood from rho*beta_I
  c      = 1.5,            # contact rate per year (adjust if given)
  a      = 0.098,         # recruitment per year
  p      = 0.2072,        # PrEP uptake per year
  b      = 1.2,           # PrEP discontinuation per year (10-month persistence)
  gamma  = 0.5,           # I -> T per year
  mu     = 0.1,           # HIV mortality untreated per year
  s      = 0.25,          # T -> U per year
  m      = 1/35           # background mortality per year
)

rho <- 0.1  # relative infectiousness on treatment: beta_T = rho * beta_I

# Initial conditions: at-risk population ~ 6 million
N0 <- 6e6
initState <- c(
  S = N0 - 50000,  # mostly susceptible
  P = 0,
  I = 50000,       # have this be one of the theta parameters (so we can vary it) (for when we are incorporating all of the data)
  T = 0,
  U = 0
)


#pre and post covid fittings

## Metropolis-Hastings Functions

logPosteriorMH <- function(MHparams) {
  beta_I <- MHparams[["beta_I"]]
  
  # fill in beta_T each time
  theta_fixed_use <- theta_fixed
  theta_fixed_use["beta_T"] <- rho * beta_I
  
  logPosterior(
    time_years = time_years,
    IncData    = IncData,
    theta_MH   = c(beta_I = beta_I),
    theta_fixed = theta_fixed_use,
    initState  = initState
  )
}

# quick sanity check
# logPosteriorMH(c(beta_I = 0.07))

mcmcMH <- function(posterior, initTheta, proposalSD, numIterations) {
  posteriorThetaCurrent <- posterior(initTheta)
  thetaCurrent <- initTheta
  samples <- initTheta
  accepted <- 0
  
  for (i in 1:numIterations) {
    thetaProposed <- rnorm(
      n    = length(thetaCurrent),
      mean = thetaCurrent,
      sd   = proposalSD
    )
    names(thetaProposed) <- names(thetaCurrent)
    
    posteriorThetaProposed <- posterior(thetaProposed)
    logAcceptance <- posteriorThetaProposed - posteriorThetaCurrent
    randNum <- runif(1)
    
    if (randNum < exp(logAcceptance)) {
      thetaCurrent <- thetaProposed
      posteriorThetaCurrent <- posteriorThetaProposed
      accepted <- accepted + 1
    }
    
    samples <- c(samples, thetaCurrent)
    cat("iteration:", i,
        "beta_I:", thetaCurrent,
        "acceptance rate:", accepted / i, "\n")
  }
  
  samples
}


set.seed(47)

mcmcTrace <- mcmcMH(
  posterior    = logPosteriorMH,
  initTheta    = c(beta_I = 0.5),   # initial guess
  proposalSD   = c(beta_I = 0.005),           # tune for ~20–40% acceptance
  numIterations = 2000
)

trace_mat <- matrix(mcmcTrace, ncol = 1, byrow = TRUE)
trace     <- mcmc(trace_mat, start = 1)
colnames(trace) <- "beta_I"

plot(trace)
summary(trace)

beta_I_samples <- as.numeric(trace)
beta_I_hat     <- mean(beta_I_samples)
beta_I_hat


## 9. Compute R0 from posterior mean beta_I ------------------------

S0 <- initState["S"]
N0 <- sum(initState)

gamma <- theta_fixed["gamma"]; mu <- theta_fixed["mu"]
s     <- theta_fixed["s"];     m  <- theta_fixed["m"]
c_c   <- theta_fixed["c"]

R0_hat <- (c_c * S0 / N0) * beta_I_hat * (
  1 / (gamma + mu + m) +
    rho * gamma / ((gamma + mu + m) * (s + m))
)

R0_hat

############################################################
## PLOT – MODEL-PREDICTED INCIDENCE VS ACTUAL INCIDENCE  ##
############################################################

# 1. Build final parameter vector using posterior beta_I
theta_final_pre <- theta_fixed
theta_final_pre["beta_T"] <- rho * beta_I_hat
theta_final_pre["beta_I"] <- beta_I_hat

# 2. Solve the ODE over your observation period (in YEARS)
traj_final <- data.frame(ode(
  y     = initState,
  times = time_years,
  func  = HIV_ode,
  parms = theta_final_pre,
  method = "ode45"
))

# 3. Compute incidence for plotting: lambda(t)*S(t)*(1/12)
N_traj <- traj_final$S + traj_final$P + traj_final$I + traj_final$T + traj_final$U

lambda_traj <- theta_final_pre["c"] *
  (theta_final_pre["beta_I"] * traj_final$I +
     theta_final_pre["beta_T"] * traj_final$T) / N_traj

model_inc_monthly <- lambda_traj * traj_final$S * (1/12)

# 4. Build a plotting data frame
plot_df <- data.frame(
  Month = hiv_inc_pre$Month,
  Observed = hiv_inc_pre$Cases,
  Model = model_inc_monthly
)

# 5. Plot
library(ggplot2)

ggplot(plot_df, aes(x = Month)) +
  geom_point(aes(y = Observed), color = "black", size = 2) +
  geom_line(aes(y = Model), color = "blue", linewidth = 1.1) +
  labs(
    title = "Observed vs Model-Predicted HIV Incidence (Posterior Mean beta i)",
    y = "Monthly HIV Incidence",
    x = "Month"
  ) +
  theme_minimal(base_size = 14)

theta_final_pre



# Post COVID Period

## Data cleaning 

#current variable names are hard to use
time_months <- hiv_inc_post$Month

#rounding trend cases so they aren't in decimal form (not realistic)
hiv_inc_post$Cases <- round(hiv_inc_post$Cases, digits = 0)
IncData <- hiv_inc_post$Cases
numPoints <- length(time_months)

#converting to years, since our other parameter estimations are in years
time_years <- time_months / 12


## Coding the SPITU ODE (same ODE as before)
HIV_ode <- function(time, state, theta){
  
  #states
  S <- state["S"]
  P <- state["P"]
  I <- state["I"]
  T <- state["T"]
  U <- state["U"]
  
  N <- S + P + I + T + U
  
  #force of infection
  lambda <- theta["c"] * (theta["beta_I"] * I + theta["beta_T"] * T) / N
  
  #ODEs --> in year units
  dS <- theta["a"] * N - lambda * S - theta["m"] * S - theta["p"] * S + theta["b"] * P
  dP <- theta["p"] * S - theta["b"] * P - theta["m"] * P
  dI <- lambda * S - theta["gamma"] * I - theta["mu"] * I - theta["m"] * I
  dT <- theta["gamma"] * I  - theta["s"] * T - theta["m"] * T
  dU <- theta["s"] * T - theta["m"] * U
  
  list(c(dS, dP, dI, dT, dU))
}


## Coding likelihood functions
#only coding Prior for estimated parameters
logPrior <- function(theta_MH) {
  beta_I <- theta_MH[["beta_I"]]
  
  logPriorbeta_I <- dlnorm(beta_I, meanlog = log(0.5), sdlog = 1, log = TRUE)
  
  return(logPriorbeta_I)
}

# likelihood function for single data point
pointLogLike <- function(i, IncData, IncModel) {
  # Incidence is observed through a Poisson process.
  poissonLike <- dpois(x=IncData[i], lambda=IncModel[i], log=TRUE)
  if (is.na(poissonLike)) {
    return(-Inf)
  } else {
    return(poissonLike)
  }
}

# Likelihood function for all data points
trajLogLike <- function(time_years, IncData, theta, initState) {
  # Solve ODE at the observation times (in YEARS)
  traj <- data.frame(ode(
    y = initState,
    times = time_years,
    func  = HIV_ode,
    parms = theta,
    method = "ode45"
  ))
  
  # Compute modelled incidence for each observation month:
  # lambda(t) * S(t) * (delta t), with delta t = 1/12 year
  N        <- traj$S + traj$P + traj$I + traj$T + traj$U
  lambda_t <- theta["c"] * (theta["beta_I"] * traj$I + theta["beta_T"] * traj$T) / N
  IncModel <- lambda_t * traj$S * (1/12)   # expected new cases per *month*
  
  logLike <- 0
  for (i in seq_along(IncData)) {
    logLike <- logLike + pointLogLike(i, IncData, IncModel)
  }
  logLike
}

# Posterior function

logPosterior <- function(time_years, IncData, theta_MH, theta_fixed, initState) {
  theta <- c(theta_MH, theta_fixed)
  lp    <- logPrior(theta_MH)
  ll    <- trajLogLike(time_years, IncData, theta, initState)
  return(lp + ll)
}


## Setting fixed parameters

theta_fixed <- c(
  # Infectiousness on treatment: beta_T = rho * beta_I, rho fixed
  beta_T = NA_real_,      # will be filled in inside likelihood from rho*beta_I
  c      = 1.5,            # contact rate per year (adjust if given)
  a      = 0.098,         # recruitment per year
  p      = 0.2072,        # PrEP uptake per year
  b      = 1.2,           # PrEP discontinuation per year (10-month persistence)
  gamma  = 0.5,           # I -> T per year
  mu     = 0.1,           # HIV mortality untreated per year
  s      = 0.25,          # T -> U per year
  m      = 1/35           # background mortality per year
)

rho <- 0.1  # relative infectiousness on treatment: beta_T = rho * beta_I

# Initial conditions: at-risk population ~ 6 million
N0 <- 6e6
initState <- c(
  S = N0 - 37000,  # mostly susceptible
  P = 0,
  I = 37000,       # have this be one of the theta parameters (so we can vary it) (for when we are incorporating all of the data)
  T = 0,
  U = 0
)


#pre and post covid fittings

## Metropolis-Hastings Functions

logPosteriorMH <- function(MHparams) {
  beta_I <- MHparams[["beta_I"]]
  
  # fill in beta_T each time
  theta_fixed_use <- theta_fixed
  theta_fixed_use["beta_T"] <- rho * beta_I
  
  logPosterior(
    time_years = time_years,
    IncData    = IncData,
    theta_MH   = c(beta_I = beta_I),
    theta_fixed = theta_fixed_use,
    initState  = initState
  )
}

# quick sanity check
# logPosteriorMH(c(beta_I = 0.07))

mcmcMH <- function(posterior, initTheta, proposalSD, numIterations) {
  posteriorThetaCurrent <- posterior(initTheta)
  thetaCurrent <- initTheta
  samples <- initTheta
  accepted <- 0
  
  for (i in 1:numIterations) {
    thetaProposed <- rnorm(
      n    = length(thetaCurrent),
      mean = thetaCurrent,
      sd   = proposalSD
    )
    names(thetaProposed) <- names(thetaCurrent)
    
    posteriorThetaProposed <- posterior(thetaProposed)
    logAcceptance <- posteriorThetaProposed - posteriorThetaCurrent
    randNum <- runif(1)
    
    if (randNum < exp(logAcceptance)) {
      thetaCurrent <- thetaProposed
      posteriorThetaCurrent <- posteriorThetaProposed
      accepted <- accepted + 1
    }
    
    samples <- c(samples, thetaCurrent)
    cat("iteration:", i,
        "beta_I:", thetaCurrent,
        "acceptance rate:", accepted / i, "\n")
  }
  
  samples
}


set.seed(47)

mcmcTrace <- mcmcMH(
  posterior    = logPosteriorMH,
  initTheta    = c(beta_I = 0.5),   # initial guess
  proposalSD   = c(beta_I = 0.005),           # tune for ~20–40% acceptance
  numIterations = 2000
)

trace_mat <- matrix(mcmcTrace, ncol = 1, byrow = TRUE)
trace     <- mcmc(trace_mat, start = 1)
colnames(trace) <- "beta_I"

plot(trace)
summary(trace)

beta_I_samples <- as.numeric(trace)
beta_I_hat     <- mean(beta_I_samples)
beta_I_hat


## 9. Compute R0 from posterior mean beta_I ------------------------

S0 <- initState["S"]
N0 <- sum(initState)

gamma <- theta_fixed["gamma"]; mu <- theta_fixed["mu"]
s     <- theta_fixed["s"];     m  <- theta_fixed["m"]
c_c   <- theta_fixed["c"]

R0_hat <- (c_c * S0 / N0) * beta_I_hat * (
  1 / (gamma + mu + m) +
    rho * gamma / ((gamma + mu + m) * (s + m))
)

R0_hat


## Estimated R0

############################################################
## PLOT – MODEL-PREDICTED INCIDENCE VS ACTUAL INCIDENCE  ##
############################################################

# 1. Build final parameter vector using posterior beta_I
theta_final_post <- theta_fixed
theta_final_post["beta_T"] <- rho * beta_I_hat
theta_final_post["beta_I"] <- beta_I_hat

# 2. Solve the ODE over your observation period (in YEARS)
traj_final <- data.frame(ode(
  y     = initState,
  times = time_years,
  func  = HIV_ode,
  parms = theta_final_post,
  method = "ode45"
))

# 3. Compute incidence for plotting: lambda(t)*S(t)*(1/12)
N_traj <- traj_final$S + traj_final$P + traj_final$I + traj_final$T + traj_final$U

lambda_traj <- theta_final_post["c"] *
  (theta_final_post["beta_I"] * traj_final$I +
     theta_final_post["beta_T"] * traj_final$T) / N_traj

model_inc_monthly <- lambda_traj * traj_final$S * (1/12)

# 4. Build a plotting data frame
plot_df <- data.frame(
  Month = hiv_inc_post$Month,
  Observed = hiv_inc_post$Cases,
  Model = model_inc_monthly
)

# 5. Plot
library(ggplot2)

ggplot(plot_df, aes(x = Month)) +
  geom_point(aes(y = Observed), color = "black", size = 2) +
  geom_line(aes(y = Model), color = "blue", linewidth = 1.1) +
  labs(
    title = "Observed vs Model-Predicted HIV Incidence (Posterior Mean beta i)",
    y = "Monthly HIV Incidence",
    x = "Month"
  ) +
  theme_minimal(base_size = 14)

theta_final_post




