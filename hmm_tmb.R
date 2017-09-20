# HMM with TMB
library(TMB)
compile("hmm_tmb.cpp")
dyn.load(dynlib("hmm_tmb"))

SimulatePoHmm <- function(n, lambda, tpm, n.states, delta) {
 s <- numeric(n)
 state.space <- 1:n.states
 s <- sample(state.space, 1, prob = delta)
 for (t in 2:n) s[t] <- sample(state.space, 1, prob = tpm[s[t - 1], ])
 data <- rpois(n, lambda = lambda[s])
 return(data)
}

ConvertN2W <- function(lambda, tpm, n.states) { 
  # Poisson means first, log link 
  wpar <- log(c(lambda[1], diff(lambda)))
  # Transition probabilities, logit link
  tr_tpm <- log(tpm / diag(tpm))
  wpar <- c(wpar, as.vector(tr_tpm[!diag(n.states)]))
  return(wpar)
}

ConvertW2N <- function(wpar, n.states) {
  # Poisson means first, log link
  lambda <- exp(wpar[1:n.states])
  lambda <- cumsum(lambda)
  # Tpm, logit link
  tpm <- diag(n.states)
  tpm[!diag(n.states)] <- exp(wpar[-(1:n.states)])
  tpm <- tpm / apply(tpm, 1, sum)
  delta <- solve(t(diag(n.states) - tpm  + 1), rep(1, n.states))
  return(list(lambda = lambda, tpm = tpm, delta = delta))
}

FitPoHmm <- function(data, n.states, ini.lambda, ini.tpm) {
  dat <- list(data = data, n_states = n.states)
  par <- list(wpar = ConvertN2W(ini.lambda, ini.tpm, n.states))
  obj <- MakeADFun(data = dat, parameters = par, DLL = "hmm_tmb", silent = TRUE)
  obj$hessian <- FALSE
  mod <- do.call("optim", obj)
  est <- ConvertW2N(mod$par, n.states)
  sd <- sdreport(obj)  
  llk <- -mod$value
  return(list(llk = llk, est = est, sd = sd))
}
