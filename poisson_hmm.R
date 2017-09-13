# Functions to fit a Poisson HMM 

library(madness)

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
  wpar <- log(lambda) 
  # Transition probabilities, logit link
  tr_tpm <- log(tpm / diag(tpm))
  wpar <- c(wpar, as.vector(tr_tpm[!diag(n.states)]))
  return(wpar)
}

ConvertW2N <- function(wpar, n.states) {
  # Poisson means first, log link
  lambda <- exp(wpar[1:n.states])
  # Tpm, logit link
  tpm <- diag(n.states)
  tpm[!diag(n.states)] <- exp(wpar[-(1:n.states)])
  tpm <- tpm / apply(tpm, 1, sum)
  delta <- solve(t(diag(n.states) - tpm  + 1), rep(1, n.states))
  return(list(lambda = lambda, tpm = tpm, delta = delta))
}

CalcNegLlk <- function(wpar, data, n.states) {
 par <- ConvertW2N(wpar, n.states)
 lambda <- par$lambda
 tpm <- par$tpm
 delta <- par$delta
 n <- length(data)
 P <- matrix(sapply(1:n.states, FUN = function(s) {return(dpois(data, lambda[s]))}), 
             nr = n, nc = n.states)
 llk <- 0
 phi <- delta 
 sum.phi <- 0 
 for (t in 1:n) {
  phi <- phi %*% tpm * P[t, ]
  sum.phi <- sum(phi) 
  llk <- llk + log(sum.phi)
  phi <- phi / sum.phi 
 } 
 return(-llk)
}

CalcMadNegLlk <- function(wpar, data, n.states) {
  llk <- madness(val = 0, dvdx = matrix(0, nr = 1, nc = length(wpar)), 
                             vtag = "llk", xtag = "wpar")
  par <- ConvertW2N(wpar, n.states)
  lambda <- par$lambda
  tpm <- par$tpm
  delta <- par$delta
  n <- length(data)
  P <- matrix(sapply(1:n.states, FUN = function(s) {return(dpois(data, lambda[s]))}), 
              nr = n, nc = n.states)
  g = "wpar")
  phi <- delta 
  sum.phi <- 0 
  for (t in 1:n) {
    phi <- phi %*% tpm * P[t, ]
    sum.phi <- sum(phi) 
    llk <- llk + log(sum.phi)
    phi <- phi / sum.phi 
  } 
  llk <- -llk
  return(to_objective(llk))
}

FitPoHmm <- function(data, n.states, ini.lambda, ini.tpm, mad = FALSE) {
  ini.par <- ConvertN2W(ini.lambda, ini.tpm, n.states) 
  if (mad) {
    mod <- nlm(CalcMadNegLlk, ini.par, data = data, n.states = n.states, hessian = TRUE)
  } else {
    mod <- nlm(CalcNegLlk, ini.par, data = data, n.states = n.states, hessian = TRUE)
  }
  est <- ConvertW2N(mod$estimate, n.states)
  V <- solve(mod$hessian)  
  llk <- -mod$minimum
  return(list(llk = llk, est = est, V = V))
}



