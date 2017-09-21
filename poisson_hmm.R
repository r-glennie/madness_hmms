# Functions to fit a Poisson HMM

# Libaries and Compiled Code
## RcppArmadillo 
library(Rcpp)
library(RcppArmadillo)
sourceCpp("hmm_arma.cpp")
## Template Model Builder 
library(TMB)
compile("hmm_tmb.cpp")
dyn.load(dynlib("hmm_tmb"))

#' Simulate Poisson HMM 
#'
#' @param n number of samples 
#' @param lambda means of Poisson for each state
#' @param tpm transition probability matrix 
#' @param n.states number of states
#' @param delta initial distribution 
#'
#' @return vector of observed counts 
SimulatePoHmm <- function(n, lambda, tpm, n.states, delta) {
 s <- numeric(n)
 state.space <- 1:n.states
 s <- sample(state.space, 1, prob = delta)
 for (t in 2:n) s[t] <- sample(state.space, 1, prob = tpm[s[t - 1], ])
 data <- rpois(n, lambda = lambda[s])
 return(data)
}


#' Convert parameters from natural to working scale 
#'
#' @param lambda means of Poisson for each state 
#' @param tpm transition probability matrix
#' @param n.states number of states 
#'
#' @return working parameters 
ConvertN2W <- function(lambda, tpm, n.states) { 
  # Poisson means first, log link 
  wpar <- log(c(lambda[1], diff(lambda))) 
  # Transition probabilities, logit link
  tr_tpm <- log(tpm / diag(tpm))
  wpar <- c(wpar, as.vector(tr_tpm[!diag(n.states)]))
  return(wpar)
}

#' Convert parameters from working to natural scale 
#'
#' @param wpar working parameters 
#' @param n.states number of states 
#'
#' @return working parameters 
ConvertW2N <- function(wpar, n.states) {
  # Poisson means first, log link
  lambda <- exp(wpar[1:n.states])
  # lambda summed cumulatively to stop label-switching 
  lambda <- cumsum(lambda)
  # Tpm, logit link
  tpm <- diag(n.states)
  tpm[!diag(n.states)] <- exp(wpar[-(1:n.states)])
  tpm <- tpm / apply(tpm, 1, sum)
  delta <- try(solve(t(diag(n.states) - tpm  + 1), rep(1, n.states)), silent = TRUE)
  if ("try-error" %in% class(delta)) delta <- rep(1, n.states) / n.states
  return(list(lambda = lambda, tpm = tpm, delta = delta))
}

#' Compute HMM Negative-Log Likelihood 
#'
#' @param wpar working parameters  
#' @param data vector of observed counts 
#' @param n.states number of states 
#'
#' @return negative-log likelihood 
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
 #cat("llk: ", llk, "par:", wpar, "\n")
 return(-llk)
}

#' Fit Poisson HMM using RcppArmadillo
#'
#' @param data vector of observed counts 
#' @param n.states number of states 
#' @param ini.lambda initial values of Poisson means for each state
#' @param ini.tpm initial transition probability matrix 
#'
#' @return list of log-likelihood and estimates for fitted model 
FitPoHmmArma <- function(data, n.states, ini.lambda, ini.tpm) {
  ini.par <- ConvertN2W(ini.lambda, ini.tpm, n.states) 
  mod <- suppressMessages(optim(ini.par, C_CalcNegLlk, data = data, 
                                n_states = n.states, method = "BFGS"))
  est <- ConvertW2N(mod$par, n.states)
  llk <- -mod$value
  return(list(llk = llk, est = est))
}

#' Fit Poisson HMM using Template Model Builder 
#'
#' @param data vector of observed counts 
#' @param n.states number of states 
#' @param ini.lambda initial values of Poisson means for each state
#' @param ini.tpm initial transition probability matrix 
#'
#' @return list of log-likelihood and estimates for fitted model 
FitPoHmmTmb <- function(data, n.states, ini.lambda, ini.tpm) {
  dat <- list(data = data, n_states = n.states)
  par <- list(wpar = ConvertN2W(ini.lambda, ini.tpm, n.states))
  # create function for automatic differentiation 
  # dat must contain all DATA required by hmm_tmb.cpp 
  # par must contain all PARAMETERS required by hmm_tmb.cpp 
  obj <- MakeADFun(data = dat, parameters = par, DLL = "hmm_tmb", silent = TRUE)
  obj$hessian <- FALSE
  mod <- do.call("optim", obj)
  est <- ConvertW2N(mod$par, n.states)
  llk <- -mod$value
  return(list(llk = llk, est = est))
}

#' Fit Poisson HMM using R
#'
#' @param data vector of observed counts 
#' @param n.states number of states 
#' @param ini.lambda initial values of Poisson means for each state
#' @param ini.tpm initial transition probability matrix 
#'
#' @return list of log-likelihood and estimates for fitted model 
FitPoHmmR <- function(data, n.states, ini.lambda, ini.tpm) {
  ini.par <- ConvertN2W(ini.lambda, ini.tpm, n.states) 
  mod <- suppressMessages(optim(ini.par, CalcNegLlk, data = data, 
                                n.states = n.states, method = "BFGS"))
  est <- ConvertW2N(mod$par, n.states)
  llk <- -mod$value
  return(list(llk = llk, est = est))
}
