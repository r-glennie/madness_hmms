#' Simulate Poisson data and fit HMMs with R, RcppAramdillo and Template Model
#' Builder
#'
#' @param nsims number of simulated data sets to fit models to 
#' @param n sample size per data set 
#' @param true.n.states true number of states to simulate data from 
#' @param n.states number of states to be fit 
#' @param lambda true value of the Poisson means for each state 
#' @param tpm true transition probability matrix
#' @param delta true initial distribution 
#' @param ini.lambda initial value for Poisson means (NULL = automatically determined)
#' @param ini.tpm initial transition probability matrix (NULL = automatically determined)
#'
#' @return mean time takes to fit models using each method (R, ARMA, TMB)
RunSimulation <- function(nsims = 1, n = 10000, true.n.states = 2, n.states = 2, lambda = c(0.2, 2.1), 
                          tpm = matrix(c(0.965, 0.027, 0.035, 0.973), nr = 2), 
                          delta = NULL, ini.lambda = NULL, 
                          ini.tpm = NULL) {
  if (is.null(ini.lambda)) ini.lambda <- seq(min(lambda), max(lambda), length = n.states)
  if (is.null(ini.tpm)) {
   ini.tpm <- diag(n.states)
   diag(ini.tpm) <- 0.9
   ini.tpm[!diag(n.states)] <- 0.1 / (n.states - 1) 
  } 
  if (is.null(delta)) delta <- solve(t(diag(n.states) - tpm  + 1), rep(1, n.states))
  times <- matrix(0, nr = nsims, nc = 3)
  for (s in 1:nsims) {
    cat(s, "/", nsims, "\r")
    # simulate data 
    data <- SimulatePoHmm(n, lambda, tpm, n.states, delta)
    # fit HMM with R 
    times[s, 1] <- as.numeric(system.time(FitPoHmmR(data, n.states, ini.lambda, ini.tpm))[[3]])
    # fit HMM with RcppArmadillo
    times[s, 2] <- as.numeric(system.time(FitPoHmmArma(data, n.states, ini.lambda, ini.tpm))[[3]])
    # fit HMM with TMB 
    times[s, 3] <- as.numeric(system.time(FitPoHmmTmb(data, n.states, ini.lambda, ini.tpm))[[3]])
  }
  cat("\n")
  t <- colMeans(times)
  names(t) <- c("R", "ARMA", "TMB")
  return(t)
}