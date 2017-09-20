# Compare speeds of R and TMB
RunSimulation <- function(nsims = 1, n = 10000, n.states = 2, lambda = c(0.2, 0.2), 
                          tpm = matrix(c(0.95, 0.05, 0.05, 0.95), nr = 2), 
                          delta = NULL, ini.lambda = c(0.5, 1.0), 
                          ini.tpm = matrix(c(0.9, 0.1, 0.1, 0.9)), nr = 2) {
  if (is.null(delta)) delta <- solve(t(diag(n.states) - tpm  + 1), rep(1, n.states))
  times <- matrix(0, nr = nsims, nc = 2)
  for (s in 1:nsims) {
    cat(s, "/", nsims, "\r")
    # simulate data 
    data <- SimulatePoHmm(n, lambda, tpm, n.states, delta)
    # fit HMM with R 
    times[s, 1] <- as.numeric(system.time(FitPoHmm(data, n.states, ini.lambda, ini.tpm))[[3]])
    # fit HMM with TMB 
    times[s, 2] <- as.numeric(system.time(FitPoHmmTMB(data, n.states, ini.lambda, ini.tpm))[[3]])
  }
  cat("\n")
  t <- colMeans(times)
  names(t) <- c("R", "TMB")
  return(t)
}