# Run simulation and store how long it takes to 
# fit 999 data sets using R, RcppArmadillo and Template Model Builder. 
# 
# All three models are fit to the same 999 data sets. 
#
# Results table is saved to "results.csv" 
source("poisson_hmm.R")
source("simulation.R")

# Do simulation test for two levels of sample size 
n <- c(100, 1000)
# Fit models with 2, 4, and 6 total number of states 
n.states <- c(2, 4, 6)
results <- expand.grid(n, n.states)
colnames(results) <- c("sample.size", "n.states")
results$TMB <- results$ARMA <- results$R <- 0
nsims <- 99

for (i in 1:nrow(results)) {
  results[i, 3:5] <- RunSimulation(nsims = nsims, 
                         n = results$sample.size[i],
                         n.states = results$n.states[i])
  write.csv(results, "results.csv")
}
