# Run HMM_R

source("poisson_hmm.R")

data <- readRDS("data.Rds")
ini.lambda <- c(0.5, 1.0)
ini.tpm <- matrix(c(0.9, 0.1, 0.1, 0.9), nc = 2)
n.states <- 2 
mod <- FitPoHmm(data, n.states, ini.lambda, ini.tpm)

