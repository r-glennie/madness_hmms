# Run HMM with TMB
source("hmm_tmb.R")

data <- readRDS("data.Rds")
ini.lambda <- c(0.5, 1.0)
ini.tpm <- matrix(c(0.9, 0.1, 0.1, 0.9), nc = 2)
n.states <- 2
dat <- list(data = data, n_states = n.states)
par <- list(wpar = ConvertN2W(ini.lambda, ini.tpm, n.states))
obj <- MakeADFun(data = dat, parameters = par, DLL = "hmm_tmb")

#mod <- FitPoHmm(data, n.states, ini.lambda, ini.tpm)
