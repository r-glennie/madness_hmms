# Run simulation 
source("poisson_hmm.R")
source("simulation.R")

times <- RunSimulation(nsims = 10, n = 100)
print(times)

