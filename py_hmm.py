# -*- coding: utf-8 -*-
"""
Simple Poisson Hidden Markov Model 

@author: Richard Glennie 
"""
import numpy as np 
import time
from scipy.stats import poisson
from scipy.optimize import minimize

def simulate_po_hmm(n, lamb, tpm, n_state):
    delta = np.linalg.solve(np.eye(n_state) - tpm + 1, np.ones(n_state))
    s = np.zeros(n, dtype="int")
    state_space = np.arange(n_state)
    s[0] = np.random.choice(state_space, size = 1, p = delta)
    for t in np.arange(1, n):
        s[t] = np.random.choice(state_space, size = 1, p = tpm[s[t - 1],])
    data = np.random.poisson(lam = lamb[s], size = n)
    return(data)
    

def convert_n2w(lamb, tpm, n_state): 
    wpar = np.log(np.concatenate(([lamb[0]], np.diff(lamb))))
    tr_tpm = np.log(tpm / tpm.diagonal())
    I = np.eye(n_state)
    wpar = np.concatenate((wpar, tr_tpm[I == 0].ravel()))   
    return(wpar)
    
def convert_w2n(wpar, n_state): 
    lamb = np.exp(wpar[0:n_state])
    lamb = np.cumsum(lamb)
    tpm = np.eye(n_state)
    tpm[tpm == 0] = np.exp(wpar[n_state:])
    tpm = tpm / tpm.sum(axis = 0)
    delta = np.linalg.solve(np.eye(n_state) - tpm + 1, np.ones(n_state))
    return([lamb, tpm, delta])
    
def calc_negllk(wpar, data, n_state): 
    par = convert_w2n(wpar, n_state)
    lamb = par[0]
    tpm = par[1]
    delta = par[2]
    n = len(data)
    P = np.ones((n, n_state))   
    for s in np.arange(n_state):      
        P[:,s] = poisson.pmf(data, lamb[s])
    llk = 0 
    phi = delta
    sum_phi = 1
    for t in np.arange(n): 
        phi = np.dot(phi, tpm)
        phi = phi * P[t,]
        sum_phi = np.sum(phi)
        llk = llk + np.log(sum_phi)
        phi = phi / sum_phi 
    return(-llk)    

def fit_hmm(data, n_state, ini_lamb, ini_tpm): 
    ini_par = convert_n2w(ini_lamb, ini_tpm, n_state) 
    mod = minimize(calc_negllk, ini_par, method = "bfgs", 
                   args = (data, n_state))
    est = convert_w2n(mod.x, n_state)
    llk = -mod.fun
    return([llk, est])
    
def do_simulation(nsims, n, n_state, lamb, tpm):
    t = np.zeros(nsims)
    for sim in np.arange(nsims):
        print("Simulation %i" % sim)
        dat = simulate_po_hmm(n, lamb, tpm, n_state)
        start = time.clock()
        mod = fit_hmm(dat, n_state, lamb, tpm)
        end = time.clock()
        t[sim] = end - start
    return(np.mean(t))
        

n_state = 2 
lamb = np.zeros(2)
lamb[0] = 0.2 
lamb[1] = 2.1
tpm = np.eye(2)
tpm[tpm == 1] = 0.95
tpm[tpm != 0.95] = 0.05

sim = do_simulation(1, 1000, n_state, lamb, tpm)
print(sim)

