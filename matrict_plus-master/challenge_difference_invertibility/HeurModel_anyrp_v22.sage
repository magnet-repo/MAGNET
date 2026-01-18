# MatRiCT+ invertibility bound heuristic model computation
# Ron Steinfeld, 29 July 2021
# Please cite the MatRiCT+ paper accepted at IEEE S&P 2022 conference (full version available at https://eprint.iacr.org/2021/545) when using scripts

import numpy as np
import time
import csv
from sage.modules.free_module_integer import IntegerLattice

# Parameters
d = 256 # ring degree
#q = 4293933569		# modulus used for fully-splitting ring with d = 256
#q = 4293925633		# modulus used for fully-splitting ring with d = 128
#q = 167770241	# small modulus q
q = 3489661057 #  large modulus qhat = 2^32-2^29-2^28+2^7+1
ell = 64 # number of irred factors of X^d+1 mod q (r in paper)
tl = 1 # ratio between ell and no. of support coeffs in S_i (power of 2 in [1,ell])
w = 11 # weight of each coeff subset S_i of challenges (hat(w) in paper)
ellp = 64 # total support coeffs in S_i (r in paper)

# Derived parameters
r = d//ell # degree of the irred factors of X^d+1 mod q = d/ell (delta in paper)


# Model results
sd_mu = RDF(1/sqrt(2*ellp)) # s.dev. of mu[j]
sd_S = RDF(sqrt(2*ellp)) # s.dev. of S[j] (compare to Srms in computations)
A = RDF(gamma((w+1)/2)/sqrt(pi()*ellp^w)) # w'th moment of Gaussian with var 1/(2*ellp)
EM1bias = RDF((1-1/q)*A)
EM1 = RDF(1/q + EM1bias) # Expected value of M1
eta = RDF(ellp^w * factorial(ellp-w) / factorial(ellp))
EM2 = eta * EM1 # Expected value of M2
pnoell = EM2^r # Expected value of p (without union bound factor over all irred factors)
p = ell * EM2^r # expected bound for p (with union bound factor over all irred factors)

# logs of Results
lgEM1bias = RDF(log(EM1bias)/log(2))
lgEM1 = RDF(log(EM1)/log(2))
lgeta = RDF(log(eta)/log(2))
lgEM2 = RDF(log(EM2)/log(2))
lgpnoell = RDF(log(pnoell)/log(2))
lgp = RDF(log(p)/log(2))

print ("Input Pars:")
print ("modulus (q in paper) = ", q)
print ("no. of irred factors (r in paper) = ", ell)
print ("no. of support coefficients = ", ellp)
print ("weight of each S_i set (hat(w) in paper)= ", w)

print ("Output Pars:")
print ("S[j] std dev (sd_S) = ",sd_S)
print ("hat(mu)[j] std dev (sd_mu) = ",sd_mu)
print ("log_2 M1bias (lgEM1bias) = ",lgEM1bias)
print ("log_2 M1 (lgEM1) = ",lgEM1)
print ("log_2 eta (lgeta) = ",lgeta)
print ("log_2 M2 (lgEM2) = ",lgEM2)
print ("log_2 pnoell (lgpnoell) = ",lgpnoell)
print ("log_2 p (lgp) = ",lgp)
