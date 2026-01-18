# MatRiCT+ invertibility bound computation
# Ron Steinfeld, April 2021
# Part 1: S table precomp
# Please cite the MatRiCT+ paper accepted at IEEE S&P 2022 conference (full version available at https://eprint.iacr.org/2021/545) when using scripts

import numpy as np
import time
import csv
from sage.modules.free_module_integer import IntegerLattice

# Parameters
d = 256 # ring degree
q = 3489661057 # 32-bit modulus
ell = 64 # number of irred factors of X^d+1 mod q (this is the parameter r in the paper)
tl = 1 # ratio between ell and no. of support coeffs in S_i (power of 2 in [1,ell])
out_filename = 'STable_gen_q='+str(q)+'_n.csv'	# name of the file output is written

# Derived parameters
r = d//ell # degree of the irred factors of X^d+1 mod q = d/ell (def: 4)
ellp = ell//tl # no. of support coeffs in S_i

if q % (4*ell) != (2*ell+1):	# check if q is of the correct form
	sys.exit("q is not of the correct form!")

# precomputed constants
Ssize = (q-1)//(2*ellp) # no. of elements in S Table
g = mod(primitive_root(q),q) # generator of Z_q^*
z = g^(Ssize) # (2*ellp)'th primitive root of 1 mod q
tprq = (2*RDF.pi()/RDF(q))  # 2*pi/q

# Z: table of <z> group elements (half)
Z = [0] * ellp # Z table: Z[k] = z^k mod q
Z[0] = mod(1,q)
for k in range(1,ellp):
	Z[k] = z * Z[k-1]

# S Table: S[j] = 2*sum_{k=0,..,ellp-1} cos(2*pi*(g^j*z^k mod q)/q)
time0=time.time()
S = [RDF(0)] * Ssize # S table, entry indexed by j corresponds to representative g^j of Z_q^*/<z>
Crep = mod(1,q)  # current coset representative g^j of Z_q^*/<z>
Srms =  0
#I = g;
#J = RDF(I);
for j in range(0,Ssize):
	S[j] = 0
	for k in range(0,ellp):
		I = RDF(Crep * Z[k]) # I = g^j * z^k mod q = k'th element in coset containing g^j
	  #J = RDF(I)
		T = tprq * I # cos radian angle 2*pi*(g^j*z^k mod q)/q
		# print ("I = ", I)
		# print ("T = ", T)
		S[j] = S[j] + cos(T)
	S[j] = 2 * S[j]
	Srms = Srms + S[j]*S[j]
	Crep = g * Crep
	if (j % 100000) == 0:
		time1=time.time()
		print ("Finished up to iteration ",j/Ssize, ". time = ", time1-time0," seconds")
time1=time.time()
Srms = sqrt(Srms/Ssize)
print ("Completed S Table in ",time1-time0," seconds")
print ("Srms = ",Srms)

print ("Writing Table to output file...")
# Write S Table to out file
f1 = open(out_filename, 'w')
with f1:
	writefile = csv.writer(f1)
	writefile.writerow(S)
	writefile.writerow([Srms])
	writefile.writerow([q])
	writefile.writerow([g])
	writefile.writerow([ell])
f1.close()

print ("Completed Writing Table to output file")
