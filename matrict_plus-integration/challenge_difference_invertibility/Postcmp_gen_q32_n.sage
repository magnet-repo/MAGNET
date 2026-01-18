# MatRiCT+ invertibility bound computation
# Ron Steinfeld, April 2021
# Part 2: Probability Bound postcomp
# Please cite the MatRiCT+ paper accepted at IEEE S&P 2022 conference (full version available at https://eprint.iacr.org/2021/545) when using scripts

import numpy as np
import time
import csv
from csv import reader
from sage.modules.free_module_integer import IntegerLattice

# Parameters
d = 256 # ring degree
q = 3489661057 # 32-bit modulus
ell = 64 # number of irred factors of X^d+1 mod q (this is the parameter r in the paper)
tl = 1 # ratio between ell and no. of support coeffs in S_i (power of 2 in [1,ell])
w = 11  # weight of the challenges in each coeff subset S_i of challenge = wf/r (def: 16)
p = RDF(0) # prob. of coeff = 0,  prob of +/-1 = (1-p)/2
in_filename = 'STable_gen_q='+str(q)+'_n.csv'	# input filename where the table is stored
out_filename = 'outinvbnd_gen_q='+str(q)+'_n.csv'	# output filename to write the output
read_table = True	# set to 'False' if the table is already read
write_table = False	# set to 'False' if don't want to output to a table

# Derived parameters
r = d//ell # degree of the irred factors of X^d+1 mod q = d/ell (def: 4)
wf = r * w # full weight of coeff set S of challenges
ellp = ell//tl # no. of support coeffs in S_i

lgw = Integer(round(log(RDF(w))/log(2))) #  log_2(w)

if q % (4*ell) != (2*ell+1):	# check if q is of the correct form
	sys.exit("q is not of the correct form!")

# precomputed constants
Ssize = (q-1)//(2*ellp) # no. of elements in S Table
g = mod(primitive_root(q),q) # generator of Z_q^*
z = g^(Ssize) # (2*ellp)'th primitive root of 1 mod q
tprq = (2*RDF.pi()/RDF(q))  # 2*pi/q
q2ellp = (1-p)/(2*RDF(ellp))
tellpq = 2*RDF(ellp)/RDF(q)
# bits of w
wbit =  [Integer(0)] * Integer(5)
wb = w
for i in range(0,5):
	wbit[i] = wb % 2
	wb = wb // 2

# Z: table of <z> group elements (half)
Z = [0] * ell # Z table: Z[k] = z^k mod q
Z[0] = mod(1,q)
for k in range(1,ellp):
	Z[k] = z * Z[k-1]

if read_table:
# Read in from precomp file the S Table: S[j] = 2*sum_{k=0,..,ellp-1} cos(2*pi*(g^j*z^k mod q)/q)
	print ("Reading S Table from file...")
	S = [RDF(0)] * Ssize # S table, entry indexed by j corresponds to representative g^j of Z_q^*/<z>
	with open(in_filename, 'r') as fread:
		# pass the file object to reader() to get the reader object
		csv_reader = csv.reader(fread)
		# Pass reader object to list() to get a list of lists
		list_of_rows = list(csv_reader)
		S =  list_of_rows[0]
		fread.close()
	print ("Finished Reading S Table from file.")

print ("Computing M1 probability bound...")
time0=time.time()
Srms = RDF(0)
bias = RDF(0)
for j in range(0,Ssize):
	Sj = RDF(S[j])
	T = (abs(p + q2ellp * Sj))^w
	bias = bias + T
	Srms = Srms + Sj * Sj
	if (j % 100000) == 0:
		time1=time.time()
		print ("Finished up to iteration ",j/Ssize, ". time = ", time1-time0," seconds")
bias = tellpq  *  bias
lgbias = RDF(log(bias)/log(2))
M1 = RDF(1/q + bias)
lgM1 = RDF(log(M1)/log(2))
explgM1 = RDF(w*log(sqrt(1/(2*RDF(ellp))))/log(2))
Srms = sqrt(Srms/Ssize)
Srms2  = list_of_rows[1]
Srms2 = RDF(Srms2[0])
Srmsexp = sqrt(2*RDF(ellp))
eta = RDF(1)
for i in range(1,w):
	eta = eta * (ellp)/(ellp-i)
M2 = eta * M1
lgM2 = RDF(log(M2)/log(2))
time1=time.time()
# ~ print( "Completed Prob Bound Computation in ",time1-time0," seconds")
# ~ print( "Srms = ",Srms)
# ~ print( "Srms2 = ",Srms2)
# ~ print( "Srmsexp = ",Srmsexp)
# ~ print( "bias = ",bias)
# ~ print( "lgbias = ",lgbias)
# ~ print( "M1 = ", M1)
# ~ print( "lgM1 = ", lgM1)
# ~ print( "explgM1 = ", explgM1)
# ~ print( "r*lgM1 = ", RDF(r)*lgM1)
# ~ print( "M2 = ", M2)
# ~ print( "r*lgM2 = ", RDF(r)*lgM2)
print( "q    = ", q)
print( "d    = ", d)
print( "num of factors= ", ell)
print( "w of one S_i = ", w)
print( "total weight = ", wf)
print( "lgM2 = ", round(lgM2,1))
print( "logp = ", round(RDF(r)*lgM2+log(ell,2),1))
print( "ChSet size = ", round(r*(w+log(binomial(ell,w),2)),1) )
print()

if write_table:
	print( "Writing Results to output file")
	# Write S Table to out file
	f1 = open(out_filename, 'w')
	with f1:
		writefile = csv.writer(f1)
		writefile.writerow(['q',q])
		writefile.writerow(['g',g])
		writefile.writerow(['ell',ell])
		writefile.writerow(['w',w])
		writefile.writerow(['wf',wf])
		writefile.writerow(['Srms',Srms])
		writefile.writerow(['Srms2',Srms2])
		writefile.writerow(['Srmsexp',Srmsexp])
		writefile.writerow(['bias',bias])
		writefile.writerow(['lgbias',lgbias])
		writefile.writerow(['M1',M1])
		writefile.writerow(['lgM1',lgM1])
		writefile.writerow(['explgM1',explgM1])
		writefile.writerow(['r*lgM1',r*lgM1])
		writefile.writerow(['eta',eta])
		writefile.writerow(['M2',M2])
		writefile.writerow(['lgM2',lgM2])
		writefile.writerow(['r*lgM2',r*lgM2])
	f1.close()

	print( "Completed Writing Table to output file")
	


