# These scripts implement the security estimation of "Asymmetric MSIS" from https://eprint.iacr.org/2019/510.pdf 
# by extending Dilithium's scripts at https://github.com/pq-crystals/security-estimates
# Muhammed F. Esgin and Ron Steinfeld
# Please cite the MatRiCT+ paper accepted at IEEE S&P 2022 conference (full version available at https://eprint.iacr.org/2021/545) when using scripts


load("MSIS_security.py")


# the below parameters are specific to MatRiCT+
N = 11		# Ring size of the ring signature
lwerank = 5
sisrank = 5

# set the modulus and polynomial ring dimension
q = 3489661057	# 31.7 bit q
d = 256		# dim(R_q), i.e., R_q = Z_q[X]/(X^d+1)

# Set the infinity-norm bounds
B1 = ceil(2^30.3)
B2 = ceil(2^27.56)
B3 = ceil(2^22.64)
B4 = ceil(2^19.92)
B5 = ceil(2^16.6)

# ~ B2 = B3 = B4 = B5 = B1		# uncomment to run the scripts for standard MSIS problem

# Set the dimensions (over R_q) of MSIS solution parts
m1 = 1
m2 = N+1
m3 = sisrank + lwerank
m4 = 1
m5 = N+1


# ~ ********************** Results ********************************
print("sisrank  =", sisrank)
print("log2(q)  =", round(log(q,2),2))
print("log2(B1) =", round(log(B1,2),2))
print("log2(B2) =", round(log(B2,2),2))
print("log2(B3) =", round(log(B3,2),2))
print("log2(B4) =", round(log(B4,2),2))
print("log2(B5) =", round(log(B5,2),2))
print()


attack_variant = 0			# The attack variant should be in {0,1,2}. '0' is the "traditional attack". The for-loop below runs over all attacks
for attack_variant in [0,1,2]:
	# function below assumes B1>=B2>=B3>=B4>=B5
	print("***** Attack Variant", attack_variant, " *********")
	params = MSISParameterSet(d, m1+m2+m3+m4+m5, sisrank, B1, B2, B3, B4, B5, m1, m2, m3, m4, m5, q, norm="linf")
	(m_pq, b_pq, c_pq) = MSIS_summarize_attacks(params, attack_variant=attack_variant)
	print("BKZ block size =", b_pq)
	print("RHF =", round(delta_BKZ(b_pq),5))
	print("Cost=", round(c_pq,5))
	print()
