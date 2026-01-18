These scripts are modifications of Dilithium's scripts available at https://github.com/pq-crystals/security-estimates.

To run the scripts, modify "ASIS_est.sage" with the choice of parameters. The scripts estimate the security of finding a solution x with A*x = 0 mod q where,
- x is of dimension m1+m2+m3+m4+m5 over R_q (not Z_q),
- m1 polynomials in x have infinity norm bounded by B1,
- m2 polynomials in x have infinity norm bounded by B2,
- m3 polynomials in x have infinity norm bounded by B3,
- m4 polynomials in x have infinity norm bounded by B4,
- m5 polynomials in x have infinity norm bounded by B5,
- B1 >= B2 >= B3 >= B4 >= B5.

As a result, the scripts allow the setting of (up to) 5 different infinity norm bounds.

IMPORTANT NOTE: These scripts are run in SageMath version 9.0 using Python 3.7.3. As the syntax of the 'print' function in Python has been changed, the scripts may not be compatible with older versions of Python. If an error due to some 'print' call is returned, try commenting out the print commands and then manually check the values of 'bkz_block' and 'round(delta_BKZ(bkz_block),5)'. 