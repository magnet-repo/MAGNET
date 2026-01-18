from math import *
from model_BKZ import *
from proba_util import gaussian_center_weight

log_infinity = 9999

STEPS_b = 5
STEPS_m = 5
MIN_b = 300


class MSISParameterSet:
    def __init__(self, n, w, h, B1, B2, B3, B4, B5, m1, m2, m3, m4, m5, q, norm=""):
        self.n = n      # Ring Dimension
        self.w = w      # MSIS dimension
        self.h = h      # Number of equations
        # ME_COMMENT: Added below params B1, B2, ... and m1, m2, ...
        self.B1 = B1    # norm bound 1
        self.B2 = B2    # norm bound 2
        self.B3 = B3    # norm bound 3
        self.B4 = B4    # norm bound 2
        self.B5 = B5    # norm bound 3
        self.m1 = m1    # number of ring elements with norm bound 1
        self.m2 = m2    # number of ring elements with norm bound 1
        self.m3 = m3    # number of ring elements with norm bound 1
        self.m4 = m4    # number of ring elements with norm bound 1
        self.m5 = m5    # number of ring elements with norm bound 1
        self.q = q      # Modulus
        self.norm = norm


def SIS_l2_cost(q, w, h, B, b, cost_svp=svp_classical, verbose=False):
    """ Return the cost of finding a vector shorter than B with BKZ-b if it works.
    The equation is Ax = 0 mod q, where A has h rows, and w collumns (h equations in dim w).
    """
    if B>=q:
        if verbose:
            print ("Do not know how to handle B>=q in l2 norm. Concluding 0 bits of security.")
        return 0
    l = BKZ_first_length(q, h, w-h, b)
    if l > B:
        return log_infinity 
    if verbose:
        print ("Attack uses block-size %d and %d equations"%(b, h))
        print ("shortest vector used has length l=%.2f, q=%d, `l<q'= %d"%(l, q, l<q))
    return cost_svp(b)


def SIS_linf_cost(q, w, h, B1, m1, B2, m2, B3, m3, B4, m4, B5, m5, b, cost_svp=svp_classical, verbose=False, attack_variant=2):
    """ Return the cost of finding a vector with m1 first coords shorter than B1,
      "  m2 following coords shorter than B2 <= B1, m3 coords shorter than B3 <= B2,
      "  m4 coords shorter than B4 <= B3, m5 = w-(m1+m2+m3+m4) coords in infinity norm,
      "  using BKZ-b, if it works.
    The equation is Ax = 0 mod q, where A has h rows, and w columns (h equations in dim w).
    w = m1 + m2 + m3 + m4 + m5.
    """
    
    if attack_variant == 2:
        c12 = B1/B2
        c13 = B1/B3
        c14 = B1/B4
        c15 = B1/B5
        (i, j, L) = construct_BKZ_shape_randomized(q, h, w-h, b, m1, m2, m3, m4, m5, c12, c13, c14, c15)
    elif attack_variant == 1:
        c12 = 1
        c13 = 1
        c14 = 1
        c15 = 1
        (i, j, L) = construct_BKZ_shape_randomized(q, h, w-h, b, m1, m2, m3, m4, m5, c12, c13, c14, c15)
    elif attack_variant == 0:
        c12 = 1
        c13 = 1
        c14 = 1
        c15 = 1
        (i, j, L) = construct_BKZ_shape(q, h, w-h, b)	
        # ~ (i, j, L) = construct_BKZ_shape( h*log(q), h, w-h, b)	
    else:
        raise ValueError("Incorrect attack variant!")

    # ~ (i, j, L) = construct_BKZ_shape(h * log(q), h, w-h, b)
    
    
    # ***************************
    
    l = exp(L[i])
    d = j - i + 1 #Is this a bug in Dilithium script code? Should be j-i as only coords i,...j-1 (j-i coords) are in the sloping Li sequence?
    sigma = l / sqrt(j - i + 1)
    p_middle1 = gaussian_center_weight(sigma, B1)  	  # = Pr[z < B1] 
    p_middle2 = gaussian_center_weight(sigma, c12*B2) # = Pr[z*c12 < c12*B2] since we scaled the next m2 coords z of SIS vector by c12
    p_middle3 = gaussian_center_weight(sigma, c13*B3) # = Pr[z*c13 < c13*B3] since we scaled the next m3 coords z of SIS vector by c13
    p_middle4 = gaussian_center_weight(sigma, c14*B4) # = Pr[z*c14 < c14*B4] since we scaled the next m4 coords z of SIS vector by c14
    p_middle5 = gaussian_center_weight(sigma, c15*B5) # = Pr[z*c15 < c15*B5] since we scaled the next m5 coords z of SIS vector by c15
    
    p_head1 = min(1, 2.*B1 / q)
    p_head2 = min(1, 2.*B2 / q)
    p_head3 = min(1, 2.*B3 / q)
    p_head4 = min(1, 2.*B4 / q)
    p_head5 = min(1, 2.*B5 / q)

    dhead1 = min(m1, i)
    dhead2 = max(0, i-dhead1)
    dhead2 = min(dhead2, m2)
    dhead3 = max(0, i-dhead1-dhead2)
    dhead3 = min(dhead3, m3)
    dhead4 = max(0, i-dhead1-dhead2-dhead3)
    dhead4 = min(dhead4, m4)
    dhead5 = max(0, i-dhead1-dhead2-dhead3-dhead4)
    dhead5 = min(dhead5, m5)

    d1 = min(m1 - dhead1, d)
    d2 = max(0, d-d1)
    d2 = min(d2, m2 - dhead2)
    d3 = max(0, d-d1-d2)
    d3 = min(d3, m3 - dhead3)
    d4 = max(0, d-d1-d2-d3)
    d4 = min(d4, m4 - dhead4)
    d5 = max(0, d-d1-d2-d3-d4)
    d5 = min(d5, m5 - dhead5)

    log2_eps = d1 * log(p_middle1, 2) + d2 * log(p_middle2, 2) + d3 * log(p_middle3, 2) + d4 * log(p_middle4, 2) + d5 * log(p_middle5, 2) + dhead1 * log(p_head1, 2) + dhead2 * log(p_head2, 2) + dhead3 * log(p_head3, 2) + dhead4 * log(p_head4, 2) + dhead5 * log(p_head5, 2)
    log2_R = max(0, - log2_eps - nvec_sieve(b))
    
    # ***************************

    L2 = []
    for a in L:
        L2.append(a/log(2))

    if verbose:
        print ("Attack uses block-size %d and %d dimensions, with %d q-vectors"%(b, w, i))
        # ~ print ("L =",L2)
        print ("d =", d)
        print ("(d1, d2, d3, d4, d5) =",(d1, d2, d3, d4, d5))
        print ("(m1, m2, m3, m4, m5) =",(m1, m2, m3, m4, m5))
        print ("(dhead1, dhead2, dhead3, dhead4, dhead5) =",(dhead1, dhead2, dhead3, dhead4, dhead5))
        print ("log2(epsilon) = %.2f, log2 nvector per run %.2f"%(log2_eps, nvec_sieve(b)))
        print ("shortest vector used has length l=%.2f, q=%d, `l<q'= %d"%(l, q, l<q))
    return cost_svp(b) + log2_R


def SIS_optimize_attack(q, max_w, h, B1, m1, B2, m2, B3, m3, B4, m4, B5, m5, cost_attack=SIS_linf_cost, cost_svp=svp_classical, verbose=False, attack_variant=2):
    """ Find optimal parameters for a given attack
    """
    best_cost = log_infinity
    best_w = 0
    best_b = 0
    cost = 0
    sanity_check = 0

    for b in range(MIN_b, max_w, STEPS_b):
        if cost_svp(b) > best_cost:
            break
        sanity_check += 1
        for w in range(max(h+1, b+1), max_w, STEPS_m):    # Do need to exhaust w here attack is suboptimal with w=maxm  
            # ~ print("w =", w)
            cost = cost_attack(q, w, h, B1, m1, B2, m2, B3, m3, B4, m4, B5, m5, b, cost_svp=cost_svp, attack_variant=attack_variant)
            if cost<best_cost:
                best_cost = cost
                best_w = w
                best_b = b

    if sanity_check < 2:
        raise ValueError("MIN_b is too big! Choose smaller MIN_b!")

    if verbose:
        cost_attack(q, best_w, h, B1, m1, B2, m2, B3, m3, B4, m4, B5, m5, best_b, cost_svp=cost_svp, verbose=verbose, attack_variant=attack_variant)

    return (best_w, best_b, best_cost)


def check_eq(m_pc, m_pq, m_pp):
    if (m_pc != m_pq):
        print("m and b not equals among the three models")
    if (m_pq != m_pp):
        print("m and b not equals among the three models")


def MSIS_summarize_attacks(ps, attack_variant=2):
    """ Create a report on the best primal and dual BKZ attacks on an l_oo - MSIS instance
    """
    q = ps.q
    h = ps.n * ps.h
    max_w = ps.n * ps.w
    # ME_COMMENT: Added below setting of B1, B2, ... and m1, m2, ...
    B1 = ps.B1
    B2 = ps.B2
    B3 = ps.B3
    B4 = ps.B4
    B5 = ps.B5
    m1 = ps.n * ps.m1	# Multiply ps.m1 by ps.n because m1 polynomials lead to ps.n * ps.m1 coefficients over Z
    m2 = ps.n * ps.m2
    m3 = ps.n * ps.m3
    m4 = ps.n * ps.m4
    m5 = ps.n * ps.m5


    if ps.norm == "linf":
        attack = SIS_linf_cost
    elif ps.norm == "l2":
        attack = SIS_l2_cost
    else:
        raise ValueError("Unknown norm: "+ps.norm)

    (m_pc, b_pc, c_pc) = (1,1,1)	# the classical estimation is removed to save time, uncomment below to run it
    # ~ (m_pc, b_pc, c_pc) = SIS_optimize_attack(q, max_w, h, B1, m1, B2, m2, B3, m3, B4, m4, B5, m5, cost_attack=attack, cost_svp=svp_classical, verbose=True)
    (m_pq, b_pq, c_pq) = SIS_optimize_attack(q, max_w, h, B1, m1, B2, m2, B3, m3, B4, m4, B5, m5, cost_attack=attack, cost_svp=svp_quantum, verbose=True, attack_variant=attack_variant)
    (m_pp, b_pp, c_pp) = (1,1,1)	# the best plausible estimation is removed to save time, uncomment below to run it
    # ~ (m_pp, b_pp, c_pp) = SIS_optimize_attack(q, max_w, h, B1, m1, B2, m2, B3, m3, B4, m4, B5, m5, cost_attack=attack, cost_svp=svp_plausible, verbose=False)
    

    # ~ check_eq(m_pc, m_pq, m_pp)
    # ~ check_eq(b_pc, b_pq, b_pp)

    print("SIS: dim=%d & blocksize=%d & cost_pq=%d"%(m_pq, b_pq, int(floor(c_pq))))

    # ~ return (b_pq, int(floor(c_pc)), int(floor(c_pq)), int(floor(c_pp)))
    return (m_pq, b_pq, c_pq)
