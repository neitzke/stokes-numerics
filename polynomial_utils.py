"""Utilities for working with polynomial differentials"""

# Polynomial is always represented by its coefficient vector, starting with constant term

from __future__ import absolute_import
import numpy as np
import cmath
from math import factorial
from scipy.integrate import quad
from scipy.special import binom

def cpow(z,a):
    return cmath.exp(a*cmath.log(z))

_third = 1.0/3.0

def ccbrt(z):
    if z == 0:
        return 0
    else:
        return cpow(z,_third)

def poly(coefs,z):
    """Evaluate a polynomial at z given its coefficient vector"""
    result = 0j
    for c in reversed(coefs):
        result = result*z + c
    return result

def derivcoefs(coefs):
    primecoefs = []
    for k in range(len(coefs)-1):
        primecoefs.append((k+1)*coefs[k+1])
    return primecoefs

def pos_pow(z,k):
    """z^k if k is positive, else 0.0"""
    if k>=0:
        return z**k
    else:
        return 0.0

def translate_poly(coefs,a):
    """Given coefficients of p(z), return those of p(z+a)"""
    N = len(coefs)
    return [ sum(coefs[k]*binom(k,j)*pos_pow(a,k-j) for k in range(N)) for j in range(N) ]

def make_centered(coefs):
    """Given list of polynomial differentials where the top one is a MONIC polynomial, 
    return a list of differentials where the top one is centered, equivalent by translation"""

    topcoefs = coefs[-1]

    a = topcoefs[-1]
    b = topcoefs[-2]
    translationdistance = -b / (a*float(len(topcoefs)-1))
    return [translate_poly(item, translationdistance) for item in coefs]

def equivalent_quasimonic_cyclic(K,coef):
    """Transform a single K-differential to one where abs(leading coef)==1"""
    return equivalent_quasimonic([ [] ]*(K-1) + [coef])[-1]

def equivalent_quasimonic_centered_cyclic(K,coef):
    """Transform a single K-differential to one where abs(leading coef)==1 and
    root sum is zero"""
    return equivalent_quasimonic_centered([ [] ]*(K-1) + [coef])[-1]

def equivalent_quasimonic(coefs):
    """Transform list of polynomial differentials to one where the
       highest-order differential will have top coef of unit magnitude"""

    # first we have to find the rescaling factor
    # top order differential will be a K-differential
    topK = len(coefs)+1
    topdiffcoefs = coefs[-1]
    topdegree = len(topdiffcoefs)-1
    L = np.log(abs(topdiffcoefs[-1])) / float(topdegree+topK)

    # now transform all the differentials
    newcoefs = []
    for n in range(len(coefs)):
        K = n+2
        curcoefs = coefs[n]
        newcoefs.append( [ np.exp(-(k+K)*L)*a for k,a in enumerate(curcoefs) ] )

    return newcoefs

def equivalent_quasimonic_centered(coefs):
    """Transform list of polynomial differentials to an equivalent list where 
       the highest-order differential will be quasi-monic and centered"""
    return make_centered(equivalent_quasimonic(coefs))

def fujiwara_root_bound(coefs):
    """Return radius of a disk guaranteed to contain all of the roots of a polynomial"""
    an = coefs[-1]
    X = [ abs(a/an)**(1.0 / (k+1)) for k,a in enumerate(reversed(coefs[:-1])) ]
    return 2.0*max( X[:-1] + [ 0.5**(1.0/(len(coefs)-1))*X[-1] ] )

def quadratic_roots(coefs):
    c,b,a = [complex(x) for x in coefs]
    D = (b*b - 4*a*c)
    x1 = 0.5*(-b + np.sqrt(D))/a
    x2 = 0.5*(-b - np.sqrt(D))/a
    return x1,x2

def cubic_roots(coefs):
    """Cardano's formula; based on http://stackoverflow.com/questions/39474254/"""
    d,c,b,a = [complex(x) for x in coefs]
    ihs3 = 0.8660254037844386j 
    Q = (3*a*c - b**2)/ (9*a**2)
    R = (9*a*b*c - 27*a**2*d - 2*b**3) / (54 * a**3)
    D = Q**3 + R**2
    E = -b / (3*a)

    S = ccbrt(R + np.sqrt(D))
    T = ccbrt(R - np.sqrt(D))

    x1 = S + T + E
    x2 = -0.5*(S+T) + E + ihs3*(S-T)
    x3 = -0.5*(S+T) + E - ihs3*(S-T)
    return x1,x2,x3

def root_enclosing_radius(coefs):
    """Radius of disk centered at 0 guaranteed to contain the roots; optimal for degree 2 and 3"""
    if len(coefs)==3:
        return max(abs(x) for x in quadratic_roots(coefs))
    elif len(coefs)==4:
        return max(abs(x) for x in cubic_roots(coefs))
    else:
        return fujiwara_root_bound(coefs)
    
def radial_seg_len(K,coefs,theta,r0,r1):
    """Length of radial segment with respect to abs(poly)^(1/K) metric"""
    f = lambda t:pow(abs(poly(coefs,np.exp(1j*theta)*t)),1.0/K)
    return quad(f,r0,r1)[0]

def min_annular_sep(K,coefs,r0,r1,theta_samples=100):
    """Approximate distance between |z|=r0 and |z|=r1 in metric abs(poly)^(1/3)"""
    return min(radial_seg_len(K,coefs,t,r0,r1) for t in np.arange(0,2*np.pi,2*np.pi/theta_samples))

def pde_error_model(K,coefs,r,nmesh,exp_coef=1.0,quad_coef=1.0,r0=None):
    """Estimate error in harmonic metric as sum of exponential contribution (from boundary values) and quadratic contribution (from mesh size)"""
    if r0 != None:
        rmin = r0
    else:
        rmin = root_enclosing_radius(coefs)
    d = min_annular_sep(K,coefs,rmin,r)
    if K==2:
        d1f_zero = 16.0
    else:
        d1f_zero = 12.0
    exp_term=np.exp(-np.sqrt(d1f_zero)*d)/np.sqrt(d)
    quad_term=(float(r)/nmesh)**2
    return exp_coef*exp_term + quad_coef*quad_term

def approx_best_rmax(K,coefs,nmesh,exp_coef=1.0,quad_coef=1.0,min_sep=0.5,max_sep=100.0,step=0.1,subdiv=20):
    """Approximate the best RMAX value based on the pde_error_model"""
    r0 = root_enclosing_radius(coefs)
    eps = lambda r:pde_error_model(K,coefs,r,nmesh,exp_coef=exp_coef,quad_coef=quad_coef,r0=r0)
    step = float(step)
    
    r = r0 + min_sep
    err = eps(r)
    while r < max_sep:
        if eps(r + step) < err:
            r = r + step
            err = eps(r)
        else:
            break
    ir0 = max(r-step,min_sep)
    ir1 = min(r+step,max_sep)
    a = np.arange(ir0,ir1,step/subdiv)
    epsvec = [ eps(r) for r in a ]
    mi = np.argmin(epsvec)
    return a[mi]

def partitions(n, I=1):
    yield (n,)
    for i in range(I, n//2 + 1):
        for p in partitions(n-i, i):
            yield (i,) + p

def expand_root(K,coefs):
    """Get the expansion of (-P(z))^(1/K) around large z,
       up to and including order 1/z, assuming P(z) monic
       function returns ([a_0, a_1, ..., a_m], clead) 
       where the expansion is (-P(z))^(1/K) ~ (clead)**(1/K) * (z)**(n/K) * (sum_k a_k z^(-k))
       and m = floor(n/K + 1), clead = leading coeff of -P(z)
    """
    out = []
    n = len(coefs) - 1
    m = n//K + 1
    clead = -coefs[-1]
    # we divide -P by clead under the radical
    rescaledcoefs = [-c/clead for c in coefs]
    for k in range(m+1):
        # work on the coefficient of z^(-k)
        if k == 0:
            expcoef = 1
        else:
            expcoef = 0
            # each partition of k gives a contribution to z^(-k)
            for p in partitions(k):
                # overall factor depending on number of pieces in partition
                curterm = binom(1.0/K, len(p))*factorial(len(p))
                for a in range(k+1):
                    # now extract the part coming from z^(-a) terms
                    # if p contains c copies of z^(-a) then get a factor (rescaledcoefs[-a-1])**c / c! z^(-ac)
                    c = p.count(a)
                    curterm *= rescaledcoefs[-a-1]**c / float(factorial(c))
                expcoef += curterm
        out.append(expcoef)
    return out,clead
