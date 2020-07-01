from __future__ import absolute_import
import theorydata
import scipy
import polynomial_utils

PERIOD_TOL = 1e-15
PERIOD_NSTEPS = 1e9

import logging
logger = logging.getLogger(__name__)

def discriminant(K, Pcoefs, z):
    if K>=2:
        P2coefs = Pcoefs[0]
        P2 = polynomial_utils.poly(P2coefs, z)
    if K>=3:
        P3coefs = Pcoefs[1]
        P3 = polynomial_utils.poly(P3coefs, z)
    if K == 2:
        return -4*P2
    if K == 3:
        return -4*(P2**3)-27*(P3**2)
    raise NotImplementedException("Discriminant only implemented for K=2,3")

def F(K, Pcoefs, z, x):
    return x**K + sum( [x**(K-n) * polynomial_utils.poly(Pcoefs[n-2], z) for n in range(2,K+1)] )

def dFdx(K, Pcoefs, z, x):
    return K*(x**(K-1)) + sum( [(K-n)*x**(K-n-1) * polynomial_utils.poly(Pcoefs[n-2], z) for n in range(2,K+1)] )

def dFdz(K, Pcoefs, z, x):
    return sum( [x**(K-n) * polynomial_utils.poly(polynomial_utils.derivcoefs(Pcoefs[n-2]), z) for n in range(2,K+1)] )

def dxdz(K, Pcoefs, z, x):
    return -dFdz(K, Pcoefs, z, x)/dFdx(K, Pcoefs, z, x)

def turningpoint(K, Pcoefs, zinitial):
    return findroot(lambda z: discriminant(K, Pcoefs, z), zinitial)

def turningpoints(K, Pcoefs, zinitial):
    return [turningpoint(K, Pcoefs, zi) for zi in zinitial]

def spectralpoint(K, Pcoefs, z, xinitial):
    return findroot(lambda x: F(K, Pcoefs, z, x), xinitial)

def spectralpoints(K, Pcoefs, z, xinitial):
    return [spectralpoint(K, Pcoefs, z, xi) for xi in xinitial]

def findroot(F, zinitial):
    def Fcpl(data):
        x,y = data
        f = F(x+y*1j)
        return [f.real, f.imag]
    opt = scipy.optimize.root(Fcpl, [zinitial.real, zinitial.imag])
    z = opt["x"][0] + 1j*opt["x"][1]
    if abs(F(z)) > 1e-8:
        raise ArithmeticError("Root finder failed")
    return z

def evolveroots(K, Pcoefs, zinitial, zfinal, xinitial):
    # evolve sheets of spectral curve along a straight-line path

    def z(t):
        return t * zfinal + (1-t) * zinitial

    zprime = zfinal-zinitial

    def odef(t, y):
        return [zprime*dxdz(K, Pcoefs, z(t), x) for x in y]

    solver = scipy.integrate.complex_ode(odef)
    solver.set_integrator("dopri5",rtol=PERIOD_TOL,nsteps=PERIOD_NSTEPS)
    solver.set_initial_value(xinitial,0.0)
    return solver.integrate(1.0)

def lineintegral(K, Pcoefs, zinitial, zfinal, xinitial):
    # solve ODE for two quantities at once: the period and the root

    def z(t):
        return t * zfinal + (1-t) * zinitial

    zprime = zfinal-zinitial

    def odef(t, y):
        x = y[1]
        return [zprime*x, zprime*dxdz(K, Pcoefs, z(t), x)]

    solver = scipy.integrate.complex_ode(odef)
    solver.set_integrator("dopri5",rtol=PERIOD_TOL,nsteps=PERIOD_NSTEPS)
    solver.set_initial_value([0,xinitial],0.0)
    return solver.integrate(1.0)[0]

def allintegrals(K, Pcoefs, zinitial, zfinal, zbase, xbase):
    zhalfway = 0.5*zinitial + 0.5*zfinal
    xlist = evolveroots(K, Pcoefs, zbase, zhalfway, xbase)
    return [lineintegral(K, Pcoefs, zhalfway, zfinal, x)-lineintegral(K, Pcoefs, zhalfway, zinitial, x) for x in xlist]
