from __future__ import division

from __future__ import absolute_import
from __future__ import print_function

import jsonpickle

import fileinfo
import namegen
import integralequations
import theorydata
import theory
import comparisons

import numpy as np
import scipy.ndimage

from nlp_euler import NLPEuler
from nlp_fourier import NLPFourier

from extsde import ExtendedSelfDualityEquation
from polynomial_utils import poly, derivcoefs, equivalent_quasimonic_centered, equivalent_quasimonic_centered_cyclic
from functools import partial

from cmath import log,pi,exp

import logging
logger = logging.getLogger(__name__)

EPS = 1e-6

def linearsolve(metric, Pdot, method = "euler", pde_thresh = 5e-11):
    '''Solve the linearized equations around a given harmonic metric,
       for a perturbation of the polynomial P --> P + Pdot.'''
    uvec = metric.u.reshape((metric.grid.nx * metric.grid.ny,) )
    pdotvec = Pdot.reshape((metric.grid.nx * metric.grid.ny,) )

    p = metric.p
    pbar = lambda z: (p(z)).conjugate()

    def flinearreal(p,pbar,u,z,udot,pdot,pdotbar):
        ps = np.abs(p(z))**2
        inhom = - 8.0 * np.exp(-2.0 * u) * (p(z) * pdotbar).real
        return (8.0*udot*(np.exp(2.0*u) + np.exp(-2.0*u)*ps) + inhom)

    def flinearimag(p,pbar,u,z,udot,pdot,pdotbar):
        ps = np.abs(p(z))**2
        inhom = - 8.0 * np.exp(-2.0 * u) * (p(z) * pdotbar).imag
        return (8.0*udot*(np.exp(2.0*u) + np.exp(-2.0*u)*ps) + inhom)

    def dflinear(p,u,z,udot):
        ps = np.abs(p(z))**2
        return (8.0*(np.exp(2.0*u) + np.exp(-2.0*u)*ps))

    freal = lambda udot,z: flinearreal(p,pbar,uvec,z,udot,pdot=pdotvec,pdotbar=pdotvec.conjugate())
    fimag = lambda udot,z: flinearimag(p,pbar,uvec,z,udot,pdot=pdotvec,pdotbar=pdotvec.conjugate())
    d1f = lambda udot,z: dflinear(p,uvec,z,udot)

    udot0 = lambda z: (z * 0.0).real
    v0 = lambda z: (z * 0.0).real

    if method == 'euler':
        NLPSolver = NLPEuler
    elif method == 'fourier':
        NLPSolver = NLPFourier
    else:
        raise ValueError('Unknown solver method "%s"' % method)
        
    Sreal = NLPSolver(freal,d1f,udot0,metric.grid,thresh=pde_thresh,relax=1.0,linear=True)
    Simag = NLPSolver(fimag,d1f,v0,metric.grid,thresh=pde_thresh,relax=1.0,linear=True)
    udot = Sreal.u
    v = Simag.u
    F = udot - 1.0j * v
    return F

def fdcomputeG(pde_nmesh = 300, c = 1.0, rmax = 10.0, method = "euler", pde_maxiter = 8, cutoff = 1.0, pde_thresh = 5e-11, Lambda = 0.60):
    '''Compute hyperkahler metric in K=2 case, P = z^3 - Lambda z - c, by solving PDE.
       We determine the squared norm of a tangent vector perturbing P --> P + eps.'''
    # TODO: make this work beyond A1A2
    coefs = [-c,-Lambda,0,1]
    K = 2

    logger.info("Starting PDE computation of hyperkahler metric element")

    # compute adjusted coefficients
    logger.info("coefs: "+str(coefs))
    coefs = equivalent_quasimonic_centered([coefs])[0]
    logger.info("Monicized coefs: "+str(coefs))

    metric = ExtendedSelfDualityEquation(K=K, coefs=coefs, rmax=rmax, pde_nmesh = pde_nmesh, pde_thresh = pde_thresh, pde_maxiter = pde_maxiter, method = method)

    grid = metric.grid
    u = metric.u
    z = grid.zm
    Pdot = 1.0 + 0.0*z
    P = metric.p(z)

    F = linearsolve(metric = metric, Pdot = Pdot, method = method, pde_thresh = pde_thresh)

    integrand = (4.0*np.exp(-2.0*u) * ( np.abs(Pdot)**2 - (F*P*Pdot.conjugate()) )).real
    boxmeasure = grid.dx * grid.dy
    rcutoff = rmax*cutoff
    mask = np.piecewise(z, [ abs(z) < rcutoff ], [ 1.0, 0.0 ])
    # compute the integral by naive sampling
    innerpart = np.sum(integrand*mask) * boxmeasure

    # compute outer part using semiflat approximation
    # TODO: make this computation more precise
    n = len(coefs) - 1
    if n > 2:
        outerpart = 4 * pi * rcutoff**(-n+2) / (n-2)
    else:
        outerpart = 0

    G = (innerpart + outerpart).real

    integrandsf = 2 * np.abs(Pdot)**2 / np.abs(P)

    logger.info("Computed metric element: %0.8E" % G)

    return {"G": G, "integrand": integrand, "integrandsf": integrandsf}

def ieqcomputeG(c = 1.0, R = 1, eps = EPS, steps = integralequations.IEQ_STEPS, Lambda = 0.60, leadingapprox = False, nterms = 0):

    logger.info("Starting IEQ computation of hyperkahler metric element")

    # TODO: make this work beyond A1A2

    def A1A2name(c, Lambda):
        return "A1A2_cr=%f_ci=%f_Lambda=%f" % (c.real,c.imag,Lambda)

    theoryname0 = A1A2name(c,Lambda)
    theoryname1 = A1A2name(c+eps,Lambda)
    theoryname2 = A1A2name(c+eps*1j,Lambda)

    if leadingapprox:
        cluster0 = integralequations.getApproxCluster(theoryname = theoryname0, R = R, theta = 0, nterms = nterms)
        cluster1 = integralequations.getApproxCluster(theoryname = theoryname1, R = R, theta = 0, nterms = nterms)
        cluster2 = integralequations.getApproxCluster(theoryname = theoryname2, R = R, theta = 0, nterms = nterms)
    else:
        xar0 = integralequations.computeXar(theoryname = theoryname0, R=R, steps = steps)
        xar1 = integralequations.computeXar(theoryname = theoryname1, R=R, steps = steps)
        xar2 = integralequations.computeXar(theoryname = theoryname2, R=R, steps = steps)
        cluster0 = xar0.getCluster()
        cluster1 = xar1.getCluster()
        cluster2 = xar2.getCluster()
    deltax1 = [log(abs(x))-log(abs(x0)) for x,x0 in zip(cluster1,cluster0)]
    deltax2 = [log(abs(x))-log(abs(x0)) for x,x0 in zip(cluster2,cluster0)]

    norm = -(deltax1[0]*deltax2[1] - deltax1[1]*deltax2[0]).real
    tangentnorm = norm / abs(eps)**2

    logger.info("Computed metric element: %0.8E" % tangentnorm)

    return {"cluster0": cluster0, "cluster1": cluster1, "cluster2": cluster2, "deltax1": deltax1, "deltax2": deltax2, "G": tangentnorm, "eps": eps, "theoryname0": theoryname0, "theoryname1": theoryname1, "theoryname2": theoryname2}

def comparemetrics(c = 1, pde_nmesh = 300, eps = EPS, steps = 16384, fdmethod = "euler", pde_thresh = 5e-11, Lambda = 0.60):
    ofd = fdcomputeG(c = c, method = fdmethod, pde_nmesh = pde_nmesh, pde_thresh = pde_thresh, Lambda = Lambda)
    oieq = ieqcomputeG(c = c, eps = eps, steps = steps, Lambda = Lambda)

    fdG = ofd["G"]
    ieqG = oieq["G"]

    return {"ieq": ieqG, "fd": fdG, "delta": ieqG - fdG, "reldiff": 2*(ieqG-fdG) / (ieqG+fdG)}

class fdmetricdata:

    def __init__(self, c, Lambda, pde_nmesh = 300, rmax = 10.0, method = "euler", pde_maxiter = 8, cutoff = 1.0, pde_thresh = 5e-11):
        self.pde_nmesh = pde_nmesh
        self.c = c
        self.Lambda = Lambda
        self.rmax = rmax
        self.method = method
        self.pde_maxiter = pde_maxiter
        self.cutoff = cutoff
        self.pde_thresh = pde_thresh

        data = fdcomputeG(pde_nmesh = pde_nmesh, c = c, Lambda = Lambda, rmax = rmax, method = method, pde_maxiter = pde_maxiter, cutoff = cutoff, pde_thresh = pde_thresh)
        self.G = data["G"]

    def save(self, filename):
        logger.info("Saving fdmetric data to %s" % filename)
        frozen = jsonpickle.encode(self, warn=True)

        with open(filename, "w") as target:
            target.write(frozen)

    @staticmethod
    def load(filename):
        """Load metric data from a file, using jsonpickle.

           arguments: filename (str)
           returns: xarray"""

        logger.info("Loading metric data from %s" % filename)
        with open(filename, "r") as target:
            frozen = target.read()
        metricdata = jsonpickle.decode(frozen)

        return metricdata


def fdComputeAndSaveG(pde_nmesh = 300, c = 1.0, rmax = 10.0, method = "euler", pde_maxiter = 8, cutoff = 1.0, pde_thresh = 5e-11, Lambda = 0.60):
    data = fdmetricdata(pde_nmesh = pde_nmesh, c = c, rmax = rmax, method = method, pde_maxiter = pde_maxiter, cutoff = cutoff, pde_thresh = pde_thresh, Lambda = Lambda)

    filename = namegen.fd_metric_filename(fullpath = True, c = c, Lambda = Lambda)
    data.save(filename)


class ieqmetricdata:

    def __init__(self, c, Lambda, eps = EPS, steps = integralequations.IEQ_STEPS):
        self.c = c
        self.Lambda = Lambda
        self.eps = eps
        self.steps = steps

        data = ieqcomputeG(c = c, Lambda = Lambda, eps = eps, steps = steps)
        self.G = data["G"]

    def save(self, filename):
        logger.info("Saving ieqmetric data to %s" % filename)
        frozen = jsonpickle.encode(self, warn=True)

        with open(filename, "w") as target:
            target.write(frozen)

    @staticmethod
    def load(filename):
        """Load metric data from a file, using jsonpickle.

           arguments: filename (str)
           returns: xarray"""

        logger.info("Loading metric data from %s" % filename)
        with open(filename, "r") as target:
            frozen = target.read()
        metricdata = jsonpickle.decode(frozen)

        return metricdata


def ieqComputeAndSaveG(c, Lambda, eps = EPS, steps = integralequations.IEQ_STEPS):
    data = ieqmetricdata(c = c, Lambda = Lambda, eps = eps, steps = steps)

    filename = namegen.ieq_metric_filename(fullpath = True, c = c, Lambda = Lambda)
    data.save(filename)

def computeAndSaveGs(c, Lambda, eps = EPS, steps = integralequations.IEQ_STEPS, rmax = 10.0, pde_nmesh = 1000, fdmethod = "euler", pde_maxiter = 8, cutoff = 1.0, pde_thresh = 5e-11):
    ieqComputeAndSaveG(c = c, Lambda = Lambda, eps = eps, steps = steps)
    fdComputeAndSaveG(c = c, Lambda = Lambda, rmax = rmax, pde_nmesh = pde_nmesh, method = fdmethod, pde_maxiter = pde_maxiter, cutoff = cutoff, pde_thresh = pde_thresh)
