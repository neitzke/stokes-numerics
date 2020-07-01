"""Compute vertices of polygons (ideal in H^2 or convex in RP^2)
corresponding to real cyclic Higgs bundles of rank 2 or 3 and save to
files for later comparison with integral equation results.

""" 

from __future__ import division
from __future__ import absolute_import
from __future__ import print_function

import numpy as np
import os
from numpy import log, sqrt
import jsonpickle
import jsonpickle.ext.numpy as jsonpickle_numpy
jsonpickle_numpy.register_handlers()
import random
import copy
import scipy
import scipy.integrate
import warnings
import time
import theorydata

from planarode import planarODE
from namegen import fd_filename
from squaregrid import SquareGrid
from extsde import ExtendedSelfDualityEquation
from functools import partial
from polynomial_utils import poly, equivalent_quasimonic_centered, root_enclosing_radius, approx_best_rmax, pde_error_model, expand_root
from cmath import pi, sin, cos, exp, phase
from codeversion import codeversion

import logging
logger = logging.getLogger(__name__)

import fileinfo
import copy

#----------------------------------------------------------------------
# Default values for key parameters in the framedata calculations
PDE_NMESH=511        # PDE grid size
PDE_THRESH=1e-9
PDE_MAXITER=5000
ODE_NSTEPS=500000    # A maximum; should be significantly larger than 1/ODE_RSTEP
ODE_RSTEP=1e-4       # Initial step size may be adjusted by the ODE integrator
ODE_THRESH=1e-14
ODE_RMAXFACTOR=1.0
PDE_NMESH_LIST=[511,1023,2047]  # Used in comparisons.compareClusterRichardson
#----------------------------------------------------------------------


def numint(f, r0):
    def g(t):
        return (f(r0*exp(t))*r0*exp(t)).real
    # TODO: fix these constants in more scientific fashion
    eps = 1e-3
    L = 12.0
    samples = [g(n*eps) for n in range(int(L/eps))]
    return sum(samples)*eps

class framedata:
# members:
#   K
#   Pcoefs
#   coefs
#   degree
#   R
#   theta
#   params: dictionary storing pde_nmesh, rmax, ode_thresh, ode_rstep, ode_nsteps
#   version
#   vertexlist
#   framelist
#   metric (can be None)
#   l2norm
#   unregl2norm
#   regularizationterm

# computes the polygon corresponding to polynomial R^K exp(-K i theta) P(z)
# where P(z) has coefficient vector coefs (constant term first), e.g.
#   coefs = [1,0,2] for 2z^2 + 1
#   coefs = [-2,0,-3,1] for z^3 - 3z^2 - 2

    def __init__(self, K = None, Pcoefs = None, R = 1.0, theta = 0.0, make_monic = True, pde_thresh = PDE_THRESH, pde_nmesh = PDE_NMESH, ode_rstep = ODE_RSTEP, ode_nsteps = ODE_NSTEPS, ode_thresh = ODE_THRESH, ode_rmaxfactor = ODE_RMAXFACTOR, theoryname = None, storemetric = False, rmax = None, rmaxfactor = 1.0, method = None, pde_maxiter = PDE_MAXITER, oper = False, theory = None, absh = 1.0):
        # self.K = K

        # TODO: change everything over to use theory object when it's given, get
        # rid of unused arguments above
        if theory is None and theoryname is not None:
            theory = theorydata.getdata(theoryname)
        self.theory = theory

        # TODO: Investigate cleaner solution than deepcopy in setting Pcoefs, coefs0
        # It is currently used to avoid jsonpickle.encode(make_refs=False)
        # writing strings instead of lists; alternative
        # jsonpickle.encode(make_refs=True) makes non-human-readable JSON and
        # seems potentially fragile for long-term serialization
        Pcoefs = copy.deepcopy(theory["Pcoefs"])
        self.K = theory["K"]
        K = self.K
        self.Pcoefs = Pcoefs

        # these two lines use just the K-differential part of Pcoefs
        self.coefs0 = copy.deepcopy(Pcoefs[K-2]) # As given, before scaling and rotation
        self.degree = len(self.coefs0) - 1

        self.theta = theta
        self.version = codeversion
        self.theoryname = theoryname
        self.dataversion = 13
        self.method = method
        self.make_monic = make_monic
        self.oper = oper
        self.real = not oper

        # start the clock
        start_time = time.time()

        if K == 2 or (K == 3 and Pcoefs[0]==[0]):
            self.onlyonedifferential = True
        else:
            self.onlyonedifferential = False

        logger.info("Raw coefs: "+str(self.Pcoefs))

        if oper:
            self.Pcoefsadjusted = [[(a/absh**n)*exp(-1j*n*theta) for a in Pcoefs[n-2]] for n in range(2,len(Pcoefs)+2)]
            self.absh = absh
        else:
            self.Pcoefsadjusted = [[a*(R**n)*exp(-1j*n*theta) for a in Pcoefs[n-2]] for n in range(2,len(Pcoefs)+2)]
            self.R = R
        logger.info("Rescaled coefs: "+str(self.Pcoefsadjusted))

        if make_monic:
            self.Pcoefsadjusted = equivalent_quasimonic_centered(self.Pcoefsadjusted)
            logger.info("Monicized coefs: "+str(self.Pcoefsadjusted))

        # decide which rmax to use
        if not oper:
            if rmax is None:
                if self.onlyonedifferential:
                    rmax = rmaxfactor*approx_best_rmax(self.K, self.Pcoefs[K-2], pde_nmesh)
                else:
                    rmax = 8.0
            logger.info("Using rmax: %g" % rmax)
            if self.onlyonedifferential:
                logger.info("Error model prediction: %g" % pde_error_model(self.K, self.Pcoefs[K-2], rmax, pde_nmesh))
        if oper:
            if rmax is None:
                # TODO: make a more scientific choice of rmax here
                if self.onlyonedifferential:
                    r = root_enclosing_radius(self.Pcoefs[K-2])
                    logger.info("Root enclosing radius: %g" % r)
                    rmax = max(r*8.0,8.0)
                else:
                    rmax = 8.0
            logger.info("Using rmax: %g" % rmax)

        self.params = {"rmax": rmax, "ode_rmaxfactor": ode_rmaxfactor, "pde_nmesh": pde_nmesh, "ode_rstep": ode_rstep, "ode_nsteps": ode_nsteps, "ode_thresh": ode_thresh, "pde_thresh": pde_thresh, "pde_maxiter": pde_maxiter}

        # get the ODE
        if not oper:
            if self.onlyonedifferential:
                # solve the PDE
                metric = ExtendedSelfDualityEquation(K=K, coefs=self.Pcoefsadjusted[K-2], method = method, pde_nmesh = pde_nmesh, pde_thresh = pde_thresh, pde_maxiter = pde_maxiter, rmax = rmax)
                # build the ODE from the solution of the PDE
                self.params["iter"] = metric.iter
                self.params["delta"] = metric.delta
                self.params["error"] = metric.error
                ODE = planarODE.buildhiggsODE(metric)
                if storemetric:
                    self.metric = metric
            else:
                raise NotImplementedError("Case of multiple differentials for Higgs bundle not implemented yet")
        if oper:
            ODE = planarODE.buildoperODEfrompolys(Pcoefs = self.Pcoefsadjusted, K = self.K, r = rmax)

        # solve the ODE to get the vertices
        self.computeverticesfromODE(ODE)

        # find total time elapsed
        totaltime = time.time() - start_time
        logger.info("Finished in %0.1f s" % totaltime)
        self.params["totaltime"] = totaltime

    def updateToCurrent(self):
        if not hasattr(self, "dataversion") or self.dataversion < 13:
            raise NotImplementedError("dataversion = %d < 13 not supported" % self.dataversion)

        old_dataversion = self.dataversion

        if old_dataversion != self.dataversion:
            logger.info("Updated framedata version from %d to %d" % (old_dataversion, self.dataversion))

    def computeverticesfromODE(self, ODE):
        rmax = self.params["rmax"]
        ode_rmaxfactor = self.params["ode_rmaxfactor"]
        ode_rstep = self.params["ode_rstep"]
        ode_nsteps = self.params["ode_nsteps"]
        ode_thresh = self.params["ode_thresh"]

        # TODO: extend this all to work with general differentials
        nvert = self.degree+self.K
        # Initial angle for integration rays
        star_angle = -(np.angle(self.coefs0[-1])-(self.K*self.theta+pi))/nvert
        self.phaselist = [star_angle + n*(2*pi/nvert) for n in range(nvert)]

        # ODE solved here along *nvert* evenly spaced rays in C
        logger.info("Solving ODE with initial angle = %0.4f*pi, r = %0.2f" % (star_angle/pi,rmax))
        self.framelist = ODE.integrate_rays(n=nvert,theta0=star_angle,r=rmax*ode_rmaxfactor,step=ode_rstep,tol=ode_thresh,nsteps=ode_nsteps,return_type='frame')

        self.buildvertexlist()

    @staticmethod
    def sortedeigensystem(m):
        # based on https://stackoverflow.com/questions/10083772/python-numpy-sort-eigenvalues
        evals, evecs = np.linalg.eig(m)
        evals_sorted = evals[abs(evals).argsort()]
        evecs_sorted = evecs[:, abs(evals).argsort()]
        return evals_sorted, evecs_sorted

    def buildvertexlist(self):
        self.vertexlist = []
        for F in self.framelist:
            w,v = self.sortedeigensystem(F)
            if abs(w[0]) < 10 and abs(w[0])>0.10:
                logger.warning("Eigenvalue %s too close to 1: full eigenvalue list %s" % (w[0],w))
            self.vertexlist.append(v[:,0])
        self.verticesnormalized = False

    def eigenvalues(self):
        evals = []
        for F in self.framelist:
            w,v = np.linalg.eig(F)
            evals.append(w)
        return evals

    def getCluster(self):
        """Return the cluster variables appropriate to this theory"""
        return computefdcluster(self.theory, self.vertexlist, realitycheck = self.real)

    def estimateODEvertexerror(self):
        '''Estimate the absolute error in each component of vertexlist based on ODE relative error goal
        Uses the linearized change in eigenvector under a matrix perturbation.
        Return array of size (nvert,K) representing component error estimates'''
        nvert = len(self.framelist)
        est_vtx_move = np.zeros(shape=(nvert,self.K),dtype='float64')
        for i,F in enumerate(self.framelist):
            w,v = self.sortedeigensystem(F)
            vdual = np.linalg.inv(v).transpose()
            dvec = np.zeros(shape=(self.K,),dtype='float64')
            for k in range(1,self.K):
                lambda_diff=abs(w[k]-w[0])
                delta_dvec = np.dot(np.abs(vdual[:,k]),np.dot(np.abs(self.params['ode_thresh']*F),np.abs(v[:,0]))) * np.abs(v[:,k]) / lambda_diff
                dvec += delta_dvec
            est_vtx_move[i] = dvec
        return est_vtx_move

    def numderivclusterA(self):
        '''Compute numerical partial derivative of the A-cluster vector as a function of each component of each vertex
        Returns array of size (nvert,K,Arank) where entry (i,j,k) is d(Acluster[k])/(dvertexlist[i][j])'''
        nvert = len(self.framelist)
        deriv_eps = 1e-6
        verts0 = np.array(self.vertexlist)
        if self.real:
            dt = 'float64'
            F = lambda L:np.real(np.array(self.theory['fdAclusterfunc'](L)))
        else:
            dt = 'complex128'
            F = lambda L:np.array(self.theory['fdAclusterfunc'](L))
        derivAcluster = np.zeros(shape=(nvert,self.K,self.theory['Arank']),dtype=dt)
        for i in range(nvert):
            for j in range(self.K):
                verts1 = verts0.copy()
                verts1[i,j] += deriv_eps 
                derivAcluster[i,j] = (F(verts1) - F(verts0)) / deriv_eps
        return derivAcluster

    def estimateODEclusterAerror(self):
        '''Estimate the absolute error in each cluster A-variable due to ode_thresh'''
        dA = self.numderivclusterA()  # complex-valued
        dvert = self.estimateODEvertexerror()  # positive real valued (already abs of estimated difference)
        est_A_move = np.tensordot(np.abs(dA),dvert,axes=((0,1),(0,1)))
        return est_A_move

    def estimateODEclusterXerror(self):
        '''Estimate the absolute error in each cluster X-variable due to ode_thresh'''
        B = np.array(self.theory["B"])
        Xbasis = np.array(self.theory["Xbasis"])
        XBT = np.dot(Xbasis,np.transpose(B))
        est_X_move = np.zeros(shape=(self.theory['Xrank'],),dtype='float64')
        est_A_move = self.estimateODEclusterAerror()
        Acluster = self.theory['fdAclusterfunc'](self.vertexlist)
        if self.real:
            Acluster = np.real(Acluster)
        Xcluster = self.getCluster()
        est_X_move = []
        for X,Xcoefv in zip(Xcluster,XBT):
            dlogasum = sum( est_A_move[j]/abs(Acluster[j]) for j in range(self.theory['Arank']) if Xcoefv[j]!=0 )
            dX = abs(X)*dlogasum
            est_X_move.append(dX)
        return np.array(est_X_move)

    def estimateODEclustererror(self):
        '''Get the absolute error estimate (due to ode_thresh) for either cluster'''
        return self.estimateODEclusterXerror()

    # save object
    def save(self, filename):
        """Save framedata to a file (with jsonpickle)"""
        logger.info("Saving data to %s" % filename)
        frozen = jsonpickle.encode(self, warn=True, make_refs=False)
        with open(filename, "w") as target:
            target.write(frozen)

    @staticmethod
    def load(filename):
        """Load an instance from a jsonpickle file"""
        logger.info("Loading data from %s" % filename)
        with open(filename, "r") as target:
            frozen = target.read()
        frozen.replace("crossratios", "framedata")
        framedata = jsonpickle.decode(frozen)
        framedata.updateToCurrent()
        return framedata

def clusterAtoX(theory, Acluster, realitycheck = True):
    Xrank = theory["Xrank"]
    if Xrank == 0:
        return []
    B = np.array(theory["B"])
    Xbasis = np.array(theory["Xbasis"])
    XBT = np.dot(Xbasis,np.transpose(B))
    if realitycheck:
        # just force the cluster variable to be negative, i.e. control sign "by hand"
        return list(-np.exp(np.dot(XBT, np.log(np.abs(Acluster)))))
    else:
        # sign here also controlled "by hand", should be done in a better way
        Xsignlist = theory["Xsigns"]
        return list(-np.exp(np.dot(XBT, np.log(Acluster)))*Xsignlist)

def computefdcluster(theory, vertexlist, theoryname = None, realitycheck = True):
    if theory is None:
        theory = theorydata.getdata(theoryname)

    Acluster = theory["fdAclusterfunc"](vertexlist)   

    if realitycheck:
        eps = 1e-4
        for A in Acluster:
            if abs(A.imag) > eps:
                logger.warning("Large imaginary part for cluster A-coordinate: %s" % A)
        Acluster = [A.real for A in Acluster]

    return clusterAtoX(theory,Acluster,realitycheck = realitycheck)

# TODO: It is not ideal that R is a required positional arg, but then ignored if oper==True
def loadFrameData(theoryname, R, theta = 0, tag = "", oper=False, absh = 1, pde_nmesh = PDE_NMESH):
    """load frame data for a specific theory and parameters from FRAMEPATH"""
    fn = fd_filename(theoryname=theoryname,R=R,theta=theta,tag=tag,oper=oper,absh=absh,fullpath=True,pde_nmesh=pde_nmesh)
    fd = framedata.load(fn)
    return fd

def loadFrameDataCluster(theoryname, R, theta = 0, tag = "", oper=False, absh = 1, pde_nmesh = PDE_NMESH):
    """load frame data for a specific theory and parameters from FRAMEPATH and return a cluster"""
    fn = fd_filename(theoryname=theoryname,R=R,theta=theta,tag=tag,oper=oper,absh=absh,fullpath=True,pde_nmesh=pde_nmesh)
    fd = framedata.load(fn)
    return fd.getCluster()


def computeFrames(theoryname, absh = 1.0, R = 1, theta = 0, storemetric = False, rmax = None, pde_nmesh = PDE_NMESH, pde_thresh = PDE_THRESH, ode_rmaxfactor = ODE_RMAXFACTOR, ode_thresh = ODE_THRESH, pde_maxiter = PDE_MAXITER, ode_nsteps = ODE_NSTEPS, ode_rstep = ODE_RSTEP, make_monic = True, method = "fourier", oper = False):
    """Compute frame/vertex data for a specific theory and parameters"""
    th = theorydata.getdata(theoryname)
    if not oper:
        logger.info("Starting to compute frames for theory %s, R = %0.3E, theta = %0.3E, pde_nmesh = %d" % (theoryname,R,theta,pde_nmesh))
    else:
        logger.info("Starting to compute frames for theory %s, oper, abs(h) = %0.8E, theta = %0.3E" % (theoryname,absh,theta))
    fd = framedata(theory = th, absh = absh, R = R, theta = theta, storemetric = storemetric, pde_nmesh = pde_nmesh, pde_thresh = pde_thresh, pde_maxiter = pde_maxiter, make_monic = make_monic, rmax = rmax, ode_rmaxfactor = ode_rmaxfactor, ode_nsteps = ode_nsteps, ode_rstep = ode_rstep, ode_thresh = ode_thresh, theoryname = theoryname, method = method, oper = oper)
    return fd

def computeAndSaveFrames(theoryname, R = 1, theta = 0, storemetric = False, pde_nmesh = PDE_NMESH, pde_maxiter = PDE_MAXITER, pde_thresh = PDE_THRESH, ode_thresh = ODE_THRESH, method = "fourier", oper = False, make_monic = True, absh = 1.0, **kwargs):
    """Compute frame/vertex data for a specific theory and parameters and save to FRAMEPATH
    Return: (framedata computed, filename written)"""
    fd = computeFrames(theoryname, R = R, theta = theta, storemetric = storemetric, pde_nmesh = pde_nmesh, pde_maxiter = pde_maxiter, pde_thresh = pde_thresh, ode_thresh = ode_thresh, method = method, oper = oper, make_monic = make_monic, absh = absh, **kwargs)
    fn = fd_filename(theoryname=theoryname, R=R, theta=theta, oper=oper, absh=absh, fullpath=True, pde_nmesh=pde_nmesh)
    fd.save(fn)
    return fd, fn
