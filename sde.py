"""Hitchin's self duality equations for the real cyclic case, rank 2
or 3, on complex plane"""

# The self-duality equation (SDE) for a cyclic higgs bundle with Higgs
# field determined by the holomorphic K-differential
#   p(z) dz^K
# can be written in the form
#   Laplacian(u) = f(u,z)
# where u is the log-density of the metric solving the SDE and where
# f(u,z) depends on K and p.  This module implements those functions
# and builds an initial guess for the solution, ultimately solving the
# equation using one of the non-linear Poisson solvers.

from __future__ import absolute_import
from __future__ import print_function
import numpy as np
from nlp_euler import NLPEuler
from nlp_fourier import NLPFourier

import logging
logger = logging.getLogger(__name__)

# Right hand sides f(u,z) and derivatives (df/du)(u,z) of SDE written
# as a nonlinear Poisson equation for rank K=2

def f_K2(p,u,z):
    ps = np.abs(p(z))**2
    return (4.0*np.exp(2.0*u) - 4.0*np.exp(-2.0*u)*ps)

def d1f_K2(p,u,z):
    ps = np.abs(p(z))**2
    return (8.0*np.exp(2.0*u) + 8.0*np.exp(-2.0*u)*ps)

# Right hand sides f(u,z) and derivatives (df/du)(u,z) of SDE written
# as a nonlinear Poisson equation for rank K=3

def f_K3(p,u,z):
    ps = np.abs(p(z))**2
    return (4.0*np.exp(u) - 4.0*np.exp(-2.0*u)*ps)

def d1f_K3(p,u,z):
    ps = np.abs(p(z))**2
    return (4.0*np.exp(u) + 8.0*np.exp(-2.0*u)*ps)

def smoothed_step(t):
    return np.piecewise(t,[t<0.0,t>1.0],[0.0,1.0,lambda x:0.5*(1.0 - np.cos(np.pi*x))])

def uinitfunc(K, req, p, z):
    '''Smoothed semiflat metric for initial guess; is exactly the semiflat metric when |z| > req'''
    ps = np.abs(p(z))**2
    rho = smoothed_step(1.0 - (np.abs(z) / req))
    if K == 3:
        return (1.0/3.0)*np.log(rho*np.exp(-ps*ps)+ps)
    if K == 2:
        return (1.0/4.0)*np.log(rho*np.exp(-ps*ps)+ps)

class SelfDualityEquation(object):
    """Solve Hitchin's self-duality equation for real cyclic case, rank 2 or 3, on complex plane"""
    def __init__(self,K,p,grid,thresh=1e-7,maxiter=5000,req=None,method='fourier'):
        """Initialize and run the solver.

        Parameters:
          K : rank of the Higgs bundle (2 or 3)
          p : Callable function, p(z) dz^K is the holomorphic
              differential appearing in the Higgs field.
          grid : SquareGrid or similar object representing a rectangular mesh
          thresh : Error goal for the PDE solver (L^2 or L^inf depending on solver)
          maxiter : raise exception if threshold not met after this many iterations
          req : Controls a radius used in setting the initial guess
                for the solution of the self-duality equations.  For
                |z| > req we set the initial guess to be exactly equal
                to |p(z)dz^K|^(2/K).  For |z| <= req we use a smoothed
                version of this.  One should have |z_i| < req < grid.r
                for all zeros z_i of p(z).  If None, used 0.95*grid.r
          method : PDE solver to use; see SUPPORTED METHODS below.

        Return:
          None

        Output class attributes:
          u : Solution to self-duality equations
          u0 : Initial guess used in solving self-duality equations
        
        Supported methods:
          'euler' : Forward Euler iteration
          'fourier' : Fourier transform quasi-linearization iteration
        """

        self.K = K
        self.p = p
        self.grid = grid
        self.thresh = thresh
        self.maxiter = maxiter
        self.method = method

        if req == None:
            self.req = 0.9*grid.r
        else:
            self.req = req
        
        if self.K==2:
            f = lambda u,z: f_K2(p,u,z)
            d1f = lambda u,z: d1f_K2(p,u,z)
        elif self.K==3:
            f = lambda u,z: f_K3(p,u,z)
            d1f = lambda u,z: d1f_K3(p,u,z)
        else:
            raise NotImplementedError('rank K must be 2 or 3')

        u0 = lambda z: uinitfunc(K,self.req,p,z)

        if method == 'euler':
            NLPSolver = NLPEuler
        elif method == 'fourier':
            NLPSolver = NLPFourier
        else:
            raise ValueError('Unknown solver method "%s"' % method)
            
        S = NLPSolver(f,d1f,u0,self.grid,thresh=self.thresh,maxiter=self.maxiter,relax=1.0,linear=False)

        self.iter = S.iter
        self.delta = S.delta
        self.error = S.error

        self.u0 = S.u0 # evaluated initial guess (not the function)
        self.u = S.u

