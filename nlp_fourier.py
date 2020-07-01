"""Fourier transform non-linear Poisson solver"""

# This module is concerned with solving the "non-linear Poisson"
# equation
#     Delta(u) = f(u,z)
# on a uniform rectangular mesh, with u = u0 on the boundary.
#
# We solve the equation by an iterative method, solving an
# approximation to the linearized equation at u_i to get u_{i+1} and
# terminating when u_{i+1} - u_i is small enough.
#
# The key feature of this solve is that we use a very coarse
# approximation of the linearization---chosen specifically so that it
# can be solved by Fourier transform methods.  The coarse
# approxmination means that each iteration makes little progress
# toward the final solution, and many iterations are necessary.
# However, the availability of efficient FFT routines means that each
# iteration is very fast, and so in many cases there is a net gain
# compared to a direct method.
#
# The exact linearized equation for v = u-u0 is
#     Delta(vdot) - d1F(v,z) vdot = F(v,z) - Delta(vdot)  (*)
# where
#     F(v,z) = f(u0+v,z) - Delta(u0)
# We rewrite (*) as
#     (Delta - A)vdot = RHS
# This is exactly solvable by Fourier methods if A is a constant
# function.
#
# To approximate a solution, we replace A = d1F(v,z) by a constant
# that is in some way representative of its values on he grid points.
# We follow the suggestion of [1] to use the "minimax" value
#
#     A = (max(d1F) + min(d1F)) / 2
#
# where max and min are taken over the grid.
#
# References
#
# [1] Concus, P. and Golub, G. H. 1973. Use of fast direct methods for
# the efficient numerical solution of nonseparable elliptic
# equations. SIAM J. Numer. Anal., 10: 1103-1103.
#
# KNOWN ISSUES:
#
# * The initialization code assumes that u_0 is harmonic in a
#   neighborhood of the boundary of the mesh.  This is not a
#   fundamental requirement of the method, but because u_0 cannot be
#   easily extended to a doubly-periodic function its Laplacian is
#   computed by a finite difference scheme rather than by FFT methods.
#   Being harmonic at the boundary allows us to simply zero out the
#   Laplacian at the edges and ignore this issue.
#
#   (Note that this assumption is satisfied for the applications to
#   the self-duality equations for which this solver was developed0).

from __future__ import absolute_import
import numpy as np
import scipy.signal
from dst2 import dst2, idst2, dst2freq
from solverexception import SolverException
import time

import logging
logger = logging.getLogger(__name__)

def _max_power_2_dividing(n):
    n = int(n)
    return n & (~(n-1))

def _suggest_sizes(n):
    if n == _max_power_2_dividing(n):
        return [n-1]
    a = np.log(n+1)/np.log(2.0)
    return [2**int(np.floor(a))-1, 2**int(np.ceil(a))-1]

def _is_bad_size_for_dst(n):
    return float(n+1) / _max_power_2_dividing(n+1) > 5.0

class NLPFourier(object):
    """Solve the system Delta(u) = f(u,z) on a uniform mesh, with u = u0
    on boundary, using Fourier transform methods.

    """
    def __init__(self,f,d1f,u0,grid,thresh=0.0000001,maxiter=5000,relax=1.0,linear=False):
        """Initialize and run the solver.

        Parameters:
          f    : f(u,z)
          d1f  : [df/du](u,z)
          u0   : initial guess and boundary conditon
          grid : SquareGrid or similar object representing a rectangular
                 mesh (use zm, nx, ny, dx, dy attributes) 
          thresh : L^2 error goal for the solver
          maxiter : raise exception if threshold not met after this many iterations
          relax : Step by (relax)*(linearized solution) at each
                  iteration; setting to less than 1.0 may enlarge
                  domain of convergence at the cost of convergence
                  speed.
          linear : is the equation to be solved actually linear?
                   (not used)

        Return:
          None

        Output class attributes:
          u  : Solution u
          u0 : Initial guess
        """

        self._t0 = time.time()

        self.f = f
        self.d1f = d1f
        self.grid = grid
        self.u0func = u0
        self.u0 = u0(self.grid.zm)
        self.thresh = thresh
        self.maxiter = maxiter
        self.relax = relax

        self.normcoef = self.grid.dx * self.grid.dy / (2.0*np.sqrt((self.grid.nx + 1)*(self.grid.ny+1)))

        self.warn_if_bad_sizes()
        
        # Capital K for unnormalized frequencies ("per index")
        # Lower k for real frequencies ("per unit x or y")
        KX, KY = dst2freq(self.u0)
        self.kx = 2*np.pi*KX/self.grid.dx  # a row
        self.ky = 2*np.pi*KY/self.grid.dy  # a column
        self.ksq = self.kx**2 + self.ky**2 # broadcasted to a 2D array now
        
        # Laplacian of initial guess and its transform
        # (These are relatively long-running computations.)
        # TODO 1: Remove implicit assumption that Delta(u0) vanishes on boundary. 
        # TODO 2: Move this finite difference stuff into its own module.
        logger.info("Computing Laplacian of initial guess and its transform")
        idx,idy = 1.0/self.grid.dx, 1.0/self.grid.dy
        lap_stencil = np.array([[0,idy**2,0],[idx**2,-2*(idx**2 + idy**2),idx**2],[0,idy**2,0]],dtype=np.float64)
        Lap_u0_raw = scipy.signal.convolve2d(self.u0,lap_stencil,mode='same')
        self.Lap_u0 = np.zeros_like(Lap_u0_raw)
        self.Lap_u0[2:-2,2:-2] = Lap_u0_raw[2:-2,2:-2]

        self._t1 = time.time()
        self._t = self._t1

        logger.info("Solving PDE: %dx%d grid, thresh=%g" % (self.grid.nx,self.grid.ny,self.thresh))
        self.u = self._iterate()
        
    def _iterate(self):
        """Fourier solver main loop"""
        vhat = np.zeros_like(self.u0)

        n = 0
        last_delta_norm = 0.0
        
        while True:
            n = n + 1

            # Compute the DST of the RHS of the inhomogeneous linearized equation
            v = idst2(vhat)
            vvec = v.reshape((self.grid.nx * self.grid.ny, ))
            Fvvec = self.F(vvec,self.grid.zv)
            Fv = Fvvec.reshape((self.grid.nx, self.grid.ny))
            Fv_hat = dst2(Fv)
            Lapv_hat = -self.ksq * vhat
            RHS_hat = Fv_hat - Lapv_hat

            # Compute the L^2 norm of the inhomogeneous term
            # ( = 0 iff we have a solution )
            residual = self.L2norm_hat(RHS_hat)

            now = time.time()
            logger.info("PDE: iter=%d L2error=%g L2delta=%g\t(%.2fs)" % (n,residual,last_delta_norm,now-self._t))
            self._t = now

            if residual < self.thresh:
                logger.info('PDE: success\t(%.2fs total; %.2fs in main loop)',now - self._t0,now - self._t1)
                break

            if np.any(np.isnan(RHS_hat)):
                # Computing RHS revealed some failure in the computation
                # (usually means the linearized solution at the previous step was bad)
                raise SolverException("NAN encountered in RHS computation (overflow or underflow?)")

            if n >= self.maxiter:
                raise SolverException("Max iterations (%d) reached without meeting error threshold %g" % (self.maxiter,self.thresh))

            # Solve a constant-coefficient approximation of the linear
            # equation in frequency space.

            # First compute the constant that approximates d1F.
            a = self.minimax(self.d1F(vvec,self.grid.zv))

            # Now compute the transform of the exact solution to this
            # constant coef problem.
            delta_vhat =  RHS_hat / (-self.ksq - a)
            last_delta_norm = self.L2norm_hat(delta_vhat)

            # Update vhat by adding this approx solution of the
            # linearization
            vhat = vhat + self.relax * delta_vhat

        self.iter = n
        self.delta = last_delta_norm
        self.error = residual
        
        return self.u0 + v

    def warn_if_bad_sizes(self):
        msgstr = '%s-size %d is a bad choice for the fourier solver; this computation will be inefficient.  Good sizes have the form (2**n)-1.  Consider using size %s instead.'
        if _is_bad_size_for_dst(self.grid.nx):
            logger.warning(msgstr % ('x',self.grid.nx, ' or '.join(str(n) for n in _suggest_sizes(self.grid.nx))))
        if _is_bad_size_for_dst(self.grid.ny):
            logger.warning(msgstr % ('y',self.grid.ny, ' or '.join(str(n) for n in _suggest_sizes(self.grid.nx))))
    
    def F(self,v,z):
        u0vec = self.u0.reshape((self.grid.nx * self.grid.ny, ))
        Lapu0vec = self.Lap_u0.reshape((self.grid.nx * self.grid.ny, ))
        return self.f(u0vec + v,z) - Lapu0vec

    def d1F(self,v,z):
        u0vec = self.u0.reshape((self.grid.nx * self.grid.ny, ))
        return self.d1f(u0vec + v,z)

    def minimax(self,m):
        """Return the average of max(m) and min(m)"""
        return 0.5*(np.max(m) + np.min(m))

    def L2norm_hat(self,m):
        """L^2 norm of a function computed from its Fourier transform coefficients"""
        return np.linalg.norm(m) * self.normcoef
