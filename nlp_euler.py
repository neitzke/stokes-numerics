"""Forward Euler non-linear Poisson solver"""

# This module is concerned with solving the "non-linear Poisson"
# equation
#     Delta(u) = f(u,z)
# on a uniform rectangular mesh, with u = u0 on the boundary.
#
# We solve the equation by an iterative method, solving the linearized
# equation at u_i to get u_{i+1} and terminating when u_{i+1} - u_i is
# small enough.
#
# Using a five-point stencil for the discrete Laplacian (i.e. centered
# finite difference), the linearization becomes a sparse linear
# system.  We solve this linear system using scipy.linalg.spsolve.
#
# This solver gives high accuracy in relatively few iterations, but
# suffers from high setup time and high memory consumption.  In
# practice it is often memory bound, e.g. requiring more than 3.5 GiB
# of RAM for a 1000x1000 mesh.
#
# KNOWN ISSUES:
#
# * The boundary conditions are not handled properly in full
#   generality; the current implementation only works when the initial
#   function u_0 satisfies the equation in a neighborhood of the
#   boundary.  That is, the difference Delta(u_0) - f(u_0,z) should
#   vanish identically near the boundary.  (This assumption holds for
#   the applications to the self-duality equations for which the
#   solver was developed.)

from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import scipy.sparse.linalg
from solverexception import SolverException
import discdiff
import time

import logging
logger = logging.getLogger(__name__)

class NLPEuler(object):
    """Solve the system Delta(u) = f(u,z) on a uniform mesh, with u = u0
    on boundary.  Use forward Euler iteration with u0 as the initial
    guess.
    """
    
    def __init__(self,f,d1f,u0,grid,thresh=0.0000001,maxiter=15,relax=0.9,linear=False):
        """Initialize and run the solver.

        Parameters:
          f    : f(u,z)
          d1f  : [df/du](u,z)
          u0   : initial guess and boundary conditon
          grid : SquareGrid or similar object representing a rectangular
                 mesh (use zm, nx, ny, dx, dy attributes) 
          thresh : Absolute error goal for the solver
          maxiter : raise exception if threshold not met after this many iterations
          relax : Step by (relax)*(linearized solution) at each
                  iteration; setting to less than 1.0 may enlarge
                  domain of convergence at the cost of convergence
                  speed.
          linear : is the equation to be solved actually linear?
                   if so, just need one iteration

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
        self.linear = linear
        if linear:
            self.relax = 1.0
        
        self.lapmat = discdiff.discdiff(self.grid, kind = "laplacian")

        self._t1 = time.time()
        self._t = self._t1

        logger.info("Solving PDE")
        self.u = self._iterate()

    def _iterate(self):
        """Forward Euler iteration loop"""

        # CONVENTIONS: 1) Whenever a grid of values is converted to a
        # flat vector, we append "vec" to the name.  For quantities
        # that are only ever manipulated as flat vectors, this does
        # not apply.
        # 2) We don't return flat vectors from public methods.
        
        # Warning: v.ravel() returns a VIEW of the array v, so making
        # a copy is necessary if we don't want to modify v when
        # accessing this flattened vector.
        uvec = np.copy(self.u0.ravel())
        
        n = 0
        logger.info('PDE: GOAL=%g' % self.thresh)
        while n < self.maxiter:
            n = n + 1
            udotvec, b = self._step(uvec)
            bnorm = max(abs(b))
            udotnorm = max(abs(udotvec))
            uvec = uvec + self.relax * udotvec

            now = time.time()
            logger.info('PDE: iter=%d error=%g delta=%g\t(%.2fs)' % (n,bnorm,udotnorm,now-self._t))
            self._t = now

            if bnorm < self.thresh:
                logger.info('PDE: success\t(%.2fs total; %.2fs in main loop)',now - self._t0,now - self._t1)
                break

            if self.linear:
                logger.info('PDE: finished solving linear equation \t(%.2fs total; %.2fs in main loop)',now - self._t0,now - self._t1)
                break
        else:
            # raise SolverException('Max iterations (%d) reached without meeting error threshold (%g)' % (self.maxiter,self.thresh))
            logger.info("PDE: reached max iterations\t(%.2fs total; %.2fs in main loop)",now - self._t0,now - self._t1)

        u = uvec.reshape( (self.grid.nx, self.grid.ny) )

        self.iter = n
        self.delta = udotnorm
        self.error = bnorm

        return u

    def _step(self,uvec):
        """Solve the linearization of Delta(u) = f(u,z) at u.  Return both the
        solution (udot, as a flat vector) and the right hand side (b)
        of the linearization equation.

        """
        # Note that the linearization is:
        #     Delta(udot) - d1f(u,z) udot = f(u,z) - Delta(u)
        # which we turn into the standard linear system
        #     A x = b
        # where
        #     A = (Delta - d1f(u,z))
        #     x = udot
        #     b = f(u,z) - Delta(u)

        A = self.lapmat - self._diag(self.d1f(uvec,self.grid.zv))
        b = self.f(uvec,self.grid.zv) - self.lapmat * uvec

        # TODO: This linear system is only nonsingular because we've
        # handled the boundary conditions poorly!  In practice it
        # seems to work well if f(0,x) = 0.  
        udotvec = scipy.sparse.linalg.spsolve(A,b) # Solves sparse linear system A x = b; returns x

        return udotvec, b
    
    def _diag(self,v):
        return scipy.sparse.spdiags(v, [0], self.grid.N, self.grid.N )
