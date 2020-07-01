"""Solve Hitchin's self-duality equation for real cyclic case, rank 2
or 3, on the complex plane and set up for the corresponding flat connection.

"""
from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import scipy.interpolate
import scipy.ndimage

from squaregrid import SquareGrid
from sde import SelfDualityEquation
from planarode import planarODE
from polynomial_utils import poly, derivcoefs
from functools import partial
from math import pi

import logging
logger = logging.getLogger(__name__)

import fileinfo
import os
import jsonpickle
import jsonpickle.ext.numpy as jsonpickle_numpy
jsonpickle_numpy.register_handlers()


class ExtendedSelfDualityEquation(SelfDualityEquation):
    """Solve Hitchin's self-duality equation for real cyclic case, rank 2
       or 3, on the complex plane.
    """
    def __init__(self,K,coefs,rmax,pde_nmesh,pde_thresh,pde_maxiter,method=None):
        """Initialize and solve the self-duality equations.

        Keyword parameters:
          See SelfDualityEquations.__init__

        Return:
          None

        Output class attributes:
          Grids of values:
            u : Solution to self-duality equations (SDE)
            u0 : Initial guess used in solving self-duality equations
            ux : du/dx
            uy : du/dy
          Bivariate interpolating functions:
            uinterp : Solution to SDE
            uxinterp : du/dx
            uyinterp : du/dy
        """

        self.K = K
        self.coefs = coefs
        self.rmax = rmax
        self.method = method
        self.degree = len(coefs)-1
        
        nvert = self.degree+self.K

        # create the grid
        self.grid = SquareGrid(rmax,pde_nmesh)

        # set p to a function which evaluates the polynomial with coefficients self.coefs
        p = partial(poly,self.coefs)

        harm_kwargs = {'thresh': pde_thresh, 'maxiter': pde_maxiter}
        # Pass method parameter only if set
        if method:
            harm_kwargs['method'] = method

        SelfDualityEquation.__init__(self,K,p,self.grid,**harm_kwargs) # call the PDE solver

        # TODO/note: We ignore the discrete laplacian module's finite difference operators
        # and use scipy.ndimage.convolve1d instead.
        self.ux = scipy.ndimage.convolve1d(self.u,[0.5,0,-0.5],axis=1) / self.grid.dx
        self.uy = scipy.ndimage.convolve1d(self.u,[0.5,0,-0.5],axis=0) / self.grid.dy

        self.uinterp = scipy.interpolate.RectBivariateSpline(self.grid.x,self.grid.y,self.u.transpose()) #,kx=2,ky=2)
        self.uxinterp = scipy.interpolate.RectBivariateSpline(self.grid.x,self.grid.y,self.ux.transpose()) #,kx=1,ky=1)
        self.uyinterp = scipy.interpolate.RectBivariateSpline(self.grid.x,self.grid.y,self.uy.transpose()) #,kx=1,ky=1)

    def _conn_data(self,z):
        pval = self.p(z)
        return np.array( [self.uinterp(z.real,z.imag),
                          self.uxinterp(z.real,z.imag),
                          self.uyinterp(z.real,z.imag),
                          pval.real,
                          pval.imag] )


    # save object
    def save(self, filename):
        """Save to a file (with jsonpickle)"""
        logger.info("Saving data to %s" % filename)
        frozen = jsonpickle.encode(self, warn=True, make_refs=False)
        with open(os.path.join(fileinfo.FRAMEPATH,filename), "w") as target:
            target.write(frozen)

    @staticmethod
    def load(filename):
        """Load an instance from a jsonpickle file"""
        logger.info("Loading data from %s" % filename)
        with open(fileinfo.FRAMEPATH+"/"+filename, "r") as target:
            frozen = target.read()
        harmonicmap = jsonpickle.decode(frozen)
        return harmonicmap
