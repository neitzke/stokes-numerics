from __future__ import absolute_import
from __future__ import print_function
import math
import cmath
import scipy
import scipy.integrate
import scipy.interpolate
import scipy.special
import scipy.fftpack
import pprint
import itertools
import numpy as np
import random
import operator
import time
import datetime
import os
import re

from codeversion import codeversion

import jsonpickle
import jsonpickle.ext.numpy as jsonpickle_numpy
jsonpickle_numpy.register_handlers()
from cmath import pi,exp,log
import namegen

import fileinfo
import theorydata
import theory

import logging
logger = logging.getLogger(__name__)


#----------------------------------------------------------------------
# Default values for key parameters in the IEQ calculations
IEQ_L = 200
IEQ_TOLERANCE = 2e-15
IEQ_DAMPING = 0.3
IEQ_STEPS = 2**17
#----------------------------------------------------------------------

class activeray:
    """
    Class for storing data about an active ray.

    members:
    id (int) = a unique id for this ray
    phase (complex) = unit-norm complex number for the position of the ray
    function (function) = ray-function of three arguments (gamma, charges, xlists), for use in the integral equations
    activedetector (function) = function of one argument (gamma), just detects whether this ray gives a contribution for X_gamma or not
    mutuallylocal (bool) = keep track of whether this ray is attached to a single BPS charge
    charge = list of charges for which we'll need to compute x on this ray
    boundaryfin (function) = function of one argument (gamma), returns a fixed quantity to be included in this ray's integrand
    """

    def __init__(self, id, phase, function, activedetector, mutuallylocal, charges, boundaryfin):
        self.id = id
        self.phase = phase
        self.function = function
        self.activedetector = activedetector
        self.mutuallylocal = mutuallylocal
        self.charges = charges
        self.boundaryfin = boundaryfin

class raydatum:
    """Class for storing new-style information about the x-functions on
    a single BPS ray."""

    def __init__(self, rayid, L, steps, tlist, xinstlists, xsflists, mutuallylocal, charges):
        self.rayid = rayid
        self.L = L
        self.steps = steps
        self.tlist = tlist
        self.xinstlists = xinstlists
        self.xsflists = xsflists
        self.mutuallylocal = mutuallylocal
        self.charges = charges


class xarray:
    """Main class for storing information about the x-functions attached to
       BPS rays.
     """

    def __init__(self, theory, raydata = None, parent = None, method = "fourier", steps = IEQ_STEPS, L = IEQ_L):
        """Creates a new xarray object from the given data. 

           If "raydata" is specified, then the xarray is populated with that data.
           (This generally happens when __init__ is called after computing an iteration.)
           Otherwise the xarray is populated with initial data corresponding to the
           0-th iteration.

           If a "parent" is specified, then the delta from the 
           parent is computed, and stored as "self.delta".
           Otherwise, no computations are done here.

           members: 
             theory (theory object)
             codeversion
             dataversion
             timestamp
             contactpotential
             contactpotentialquantum
             contactpotentialclassical
             method
             raydata (list of raydatum object)
             parentHash
             delta
             iter
        """
        self.theory = theory
        self.codeversion = codeversion
        self.dataversion = 15
        self.timestamp = datetime.datetime.utcnow()
        self.method = method

        if raydata is None:
            self.raydata = []
            self.parentHash = None
            self.delta = None
            self.iter = 0

            for ray in self.theory.activerays:
                charges = []
                xinstlists = []
                xsflists = []
                for charge in ray.charges:
                    tlist, xinstlist, xsflist = self.theory.initialfunctions(charge, L, steps = steps, phase = ray.phase)

                    charges.append(charge)
                    xinstlists.append(xinstlist)
                    xsflists.append(xsflist)

                newrayitem = raydatum(rayid = ray.id, L = L, steps = steps, tlist = tlist, xinstlists = xinstlists, xsflists = xsflists, mutuallylocal = True, charges = charges)
                self.raydata.append(newrayitem)

        else:
            self.raydata = raydata
            self.parentHash = hash(parent)
            self.delta = self.compareToXar(parent)
            self.iter = parent.iter + 1
        
    def updateToCurrent(self):
        """Update xarray to latest version of data structure, or complain if
        it's too old."""
        if not hasattr(self, "dataversion") or self.dataversion < 15:
            raise NotImplementedError("dataversion = %d < 15 not supported" % self.dataversion)

        old_dataversion = self.dataversion

        if old_dataversion != self.dataversion:
            logger.info("Updated xar version from %d to %d" % (old_dataversion, self.dataversion))

    def getRaydatum(self, rayid):
        for raydatum in self.raydata:
            if raydatum.rayid == rayid:
                return raydatum
        raise KeyError("Desired ray datum not found")

    def getrayfunction(self, rayid, gamma, part = "real"):
        raydatum = self.getRaydatum(rayid)
        xlists = [xinst+xsf for xinst, xsf in zip(raydatum.xinstlists,raydatum.xsflists)]
        return self.theory.activerays[rayid].function(gamma, raydatum.charges, xlists)

    def getboundaryfinlist(self, rayid, gamma):
        raydatum = self.getRaydatum(rayid)
        return self.theory.activerays[rayid].boundaryfin(gamma, raydatum.tlist)

    def tlist(self, charge = None, rayid = None):
        """Get the tlist for a specific charge.
           arguments: charge (list of ints)
           returns: list of floats"""
        return self.getRaydatum(rayid).tlist

    def computeNextRaydatum(self, ray, damping, method):
        """Compute next ray datum for a specific ray.
        """
        if method == "fourier":
            return self.computeNextRaydatumFourier(ray, damping = damping)
        if method == "simps":
            return self.computeNextRaydatumSimps(ray, damping = damping)
        else:
            raise NotImplementedError("Asked for non-implemented method '%s'" % method)

    def computeNextRaydatumSimps(self, ray, damping):
        currentraydatum = self.getRaydatum(ray.id)
        L = currentraydatum.L
        steps = currentraydatum.steps
        tlist = currentraydatum.tlist
        charges = currentraydatum.charges

        xinstlistsold = currentraydatum.xinstlists
        xsflists = currentraydatum.xsflists

        zetaphase = ray.phase
        zetalist = zetaphase*np.exp(tlist)

        xinstlists = []
        for charge,xinstlistold,xsflist in zip(charges,xinstlistsold,xsflists):
            xinstlist = (1-damping)*self.computexinstlist(charge, zetalist, method = "simps", onray = True) + (damping)*xinstlistold
            xinstlists.append(xinstlist)

        return raydatum(rayid = ray.id, L = L, steps = steps, tlist = tlist, xsflists = xsflists, xinstlists = xinstlists, mutuallylocal = ray.mutuallylocal, charges = charges)

    def computeNextRaydatumFourier(self, ray, damping = IEQ_DAMPING):
        currentraydatum = self.getRaydatum(ray.id)
        L = currentraydatum.L
        steps = currentraydatum.steps
        tlist = currentraydatum.tlist
        charges = currentraydatum.charges

        xinstlistsold = currentraydatum.xinstlists
        xsflists = currentraydatum.xsflists

        zetaphase = ray.phase

        xinstlists = []
        # run over all charges whose X-functions are relevant on this ray
        # for each one, need to compute new X-function
        for charge,xinstlistold in zip(charges,xinstlistsold):

            # we want to build new xinstlist; start from old one
            xinstlist = damping*xinstlistold

            # contributions to new xinst will be sum over all active rays
            # for which "activedetector(charge)" is not zero
            for equivclass in self.theory.rayequivclasses:
                rayid = equivclass[0]
                activeray = self.theory.getraybyid(rayid)
                if not activeray.activedetector(charge):
                    continue
                numinclass = len(equivclass)

                # compute Fourier transform of ray function
                rayflist = self.getrayfunction(rayid, charge)
                ftrayflist = scipy.fftpack.fft(rayflist)

                # compute Fourier transform of boundary-condition function
                boundaryfinlist = self.getboundaryfinlist(rayid, charge)
                ftboundaryfinlist = scipy.fftpack.fft(boundaryfinlist)

                # next, compute convolution kernel
                lzb = log(zetaphase/activeray.phase)
                # conventions: log(0) = -0, log(-1) = +i pi
                if abs(lzb.imag + pi) < 0.000001:
                    lzb = lzb.real + pi*1j

                # list of Fourier momenta along the line
                # TODO: check for off-by-1 error -- note endpoints in position space are L and -L
                p = -(pi*(steps-1)/L)*(scipy.fftpack.fftfreq(steps, 1).astype(np.complex_))

                # Fourier transformed kernel, computed as a function of the Fourier momentum p
                # this kernel depends on how zeta sits relative to the active ray
                # if zeta is on the ray, we take the limit as zeta approaches the ray counterclockwise
                if lzb.imag > 0:
                    # piecewise function to avoid overflow problems
                    kernel = np.piecewise(p, [p.real > 150, p.real <= 150], [lambda pp: 0.0, lambda pp: -2*pi*1j*np.exp((-pp*1.0j)*lzb)/(1 + np.exp(pi*pp))])
                else:
                    # piecewise function to avoid overflow problems
                    kernel = np.piecewise(p, [p.real < -150, p.real >= -150], [lambda pp: 0.0, lambda pp: 2*pi*1j*np.exp((-pp*1.0j)*lzb)/(1 + np.exp(-pi*pp))])

                ftconvolved = kernel * (ftrayflist + ftboundaryfinlist)

                # inverse Fourier transform to get the convolution
                convolved = scipy.fftpack.ifft(ftconvolved)

                contribution = 1.0/(4*pi*1j) * convolved

                xinstlist += (1-damping)*numinclass*contribution

            xinstlists.append(xinstlist)

        return raydatum(rayid = ray.id, L = L, steps = steps, tlist = tlist, xsflists = xsflists, xinstlists = xinstlists, mutuallylocal = ray.mutuallylocal, charges = charges)


    def computeNextIteration(self, damping, method):
        """Compute one iteration, by iteration over equivalence classes
           and calling computeNextRaydata for one charge in each class.
           arguments: damping (float), method (str) = "simps" or "fourier"
           returns a new xarray"""
        newraydata = []

        for equivclass in self.theory.rayequivclasses:
            logger.debug("Computing ray data for ray class %s" % equivclass)

            # choose the first ray in this class, use it for computing
            rayid = equivclass[0]
            ray = self.theory.getraybyid(rayid)
            # no optimizations yet in non-symmetric case
            newraydatum = self.computeNextRaydatum(ray, damping = damping, method = method)
            newraydata.append(newraydatum)

        newxarray = xarray(theory = self.theory, raydata = newraydata, parent = self, method = method)

        return newxarray

    def getCluster(self, theta = 0, absh = 1.0, method = "simps"):
        """Get cluster (list of floats or complex numbers).
           arguments: theta (float), method (str) = "simps" """
        logger.debug("Calling getCluster with method = %s" % method)
        def co(charge):
            if self.theory.oper:
                zeta = absh*exp(theta*1.0j)
            else:
                zeta = exp(theta*1.0j)
            return self.computeX(charge=charge, zeta=zeta, method = method, realpartonly = None)
        return computexarcluster(self.theory, co, theta, realitycheck = not self.theory.oper)

    def computexinstlist(self, charge, zetalist, method = "simps", realpartonly = None, onray = False):

        eps = 1e-9

        if method is None:  # use the same method we used for computing the current xarray
            method = self.method
        if method == "fourier":
            # logger.warning("fourier method not implemented for computing a single value of xinst: switched to simps instead")
            method = "simps"
        if method != "simps":
            raise NotImplementedError("Only simps method supported")
        if realpartonly is None:
            if len(zetalist) == 1:
                zeta = zetalist[0]
                realpartonly = (abs(zeta) - 1 < 1e-9) and (self.theory.symmetric) and (not self.theory.oper)
            else:
                realpartonly = False

        # TODO: remember why realpartonly is set identically to False
        realpartonly = False

        ntocompute = len(zetalist)
        xinstlist = np.zeros(ntocompute,dtype=np.complex_)

        for equivclass in self.theory.rayequivclasses:
            # pick a ray in this class
            rayid = equivclass[0]
            ray = self.theory.getraybyid(rayid)
            numinclass = len(equivclass)

            # skip this ray if its contribution would be zero as measured
            # by "rayactivedetector" function
            if ray.activedetector(charge):
                # now fetch data for this ray from the current xar (self)
                raydatum = self.getRaydatum(ray.id)
                raytlist = raydatum.tlist
                raycharges = raydatum.charges
                steps = raydatum.steps
                rayxlists = [ rayxsflist + rayxinstlist for (rayxsflist,rayxinstlist) in zip(raydatum.xsflists,raydatum.xinstlists) ]
                L = raydatum.L
                boundaryfinlist = self.getboundaryfinlist(ray.id, charge)

                rayzetalist = ray.phase*np.exp(raytlist)

                if steps % 2 == 1:  # because we didn't implement principal value here, to be safe, want steps even
                    raise ValueError("simps requires even number of steps")
                stepsize = 2.0*L/(steps-1)

                for n in range(ntocompute):
                    zeta = zetalist[n]

                    def symmetricintegrandfunction(zeta, rayzetalist, charge, raycharges, rayxlists):
                        return (2*zeta*rayzetalist / (zeta**2 - rayzetalist**2)) * (ray.function(charge, raycharges, rayxlists) + boundaryfinlist)

                    def unsymmetricintegrandfunction(zeta, rayzetalist, charge, raycharges, rayxlists):
                        return (rayzetalist - zeta)/(rayzetalist + zeta) * ray.function(charge, raycharges, rayxlists)

                    if ((abs(zeta / (abs(zeta)*ray.phase) - 1) < eps and self.theory.symmetric) or (abs(zeta / (abs(zeta)*ray.phase) + 1) < eps)) and (abs(zeta) - 1 > eps):

                        # need to use principal-value integration
                        if onray:
                            # we're deliberately trying to recompute xinst at all sample points  along the ray, so zetalist came from tlist and 
                            # the singularity occurs exactly at the n-th sampling point
                            # thus we can approximate principal value by deleting some range symmetrically around this point
                            # experiments suggest using gapwidth = 0 works best
                            # TODO: make a better approximation of the principal value

                            gapwidth = 0
                            range1end = max(n-gapwidth,0)
                            range2begin = min(n+1+gapwidth,ntocompute)
                        else:
                            # we're "accidentally" computing along the ray; this isn't well implemented yet
                            raise NotImplementedError("principal value integration not implemented except for ray iteration or |zeta| = 1")

                        rayzetalist1 = rayzetalist[:range1end]
                        rayzetalist2 = rayzetalist[range2begin:]
                        rayxlists1 = [l[:range1end]   for l in rayxlists]
                        rayxlists2 = [l[range2begin:] for l in rayxlists]

                        if self.theory.symmetric:
                            integrand1 = symmetricintegrandfunction(zeta, rayzetalist1, charge, raycharges, rayxlists1)
                            integrand2 = symmetricintegrandfunction(zeta, rayzetalist2, charge, raycharges, rayxlists2)

                            if (abs(zeta / (abs(zeta)*ray.phase) - 1) < eps):
                                residuesign = -1
                            if (abs(zeta / (abs(zeta)*ray.phase) + 1) < eps):
                                residuesign = 1

                            residuepart = (-1*pi*1j) * residuesign * ray.function(charge, raycharges, [l[n] for l in rayxlists])
                        else:
                            integrand1 = unsymmetricintegrandfunction(zeta, rayzetalist1, charge, raycharges, rayxlists1)
                            integrand2 = unsymmetricintegrandfunction(zeta, rayzetalist2, charge, raycharges, rayxlists2)

                            residuepart = (-1*pi*1j) * 2 * ray.function(charge, raycharges, [l[n] for l in rayxlists])

                        if len(integrand1)==0:
                            integrand1 = np.zeros(1)
                        if len(integrand2)==0:
                            integrand2 = np.zeros(1)


                        # do the integrals
                        imintegral = scipy.integrate.simps(integrand1.imag, x = None, dx = stepsize) + scipy.integrate.simps(integrand2.imag, x = None, dx = stepsize)
                        if realpartonly:
                            reintegral = 0.0
                        else:
                            reintegral = scipy.integrate.simps(integrand1.real, x = None, dx = stepsize) + scipy.integrate.simps(integrand2.real, x = None, dx = stepsize)

                    else:  # don't need principal-value integration
                        if self.theory.symmetric:
                            integrand = symmetricintegrandfunction(zeta, rayzetalist, charge, raycharges, rayxlists)
                        else:
                            integrand = unsymmetricintegrandfunction(zeta, rayzetalist, charge, raycharges, rayxlists)
                        residuepart = 0

                        # do the integrals
                        imintegral = scipy.integrate.simps(integrand.imag, x = None, dx = stepsize)
                        if realpartonly:
                            reintegral = 0.0
                        else:
                            reintegral = scipy.integrate.simps(integrand.real, x = None, dx = stepsize)

                    # finally evaluate the full contribution from this ray
                    integral = reintegral + imintegral*1j + residuepart

                    integralcontribution = 1/(4*pi*1j) * integral
                    xinstlist[n] += numinclass*integralcontribution

        return xinstlist

    def computexinst(self, charge, zeta, method = None, realpartonly = None):
        # just call computexinstlist with a list of a single element
        xinstlist = self.computexinstlist(charge, zetalist = [zeta], method = method, realpartonly = realpartonly)
        return xinstlist[0]

    def computex(self, charge, zeta, method = "simps", realpartonly = False):
        """Compute the x function. Just passes through to computexinst
           and then adds the semiflat part. See computexinst for
           details about arguments, methods."""
        sfvalue = self.theory.xsf(charge, zeta)
        return self.computexinst(charge = charge, zeta = zeta, method = method, realpartonly = realpartonly)+sfvalue

    def computeX(self, charge, h = None, zeta = None, method = None, realpartonly = False, twisted = True):
        """Compute the X function. Just passes through to computex
           and exponentiates. See computexinst for details about arguments,
           methods."""

        # the code deeper in the stack generally uses the name "zeta" for the parameter, whether or not we are dealing with opers;
        # but the user will want to name the parameter "h" in the oper case; we allow this just by setting zeta = h here
        if self.theory.oper and zeta is None:
            zeta = h

        if zeta is None:
            raise ArithmeticError("No argument provided for X function in computeX()")

        if twisted:
            sign = self.theory.sigma(charge)
        else:
            sign = 1
        return sign*exp(self.computex(charge = charge, zeta = zeta, method = method, realpartonly = realpartonly))

    def computexinstonray(self, charge, t, method = None):
        """Compute xinst on the BPS ray. Just passes through to computexinst
           with zeta = -Z/|Z| * exp(t).
           See computexinst for details about arguments, methods."""
        zetaphase = self.theory.Z(charge) / abs(self.theory.Z(charge))
        zeta = -zetaphase*exp(t)
        return self.computexinst(charge = charge, zeta = zeta, method = method)

    def computeXonray(self, charge, t, method = None):
        """Compute X on the BPS ray. Just passes through to computexinst
           with zeta = -Z/|Z| * exp(t).
           See computexinst for details about arguments, methods."""
        zetaphase = self.theory.Z(charge) / abs(self.theory.Z(charge))
        zeta = -zetaphase*exp(t)
        return self.computeX(charge = charge, zeta = zeta, method = method)

    def compareToXar(self, other):
        """Compute L^infinity norm of difference between this xar and
           another. For simplicity, this is only implemented when the
           two xars have same L and steps on each BPS ray. This will
           be true in the usual case of comparing iteration to its parent.

           arguments: other (xarray)

           Returns list of list of floats."""
        norms = []

        for raydatum in self.raydata:
            raydatumother = other.getRaydatum(raydatum.rayid)

            xinstlists = raydatum.xinstlists
            xinstlistsother = raydatumother.xinstlists
            L = raydatum.L
            Lother = raydatumother.L
            steps = raydatum.steps
            stepsother = raydatumother.steps

            if steps != stepsother:
                raise ArithmeticError("Trying to compare xars with different steps")
            if abs(L - Lother) > 1e-12:
                raise ArithmeticError("Trying to compare xars with different L")

            curnorms = [np.abs(xinstlist-xinstlistother).max() for xinstlist,xinstlistother in zip(xinstlists,xinstlistsother)]
            norms.append(curnorms)

        return norms

    # get maximum delta
    def maxdelta(self):
        """Get the maximum value of "delta"; this is a rough estimate of the difference
           between this xar and its parent.
           returns: float"""
        if len(self.delta) == 0:
            return 0
        else:
            return max(max(self.delta))

    def minsteps(self):
        """Get the minimum value of "steps" over all elements of the raydata list.
           returns: int"""
        stepslist = [item.steps for item in self.raydata]
        if len(stepslist) == 0:
            return 0
        else:
            return min(stepslist)

    def minL(self):
        """Get the minimum value of "L" over all elements of the raydata list.
           returns: float"""
        Llist = [item.L for item in self.raydata]
        if len(Llist) == 0:
            return 0
        else:
            return min(Llist)

    # save to file
    def save(self, filename):
        """Save an xarray to a file, using jsonpickle.
           See computeAndSaveXar() for a higher-level interface that will run the calculation
           and save to a file, choosing the filename automatically.
           arguments: filename (str)"""
        logger.info("Saving xar data to %s" % filename)
        frozen = jsonpickle.encode(self, warn=True)

        with open(filename, "w") as target:
            target.write(frozen)

    def saveclusters(self, filename, theta = 0.0, absh = 1.0):
        logger.info("Saving cluster data to %s" % filename)
        metadata = { 
            "params": self.params, 
            "codeversion": self.codeversion, 
            "theta": theta, 
            "oper": self.theory.oper, 
            "theoryname": self.theory.theoryname
        }
        if self.theory.oper:
            metadata["absh"] = absh
        else:
            metadata["R"] = self.theory.R
        cluster = {"X": self.getCluster(theta = theta, absh = absh), "A": None, "metadata": metadata}

        frozen = jsonpickle.encode(cluster, warn=True)
        with open(filename, "w") as target:
            target.write(frozen)

    @staticmethod
    def load(filename, theoryname, R):
        """Load an xarray from a file, using jsonpickle. Generally we don't call this directly.

           arguments: filename (str)
           returns: xarray"""

        logger.info("Loading xar data from %s" % filename)
        with open(filename, "r") as target:
            frozen = target.read()
        xar = jsonpickle.decode(frozen)
        xar.updateToCurrent()

        # reconstruct the theory object, since it can't be saved/loaded
        xar.theory = theory.theory(theoryname = theoryname, R = R)

        return xar

    @staticmethod
    def loadclusters(filename):
        logger.info("Loading cluster data from %s" % filename)
        with open(filename, "r") as target:
            frozen = target.read()
        clusters = jsonpickle.decode(frozen)
        return clusters

    def testSample(self, charge, x, method = None):
        """Test the fixed-point equation for a given charge, at a fixed sample point.
        The sample point is specified by x in [-1,1]. This is converted
        into a value t = -L + n*stepsize (rounding to get n), then xinst 
        is computed at t along the BPS
        ray, using computexinstonray, with a given integration method specified by "method".
        If the convergence is good this should be very close to xinstlist[n]. Both
        values are returned, along with their difference.

        arguments: charge (list of ints), n (int), method (str) 
        returns: tuple of complex numbers (x1,x2,x1-x2)
          x1 and x2 are the computed value and the current value of xinst at t""" 

        L = self.L(charge)
        xinstlist = self.xinstlist(charge)
        steps = self.steps(charge)
        stepsize = 2*L/(steps-1)
        n = int( ((x + 1.0)/(2.0))*steps )
        if not (0 <= n < steps):
            raise ValueError("testSample called with n=%d outside range [0,%d)" % (n,steps))
        t = -L + n*stepsize
        x1 = self.computexinstonray(charge, t, method = method)
        x2 = xinstlist[n]
        return (x1,x2,x1-x2)


def computexarcluster(theory, co, theta, realitycheck = True):
    """Compute a cluster (returned as list of float or complex numbers) from values of X functions
       arguments: theoryname (str), theta (float) determine the combinatorics of the cluster
                  co (function) is a function returning X for a given charge, which we use as input"""
    cluster = theory.data["xarXclusterfunc"](co, theta)

    eps = 1e-4
    if realitycheck:
        for X in cluster:
            if abs(X.imag) > eps:
                logger.warning("Large imaginary part for cluster X-coordinate: %s" % X)
        return [X.real for X in cluster]
    else:
        return cluster


def xdiff(xar1, xar2, charge):
    """Get difference between xinst for two different xars and a given charge.
       These xars need not have the same L's or steps: the difference
       is computed using linear interpolation. It is returned as a pair of lists:
       tlist (values of t where the difference was computed),
       xlist (xalues of the difference, xinst1 - xinst2)
       The number of steps used is the maximum of the numbers of steps in the input xars.
       The value of L used is the minimum of the L for the two input xars.

       arguments: xar1 (xarray), xar2 (xarray), charge (list of ints)
       returns: tuple (tlist,xlist)
         tlist: list of floats
         xlist: list of complex"""
    xinstlist1 = xar1.xinstlist(charge)
    xinstlist2 = xar2.xinstlist(charge)
    tlist1 = xar1.tlist(charge)
    tlist2 = xar2.tlist(charge)
    xinst1 = scipy.interpolate.interp1d(tlist1, xinstlist1, fill_value = 0, bounds_error = False)
    xinst2 = scipy.interpolate.interp1d(tlist2, xinstlist2, fill_value = 0, bounds_error = False)
    L1 = xar1.L(charge)
    L2 = xar2.L(charge)
    samples = max(len(tlist1),len(tlist2))
    L = min(L1,L2)
    tlist = np.linspace(-L,L,samples)
    xlist = [xinst1(t) - xinst2(t) for t in tlist]
    return (tlist,xlist)

def getApproxCluster(theoryname, R = 1, theta = 0, nterms = 0, th = None, oper = False, absh = 1.0):
    """Compute an approximate cluster (returned as list of float or complex numbers)
                  nterms (int) controls how many terms to take in the approximation;
                  nterms = 0 means just the leading exponential"""
    if th is None:
        th = theory.theory(theoryname = theoryname, R = R, oper = oper)
    if oper:
        zeta = absh*exp(theta*1j)
    else:
        zeta = exp(theta*1j)
    def cocorr(charge):
        return th.Xcorr(charge, zeta = zeta, nterms = nterms, oper = oper)
    return computexarcluster(th, cocorr, theta, realitycheck = not oper)

def semiflatXar(theory):
    xar = xarray(theory, raydata = None, method = "fourier", steps = IEQ_STEPS, L = IEQ_L)
    xars = [xar]

    return (xar,xars)

def computeXar(theoryname, R = 1, L = IEQ_L, tolerance = IEQ_TOLERANCE, steps = IEQ_STEPS, method = "fourier", damping = IEQ_DAMPING, oper = False, nosymmetric = False, maxiter = 1000, splitrays = True, iterattoleranceneeded = 5, usebc = False, failonmaxiter = True, zerorayfunctions = False):
    """Compute a single xar.
       returns: xarray"""

    if steps not in [2**n for n in range(30)]:
        raise ValueError("steps must be of the form 2**n")

    if method == "simps" and steps > 4000:
        logger.warning("Using simps method and steps = %d: this will take a long time" % steps)

    curtheory = theory.theory(theoryname = theoryname, R = R, oper = oper, nosymmetric = nosymmetric, splitrays = splitrays, usebc = usebc, zerorayfunctions = zerorayfunctions)
    if not oper:
        logger.info("Starting to compute xars for theory %s, R = %0.3E, tolerance = %0.3E, steps = %d, maxiter = %s, L = %s, damping = %s" % (curtheory.theoryname,curtheory.R,tolerance,steps,maxiter,L,damping))
    else:
        logger.info("Starting to compute xars for theory %s, oper, tolerance = %0.3E, steps = %d, maxiter = %s, L = %s, damping = %s" % (curtheory.theoryname,tolerance,steps,maxiter,L,damping))

    # initialize the xar
    xar = xarray(curtheory, raydata = None, method = method, steps = steps, L = L)

    maxdelta = 100

    totaltime = 0.0
    iterattolerance = 0
    niter = 0
    while iterattolerance < iterattoleranceneeded:

        if niter == maxiter:
            if failonmaxiter:
                raise ArithmeticError("hit maxiter = %d without reaching tolerance" % maxiter)
            else:
                logger.warning("Hit maxiter = %d without reaching tolerance: last maxdelta = %0.3E" % (maxiter,maxdelta))
                break

        if maxdelta <= tolerance:
            iterattolerance += 1
            # logger.info("%d / %d iterations at tolerance so far" % (iterattolerance,iterattoleranceneeded))
        else:
            iterattolerance = 0

        # start the clock
        start_time = time.time()

        # now run the next iteration
        nextxar = xar.computeNextIteration(damping = damping, method = method)
        xar = nextxar
        niter += 1

        # update maxdelta for next round
        maxdelta = nextxar.maxdelta()

        # display timing and other info
        elapsed = time.time() - start_time
        totaltime += elapsed
        logger.info("IEQ: iteration %d, maxdelta = %0.3E (%0.1f s)" % (niter,maxdelta,elapsed))

        # abort if maxdelta too big
        if math.isnan(maxdelta):
            raise ArithmeticError("maxdelta grew out of control")

    logger.info("Finished in %0.1f s" % totaltime)
    if not xar.theory.oper and False:  # disabled because of switch to new raydata; TODO: re-implement
        xar.computeContactPotential()

    xar.params = {
                    "tolerance": tolerance,
                    "iterattoleranceneeded": iterattoleranceneeded,
                    "damping": damping,
                    "L": L,
                    "splitrays": splitrays,
                    "method": method,
                    "maxiter": maxiter,
                    "totaltime": totaltime
    }

    return xar

def computeAndSaveXar(theoryname, R = 1.0, L = IEQ_L, tolerance = IEQ_TOLERANCE, steps = IEQ_STEPS, method = "fourier", damping = IEQ_DAMPING, oper = False, saveclusters = True, savexar = False, thetalist = [0.0], abshlist = [1.0], splitrays = True, iterattoleranceneeded = 5, maxiter = 1000):
    """Compute and save to disk a single xar and/or its clusters.
    Return: (xar computed, filename written)
    Note: If both xar and cluster files are requested, returns only the xar file."""
    xar = computeXar(theoryname = theoryname, R = R, L = L, tolerance = tolerance, steps = steps, method = method, damping = damping, oper = oper, maxiter = maxiter, iterattoleranceneeded = iterattoleranceneeded, splitrays = splitrays)
    if saveclusters:
        for theta in thetalist:
            if oper:
                for absh in abshlist:
                    fn = namegen.ie_filename(theoryname, R, oper=oper, theta=theta, clusteronly=True, fullpath=True, absh=absh)
                    xar.saveclusters(fn, theta=theta, absh=absh)
            else:
                fn = namegen.ie_filename(theoryname, R, oper=oper, theta=theta, clusteronly=True, fullpath=True)
                xar.saveclusters(fn, theta=theta)
    if savexar:
        fn = namegen.ie_filename(theoryname, R, oper=oper, clusteronly=False, fullpath=True)
        xar.save(fn)
    if not savexar and not saveclusters:
        logger.warning("computeAndSaveXar() was called with savexar=False and saveclusters=False.  NOTHING WILL BE SAVED TO DISK.")
    return xar, fn
    
def loadXar(theoryname, R = 1.0, tag = "", oper = False):
    """Load xar from file.

       NB: R gets truncated to fixed number of decimal digits (8)
       
       returns: xarray"""
    fn = namegen.ie_filename(theoryname=theoryname,R=R,tag=tag,oper=oper,fullpath=True)
    xar = xarray.load(fn, theoryname = theoryname, R = R)
    return xar

def loadXarCluster(theoryname, R = 1.0, absh = 1.0, tag = "", oper = False, theta = 0.0):
    """Load clusters from file.  Or if no cluster file exists but the xar file does, load the
    xar and compute the requested clusters."""

    # if we have a .xarcluster file then load from there
    fn = namegen.ie_filename(theoryname=theoryname,R=R,absh=absh,tag=tag,oper=oper,clusteronly=True,theta=theta,fullpath=True)
    if namegen.ie_file_exists(theoryname=theoryname,R=R,absh=absh,tag=tag, clusteronly=True, oper=oper, theta=theta):
        clusters = xarray.loadclusters(fn)
        if clusters["X"] != None:
            return clusters["X"]

    # otherwise load from .xar file and compute
    logger.info("loadXarCluster did not find \"%s\"; loading full xar and computing clusters instead." % fn)
    xar = loadXar(theoryname = theoryname, R = R, oper = oper)
    return xar.getCluster(theta = theta)
