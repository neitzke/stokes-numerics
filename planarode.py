from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import scipy.integrate
import contextlib
import warnings

from polynomial_utils import poly

import logging
logger = logging.getLogger(__name__)

#----------------------------------------------------------------------
#  BEGIN pre-computed flat connection matrices for Higgs bundles
#----------------------------------------------------------------------
#
# CONVENTIONS
#
# Infinitesimal invariants are passed around as "connection data"
# vectors, which are 5-tuples consisting of
#
# (u(z), du/dx(z), du/dy(z), Re(P(z)), Im(P(z)))
#
# For K=3
# The affine frame field is considered as a 3x3 matrix with ROWS
#   f, f_x, f_y

def coefmat_higgs_dx_K2(a):
    """dx component of connection in rank 2 as a function of the
    connection data vector 'a'

    """
    u = np.ndarray.item(a[0])
    ux = np.ndarray.item(a[1])
    uy = np.ndarray.item(a[2])
    Pre = a[3]
    Pim = a[4]

    euhalf = np.exp(0.5*u)
    eu = np.exp(u)
    emu = np.exp(-u)

    Ax = np.array([ [eu - emu*Pre, -emu*Pim - 0.5*uy],
                    [-emu*Pim + 0.5*uy, -eu + emu*Pre] ])
    return -Ax

def coefmat_higgs_dy_K2(a):
    """dy component of connection in rank 2 as a function of the
    connection data vector 'a'

    """
    u = np.ndarray.item(a[0])
    ux = np.ndarray.item(a[1])
    uy = np.ndarray.item(a[2])
    Pre = a[3]
    Pim = a[4]

    euhalf = np.exp(0.5*u)
    eu = np.exp(u)
    emu = np.exp(-u)

    Ay = np.array([ [emu*Pim, -eu - emu*Pre + 0.5*ux],
                    [-eu - emu*Pre - 0.5*ux, -emu*Pim] ])
    return -Ay

def coefmat_higgs_dx_K3(a):
    """dx component of connection in rank 3 as a function of the
    connection data vector 'a'

    """
    u = np.ndarray.item(a[0])
    ux = np.ndarray.item(a[1])
    uy = np.ndarray.item(a[2])
    Pre = a[3]
    Pim = a[4]

    euhalf = np.exp(0.5*u)
    eu = np.exp(u)
    emu = np.exp(-u)

    sqrt2 = 2**0.5
    Ax = np.array([ [0.0, sqrt2*euhalf, 0.0],
                    [sqrt2*euhalf, -Pre*emu, Pim*emu + 0.5*uy],
                    [0.0, Pim*emu - 0.5*uy, Pre*emu] ])
    return -Ax

def coefmat_higgs_dy_K3(a):
    """dy component of connection in rank 3 as a function of the
    connection data vector 'a'

    """
    u = np.ndarray.item(a[0])
    ux = np.ndarray.item(a[1])
    uy = np.ndarray.item(a[2])
    Pre = a[3]
    Pim = a[4]

    euhalf = np.exp(0.5*u)
    eu = np.exp(u)
    emu = np.exp(-u)

    sqrt2 = 2**0.5
    Ay = np.array([ [0.0, 0.0, sqrt2*euhalf],
                    [0.0, Pim*emu, Pre*emu - 0.5*ux],
                    [sqrt2*euhalf, Pre*emu + 0.5*ux, -Pim*emu] ])
    return -Ay

#----------------------------------------------------------------------
#  END pre-computed flat connection matrices for Higgs bundles
#----------------------------------------------------------------------


#----------------------------------------------------------------------
#  BEGIN pre-computed flat connection matrices for opers
#----------------------------------------------------------------------
#
# CONVENTIONS
#
# "connection data" a is tuple of complex numbers [P2(z), P3(z), ..., PK(z)]

def coefmat_oper_dx_K2(a):
    """dx component of connection in rank 2 as a function of the
    connection data vector 'a'

    """
    P2 = a[0]

    Ax = np.array([ [0, -P2],
                    [1, 0] ])
    return -Ax

def coefmat_oper_dy_K2(a):
    """dy component of connection in rank 2 as a function of the
    connection data vector 'a'

    """
    P2 = a[0]

    Ay = np.array([ [0, -P2*1j],
                    [1j, 0] ])
    return -Ay

def coefmat_oper_dx_K3(a):
    """dx component of connection in rank 3 as a function of the
    connection data vector 'a'

    """
    P2 = a[0]
    P3 = a[1]

    Ax = np.array([ [0, -0.5*P2, -P3],
                    [1, 0, -0.5*P2],
                    [0, 1, 0] ])

    return -Ax

def coefmat_oper_dy_K3(a):
    """dy component of connection in rank 3 as a function of the
    connection data vector 'a'

    """
    P2 = a[0]
    P3 = a[1]

    Ay = np.array([ [0, -0.5j*P2, -1j*P3],
                    [1j, 0, -0.5j*P2],
                    [0, 1j, 0] ])
    return -Ay


#----------------------------------------------------------------------
#  END pre-computed flat connection matrices for opers
#----------------------------------------------------------------------


def coefmat(xpart, ypart, theta, a):
    """Convert dx and dy components of a connection to the dt component
    for a ray of angle theta in the plane"""
    return np.cos(theta)*xpart(a) + np.sin(theta)*ypart(a)

def odeA(K, xpart, ypart, z, zprime, conn_data, y):
    """Convert a connection integration problem to the
    callback style used by scipy ODE solvers.
    Parameters:
        K : rank
        xpart : Coefficient of dx in connection (function conn. data -> matrix)
        ypart : Coefficient of dy in connection (function conn. data -> matrix)
        z : current value of z
        zprime : current value of zprime
        conn_data : Function z -> connection data
        y : Flat vector from current solution matrix"""
    m = abs(zprime)*coefmat(xpart, ypart, np.angle(zprime), conn_data(z))
    ya = np.array(y)
    ya.shape = (K, K)
    return m.dot(ya).ravel() # Compute Y' as A*Y, then flatten KxK matrix to K^2-vector

@contextlib.contextmanager
def _make_fatal(category):
    with warnings.catch_warnings():
        warnings.filterwarnings("error",category=category)
        yield None        

class planarODE():
    """Class representing a planar ODE."""

    # TODO: if we want to be more general about which kind of ODEs we can 
    # take, then r shouldn't be a member, it should be passed to functions
    # that want to integrate along rays

    DEFAULT_RET_TYPE='affine'

    def __init__(self,K,xpart,ypart,conn_data,r,real=False,**kwargs):
        self.K = K
        self.xpart = xpart
        self.ypart = ypart
        self.conn_data = conn_data
        self.r = r
        self.yinit = np.eye(self.K).ravel()
        self.real = real

    def _prepare_output(self,F,return_type='affine'):
        if return_type == 'affine':
            return F[0,1:] / F[0,0]
        elif return_type == 'homog':
            return F[0,:]
        elif return_type == 'frame':
            return F

    def integrate_ray(self,r=None,theta=0.0,step=0.01,fixstep=True,tol=0.00001,nsteps=50000,method="dopri5",return_type='affine'):
        if r is None:
            r = self.r

        alpha = np.exp(1j*theta)
        odef = lambda t,y: odeA(self.K, self.xpart, self.ypart, z = t*alpha, zprime = alpha, conn_data = self.conn_data, y = y)

        if self.real:
            solver = scipy.integrate.ode(odef)
        else:
            solver = scipy.integrate.complex_ode(odef)

        si_kwargs = dict()
        if fixstep:
            si_kwargs['max_step'] = 2.0*step*r

        if method=="dopri5":
            solver.set_integrator("dopri5",first_step=(step*r),rtol=tol,nsteps=nsteps,**si_kwargs)
        elif method=="adams":
            solver.set_integrator("vode",method="adams",with_jacobian=False,first_step=(step*r),rtol=tol,nsteps=nsteps,**si_kwargs)
        else:
            raise ValueError("Unknown ODE method: {}".format(method))

        solver.set_initial_value(self.yinit,0.0)
        logger.info("Calling integrator at angle %0.4f*pi, r = %0.2f" % (theta/np.pi,r))
        # ODE Solver raises UserWarning when things go wrong, and
        # typically returns bad results in that case.  We upgrade
        # this warning to an exception to prevent using such data.
        with _make_fatal(UserWarning):
            solver.integrate(r)
        return self._prepare_output(solver.y.reshape((self.K,self.K)), return_type)

    def integrate_point(self,z,**kwargs):
        return self.integrate_ray(r=abs(z),theta=np.angle(z),**kwargs)

    def integrate_rays(self,n,r=None,theta0=0.0,**kwargs):
        thetalist = np.linspace(theta0, theta0 + 2*np.pi, num=n, endpoint=False)
        return [ self.integrate_ray(r=r,theta=theta,**kwargs) for theta in thetalist ]

    def integrate_points(self,zlist,**kwargs):
        return [ self.integrate_point(z,**kwargs) for z in zlist ]

    def integrate_path(self,zfunc,zprimefunc,step=0.001,fixstep=True,tol=0.00001,method="dopri5",nsteps=50000):
        odef = lambda t,y: odeA(self.K, self.xpart, self.ypart, zfunc(t), zprimefunc(t), self.conn_data, y)
        solver = scipy.integrate.complex_ode(odef)
        si_kwargs = dict()
        if fixstep:
            si_kwargs['max_step'] = 2.0*step
        if method=="dopri5":
            solver.set_integrator("dopri5",first_step=step,rtol=tol,nsteps=nsteps,**si_kwargs)
        elif method=="adams":
            solver.set_integrator("vode",method="adams",with_jacobian=False,first_step=step,rtol=tol,nsteps=nsteps,**si_kwargs)
        else:
            raise ValueError("Unknown ODE method: {}".format(method))
        solver.set_initial_value(self.yinit,0.0)
        with _make_fatal(UserWarning):
            solver.integrate(1.0)
        return self._prepare_output(solver.y.reshape((self.K,self.K)), "frame")

    @staticmethod
    def buildoperODEfromfunctions(Pfuncs, K, r):
        def _conn_data(z):
            return [Pfunc(z) for Pfunc in Pfuncs]
        if K == 2:
            xpart, ypart = coefmat_oper_dx_K2, coefmat_oper_dy_K2
        elif K == 3:
            xpart, ypart = coefmat_oper_dx_K3, coefmat_oper_dy_K3
        else:
            raise NotImplementedError("Oper ODE only implemented for K=2 or K=3")
        return planarODE(K, xpart, ypart, _conn_data, r)

    @staticmethod
    def buildoperODEfrompolys(Pcoefs, K, r):
        Pfuncs = [(lambda z,coefs=coefs: poly(coefs,z)) for coefs in Pcoefs]
        return planarODE.buildoperODEfromfunctions(Pfuncs, K, r)

    @staticmethod
    def buildhiggsODE(sdesol, r = None):
        if r is None:
            r = 0.95 * sdesol.grid.r

        if sdesol.K == 2:
            xpart, ypart = coefmat_higgs_dx_K2, coefmat_higgs_dy_K2
        elif sdesol.K == 3:
            xpart, ypart = coefmat_higgs_dx_K3, coefmat_higgs_dy_K3
        else:
            raise NotImplementedError("Higgs bundle ODE only implemented for K=2 or K=3")
        return planarODE(K = sdesol.K, xpart = xpart, ypart = ypart, conn_data = sdesol._conn_data, r = r)
