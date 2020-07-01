'''Compute vertices of polygon (hyperbolic ideal or real projective convex) corresponding to a meromorphic cyclic Higgs budnle'''

# Imports needed to process arguments
from __future__ import absolute_import
from __future__ import print_function
import sys
import argparse
import numpy as np

# This was adapted from ngon.py from Blaschke

parser = argparse.ArgumentParser(description=__doc__)

parser.add_argument('-N',type=int,default=2,help='Rank of Higgs bundle')
parser.add_argument('coefs',type=complex,nargs='+',help='If using --coefs, list of polynomial coefficients, constant term first, complex() format, defining the holomorphic N-differential')
parser.add_argument('-r','--roots',action='store_true',help='instead of coefficients, roots are given')
parser.add_argument('-t','--theta',type=float,default=0.0,help='multply by exp(-i N theta)')
parser.add_argument('--no-monic',dest='make_monic',action='store_false',help='Do not adjust the coefficients')
parser.add_argument('-v','--verbose',action='store_true')

parser.add_argument('--homogeneous',dest='outmode',action='store_const',const='homogeneous',help='Print homogeneous coordinates',default='homogeneous')
parser.add_argument('--affine',dest='outmode',action='store_const',const='affine',help='Print affine coordinates')
parser.add_argument('--angles',dest='outmode',action='store_const',const='angles',help='Print angles of ideal vertices (N=2) ')
parser.add_argument('--affine-square',dest='outmode',action='store_const',const='affine-square',help='Print affine coordinates in square normalization (N=3)')

# TODO: Add and re-implement root-images from ngon.py

advgroup = parser.add_argument_group('advanced solver options')
advgroup.add_argument('-m','--method',default='fourier',help='Computation method (fourier or euler)')
advgroup.add_argument('--pde-nmesh',type=int,default=511)
advgroup.add_argument('--pde-thresh',type=float,default=1e-9)
advgroup.add_argument('--pde-maxiter',type=int,default=5000)
advgroup.add_argument('--ode-thresh',type=float,default=1e-10)
advgroup.add_argument('--ode-rstep',type=float,default=1e-4)
advgroup.add_argument('--radius',type=float,help='Inradius of square in which to solve the self-duality equation (default: automatic selection)',default=None)
advgroup.add_argument('--ode-nsteps',type=int,default=500000)


args = parser.parse_args()

if args.outmode=='affine-square' and args.N==2:
    raise ValueError('--affine-square only valid for N=3')

if args.outmode=='angles' and args.N==3:
    raise ValueError('--angles only valid for N=2')

# General imports
import logging
if args.verbose:
    from logconfig import logconfig
    logconfig(filename=None)
logger = logging.getLogger(__name__)
import numpy as np
from planarode import planarODE
from extsde import ExtendedSelfDualityEquation
from polynomial_utils import approx_best_rmax


# CONSTRUCT THE POLYNOMIAL

roots = None
if args.roots:
    # convert roots to polynomial
    roots = args.coefs
    from operator import mul
    from itertools import combinations
    if args.zero:
        kmin = 1
        roots.append(-1.0*sum(roots))
    else:
        kmin = 0
    degree = len(roots)
    coefs = []
    for k in range(degree,kmin,-1):
        c = 0
        for subset in combinations(roots,k):
            c += reduce(mul, subset)
        if k%2 == 1:
            c = -c
        coefs.append(c)
    if args.zero:
        coefs.append(0.0)
    coefs.append(1.0)
else:
    # use coefficients as given
    coefs = args.coefs
    degree = len(coefs)-1

# Apply angle
coefs = [ np.exp(-1j*args.N*args.theta)*a for a in coefs]


# COMPUTE THE SELF DUAL METRIC

logger.info('Computing self dual metric for N={}, coefs={}'.format(args.N,str(coefs)))

if args.radius == None:
    args.radius = approx_best_rmax(args.N,coefs,args.pde_nmesh)
    logger.info('Using radius %f\n' % args.radius)

metric = ExtendedSelfDualityEquation(args.N, coefs, args.radius, args.pde_nmesh, args.pde_thresh, args.pde_maxiter, args.method)


# COMPUTE THE VERTICES (HOMOGENEOUS COORDINATES)
# mostly copied from framedata.framedat.computeverticesfromODE

ODE = planarODE.buildhiggsODE(metric)
nvert = degree+args.N
# Initial angle for integration rays
star_angle = -(np.angle(coefs[-1])-np.pi)/nvert
frames = ODE.integrate_rays(n=nvert,theta0=star_angle,r=args.radius,step=args.ode_rstep,tol=args.ode_thresh,nsteps=args.ode_nsteps,return_type='frame')
vertices = []
for k,F in enumerate(frames):
    w,v = np.linalg.eig(F)
    i = np.argmin(abs(w))
    if abs(w[i]) < 10 and abs(w[i])>0.10:
        logger.warning("Eigenvalue %s too close to 1: full eigenvalue list %s" % (w[i],w))
    subdom = v[:,i]
    imagfrac = np.linalg.norm(np.imag(subdom)) / np.linalg.norm(subdom)
    if imagfrac > 1e-6:
        raise ValueError('Vertex {} has unexpected imaginary part; v={}'.format(k,str(subdom)))
    subdom = np.real(subdom)
    vertices.append(subdom)

if args.outmode == 'homogeneous':
    for v in vertices:
        print(','.join(str(x) for x in v))
elif args.outmode == 'affine':
    for v in vertices:
        a = v[:-1] / v[-1]
        print(','.join(str(x) for x in a))
elif args.outmode == 'affine-square':
    raise NotImplementedError('Affine square not implemented yet')
elif args.outmode == 'angles':
    for v in vertices:
        print(2.0*np.arctan(v[0]/v[1]))
        
