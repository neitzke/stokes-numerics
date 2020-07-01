'''Compute vertices of polygon (hyperbolic ideal or real projective convex) corresponding to a meromorphic cyclic Higgs bundle'''

# Requires:
#    mathematica

# Imports needed to process arguments
from __future__ import absolute_import
from __future__ import print_function
import sys
import argparse

# This was adapted from ngon.py from Blaschke

parser = argparse.ArgumentParser(description=__doc__)

parser.add_argument('-N',type=int,default=2,help='Rank of Higgs bundle')
parser.add_argument('coefs',type=complex,nargs='+',help='If using --coefs, list of polynomial coefficients, constant term first, complex() format, defining the holomorphic N-differential')
parser.add_argument('-r','--roots',action='store_true',help='instead of coefficients, roots are given')
parser.add_argument('-t','--theta',type=float,default=0.0)
parser.add_argument('-R',type=float,default=1.0)
parser.add_argument('--cutoff',type=float,default=50.0,help='network tracing cutoff')
#parser.add_argument('--no-monic',dest='make_monic',action='store_false',help='Do not adjust the coefficients')
parser.add_argument('-v','--verbose',action='store_true')
parser.add_argument('--polystring',default=None)
parser.add_argument('--resolution',type=int,default=512)
parser.add_argument('--radius',type=float,default=1.1,help='radius of polygon view')
parser.add_argument('--network-radius',type=float,default=3,help='radius of network view')
parser.add_argument('--network-rstep',type=float,default=0.01)
parser.add_argument('-o','--output',help='output filename')
parser.add_argument('--show-theta',action='store_true',help='show theta in title even if it is zero',default=None)
parser.add_argument('--R-fmt',default='%1.2f')
parser.add_argument('--theta-fmt',default='%1.2f')
parser.add_argument('--color-vertices',action='store_true',help='color vertex corresponding to first two rays red,green')

# TODO: Add and re-implement root-images from ngon.py

advgroup = parser.add_argument_group('advanced solver options')
advgroup.add_argument('-m','--method',default='fourier',help='Computation method (fourier or euler)')
advgroup.add_argument('--pde-nmesh',type=int,default=511)
advgroup.add_argument('--pde-thresh',type=float,default=1e-9)
advgroup.add_argument('--pde-maxiter',type=int,default=5000)
advgroup.add_argument('--ode-thresh',type=float,default=1e-10)
advgroup.add_argument('--ode-rstep',type=float,default=1e-4)
advgroup.add_argument('--solver-radius',type=float,help='Inradius of square in which to solve the self-duality equation (default: automatic selection)',default=None)
advgroup.add_argument('--ode-nsteps',type=int,default=500000)


args = parser.parse_args()

# General imports
import logging
if args.verbose:
    from logconfig import logconfig
    logconfig(filename=None)
logger = logging.getLogger(__name__)
import numpy as np
from planarode import planarODE
from extsde import ExtendedSelfDualityEquation
from polynomial_utils import approx_best_rmax, equivalent_quasimonic_centered_cyclic
from latexutil import pretty_polynomial
import os
import time
from networkplot import plot
from subprocess import check_call

try:
    from subprocess import DEVNULL
except ImportError:
    import os
    DEVNULL = open(os.devnull, 'wb')

if not args.output:
    args.output = 'framevis%d.pdf' % int(time.time())

if args.output.endswith('.png'):
    check_call(['convert','--version'],stdout=DEVNULL)

if args.show_theta == None:
    args.show_theta = (args.theta != 0)
    
def tempnam(prefix='',suffix=''):
    import uuid
    return os.path.join('/tmp',prefix + uuid.uuid4().hex + suffix)


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
    coefs0 = []
    for k in range(degree,kmin,-1):
        c = 0
        for subset in combinations(roots,k):
            c += reduce(mul, subset)
        if k%2 == 1:
            c = -c
        coefs0.append(c)
    if args.zero:
        coefs0.append(0.0)
    coefs0.append(1.0)
else:
    # use coefficients as given
    coefs0 = args.coefs
    degree = len(coefs0)-1

coefs = [ (args.R*np.exp(-1j*args.theta))**(args.N) * x for x in coefs0 ]

if args.R > 0:
    mcoefs = equivalent_quasimonic_centered_cyclic(args.N,coefs)
else:
    mcoefs = [ 0.0 ] * degree + [ 1.0 ]

if args.polystring == None:
    if any( a.imag != 0 for a in coefs0 ):
        raise ValueError('Can\'t render polynomial string for complex coefs; pass --polystring')
    args.polystring = pretty_polynomial([x.real for x in coefs0])

nvert = degree+args.N
# Initial angle for integration rays
star_angle = -(np.angle(mcoefs[-1])-np.pi)/nvert



# MATPLOTLIB SETUP

import hypideal
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
rc('figure',figsize=(16,9))
rc('figure',titlesize=30)
rc('axes',titlesize=15)
rc('xtick',labelsize=15)
rc('ytick',labelsize=15)

fig = plt.figure() #figsize=(16,9),dpi=120)
if args.show_theta:
    fig.suptitle('$\\phi_{N} = R^{N} e^{{-{N} i \\theta}} ({pstr}) \\;\\;\\;\\; R = \\texttt{{{R}}} \\;\\;\\;\\; \\theta = \\texttt{{{theta}}}$'.format(N=args.N,pstr=args.polystring,R=(args.R_fmt % args.R),theta=(args.theta_fmt % args.theta)))
else:
    fig.suptitle('$\\phi_{N} = R^{N} ({pstr}) \\;\\;\\;\\; R = \\texttt{{{R}}}$'.format(N=args.N,pstr=args.polystring,R=(args.R_fmt % args.R)))

# Axis for network plot
nax = plt.subplot('121')
nax.axis('off')
nax.axis('off')
nax.set_xticks([])
nax.set_yticks([])
nax.set_aspect(1)
nax.set_ylim(-args.network_radius,args.network_radius)
nax.set_xlim(-args.network_radius,args.network_radius)

# Axis for polygon plot
pax = plt.subplot('122')
pax.axis('off')
pax.set_xticks([])
pax.set_yticks([])
pax.set_ylim(-args.radius,args.radius)
pax.set_xlim(-args.radius,args.radius)
pax.set_aspect(1)



# NETWORK PLOTTING

logger.info('Drawing spectral network')

network = None
try:
    network = plot(cutoff=args.cutoff, rmax=args.network_radius, rstep=args.network_rstep, rank=args.N, coefs=mcoefs)
except Exception as e:
    logger.error('Unable to compute spectral network (mathematica missing?): {}'.format(e))

if network:
    for T in network.trajectories:
        nax.add_patch(patches.Polygon([ [ z.real, z.imag ] for z in T ], closed=False, facecolor='none', linewidth=2.5, edgecolor='black',zorder=10))
    nax.scatter( [z.real for z in network.turningpoints],
                 [z.imag for z in network.turningpoints],
                 marker='X',
                 s=10**2,
                 color='orange',
                 zorder=100)

if args.color_vertices:
    t0 = star_angle
    t1 = star_angle + 2*np.pi/nvert
    t2 = star_angle + 4*np.pi/nvert
    thw = 0.5 * np.pi/nvert
    r = args.network_radius*1.95
    degrees = 180.0 / np.pi
    nax.add_patch(patches.Arc([0,0],r,r,
                              theta1=degrees*(t0 - thw),
                              theta2=degrees*(t0 + thw),
                              linewidth=3,
                              color='red'))
    nax.add_patch(patches.Arc([0,0],r,r,
                              theta1=degrees*(t1 - thw),
                              theta2=degrees*(t1 + thw),
                              linewidth=3,
                              color='green'))
    
# COMPUTE THE SELF DUAL METRIC

logger.info('Computing self dual metric for N={}, coefs={}, mcoefs={}'.format(args.N,str(coefs),str(mcoefs)))

if args.solver_radius == None:
    args.solver_radius = approx_best_rmax(args.N,mcoefs,args.pde_nmesh)
    logger.info('Using radius %f\n' % args.radius)

metric = ExtendedSelfDualityEquation(args.N, mcoefs, args.solver_radius, args.pde_nmesh, args.pde_thresh, args.pde_maxiter, args.method)


# COMPUTE THE VERTICES (HOMOGENEOUS COORDINATES)
# mostly copied from framedata.framedate.computeverticesfromODE

ODE = planarODE.buildhiggsODE(metric)
logger.info('Computing subdominant solutions along %d rays' % nvert)
frames = ODE.integrate_rays(n=nvert,theta0=star_angle,r=args.solver_radius,step=args.ode_rstep,tol=args.ode_thresh,nsteps=args.ode_nsteps,return_type='frame')
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


# DRAW PLOYGON

logger.info('Drawing polygon')

CIRCLE_COLOR='#808080'
VERTEX_COLOR='#0000FF'
POLY_COLOR='#B0B0FF'

if args.N == 2:
    # Draw unit circle
    pax.add_patch(patches.Circle((0,0),radius=1,fc='none',color='gray',lw=2,zorder=1))

    P = hypideal.H2IdealPolygon(np.pi + 2.0*np.arctan(v[0]/v[1]) for v in vertices)
    pax.add_patch(P.patch(fc=POLY_COLOR,lw=0,zorder=5))
    pc = PatchCollection(P.vertex_patches(radius=0.02),facecolor=VERTEX_COLOR,zorder=100)
    pax.add_collection(pc)
    if args.color_vertices:
        V = list(P.vertices())
        pax.add_patch(V[0].patch(radius=0.02, facecolor='red', zorder=110))
        pax.add_patch(V[1].patch(radius=0.02, facecolor='green', zorder=110))
else:
    affines = [ v[1:] / v[0] for v in vertices ]
    pax.add_patch(patches.Polygon(affines,fc=POLY_COLOR,lw=0,zorder=5))
    pc = PatchCollection( [ patches.Circle(a,radius=0.02) for a in affines ],facecolor=VERTEX_COLOR,zorder=100)
    pax.add_collection(pc)
    if args.color_vertices:
        pax.add_patch(patches.Circle(affines[0],radius=0.02, facecolor='red', zorder=110))
        pax.add_patch(patches.Circle(affines[1],radius=0.02, facecolor='green', zorder=110))


logger.info('Layout and typesetting')
plt.tight_layout()
logger.info('Writing "{}"'.format(args.output))
plt.savefig(args.output)
