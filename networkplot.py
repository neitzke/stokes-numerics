'''Python interface to Neitzke's swn-plotter'''
from __future__ import absolute_import
from __future__ import print_function
import os
import sys
import argparse
import time
import subprocess
import uuid
import json
import tempfile

try:
    from subprocess import DEVNULL
except ImportError:
    DEVNULL = open(os.devnull, 'wb')

scriptdir=os.path.dirname(os.path.realpath(__file__))
with open(os.path.join(scriptdir,"swnpipe.m"),"rt") as templatefile:
    _template = templatefile.read()

def cplxstr(z):
    x = z.real
    y = z.imag
    return '(%.4f + I*%.4f)' % (x,y)

def polystr(coefs):
    s = ''
    for c in reversed(coefs):
        if s:
             s = cplxstr(c) + ' + z*( ' + s + ')'
        else:
            s = cplxstr(c)
    return s

class plot:
    def __init__(self,cutoff=100,rmax=3,rstep=0.05,rank=2,theta=0.0,coefs=[-1,0,0,1],secondary_coefs=None,verbose=False,xvfb=False):
        self.outfn = os.path.join(tempfile.gettempdir(),uuid.uuid4().hex + '.json')
        self.cutoff = cutoff
        self.theta = theta
        self.rank = rank
        self.coefs = coefs
        self.secondary_coefs = secondary_coefs
        self.rmax = rmax
        self.rstep = rstep

        if rank==2 or not secondary_coefs or not any(secondary_coefs):
            # The swn plotter treats exact zero and numerical zero differently
            # so it needs to be forced to use "exact zero" behavior if the user
            # specified a secondary polynomial but all coefs are zero
            secondary_polystr = '0'
        else:
            secondary_polystr = polystr(self.secondary_coefs)
        
        self.repvars = {
            '$CUTOFF$': self.cutoff,
            '$ANGLE$' : self.theta,
            '$RANK$' : self.rank,
            '$OUTFN$' : self.outfn,
            '$POLY$' : polystr(self.coefs),
            '$SECONDARYPOLY$' : secondary_polystr,
            '$RMAX$' : self.rmax,
            '$RSTEP$' : self.rstep,
        }
        self.script = _template
        for k in self.repvars:
            self.script = self.script.replace(k, str(self.repvars[k]))

        cmdlst = ['math']
        if xvfb:
            cmdlst = ['xvfb-run'] + cmdlst
        self.p = subprocess.Popen(cmdlst,stdin=subprocess.PIPE,stdout=DEVNULL,
                                  stderr=subprocess.STDOUT)
        self.p.communicate(self.script.encode('ascii'))
    
        if self.p.returncode != 0:
            raise ChildProcessError('Nonzero return code from mathematica')

        if not os.path.exists(self.outfn):
            raise FileNotFoundError('JSON output from mathematica call not found')
        try:
            with open(self.outfn,'rt') as infile:
                self.raw_output = json.load(infile)
            self.turningpoints = [ complex(*t) for t in self.raw_output['turningpoints'] ]
            self.trajectories = [ [ complex(*t) for t in x ] for x in self.raw_output['trajectories'] ]
        except json.decoder.JSONDecodeError:
            pass
        finally:
                os.unlink(self.outfn)


def main():
    '''Compute vertices of polygon (hyperbolic ideal or real projective convex) corresponding to a meromorphic cyclic Higgs bundle'''

    import sys
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('-N',type=int,default=2,help='Rank of Higgs bundle')
    parser.add_argument('coefs',type=complex,nargs='+',help='List of polynomial coefficients, constant term first, complex() format, with N-differential first, followed by (N-1)-differential')
    parser.add_argument('-d','--degree0',type=int,default=None,help='Degree of N-differential; used to split coefs into N and N-1 differentials.  Valid only for rank 3.')
    parser.add_argument('-t','--theta',type=float,default=0.0)
    parser.add_argument('-R',type=float,default=1.0)
    parser.add_argument('--cutoff',type=float,default=50.0,help='network tracing cutoff')
    parser.add_argument('-v','--verbose',action='store_true')
    parser.add_argument('--network-radius',type=float,default=3,help='radius of network view')
    parser.add_argument('--network-rstep',type=float,default=0.05)
    parser.add_argument('-o','--output',help='output filename')
    parser.add_argument('--color-vertices',action='store_true',help='color the firs two sectors')

    args = parser.parse_args()

    # General imports
    import logging
    if args.verbose:
        from logconfig import logconfig
        logconfig(filename=None)
    logger = logging.getLogger(__name__)

    import numpy as np
    from latexutil import pretty_polynomial
    import os
    import time
    from networkplot import plot

    if not args.output:
        args.output = 'network%d.pdf' % int(time.time())

    def tempnam(prefix='',suffix=''):
        import uuid
        return os.path.join('/tmp',prefix + uuid.uuid4().hex + suffix)


    # CONSTRUCT THE POLYNOMIAL

    if args.degree0 != None:
        if args.N == 2:
            raise ValueError('degree0 can only be specified if rank is 3')
        coefs0 = args.coefs[:args.degree0+1]
        coefs1 = args.coefs[args.degree0+1:]
        logger.info('Both differentials given:')
        logger.info('\tcubic:  \t' + str(coefs0))
        logger.info('\tquadratic:\t' + str(coefs1))
    else:
        coefs0 = args.coefs
        coefs1 = None
    
    degree = len(coefs0)-1

    scaled_coefs0 = [ (args.R*np.exp(-1j*args.theta))**(args.N) * x for x in coefs0 ]
    if coefs1 != None:
        scaled_coefs1 = [ (args.R*np.exp(-1j*args.theta))**(args.N-1) * x for x in coefs1 ]
    else:
        scaled_coefs1 = None

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
    rc('figure',figsize=(9,9))
    rc('figure',titlesize=30)
    rc('axes',titlesize=15)
    rc('xtick',labelsize=15)
    rc('ytick',labelsize=15)

    fig = plt.figure() #figsize=(16,9),dpi=120)

    # Axis for network plot
    nax = plt.subplot('111')
    nax.axis('off')
    nax.axis('off')
    nax.set_xticks([])
    nax.set_yticks([])
    nax.set_aspect(1)
    nax.set_ylim(-args.network_radius,args.network_radius)
    nax.set_xlim(-args.network_radius,args.network_radius)

    # NETWORK PLOTTING

    logger.info('Drawing spectral network')
    
    network = None
    try:
        network = plot(cutoff=args.cutoff, rmax=args.network_radius, rstep=args.network_rstep, rank=args.N, coefs=scaled_coefs0, secondary_coefs=scaled_coefs1)
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

    logger.info('Layout and typesetting')
    plt.tight_layout()
    plt.savefig(args.output)

if __name__=='__main__':
    main()
