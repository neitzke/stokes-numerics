'''Make matplotlib figures showing heatmaps'''

import sys
import os
# Allow imports from parent directory
sys.path.insert(1, os.path.join(sys.path[0], '..'))

import logging
if __name__ == '__main__':
    import logconfig
    logconfig.logconfig(filename=None)
    logconfig.loglevel(logging.WARN)
logger = logging.getLogger(__name__)
import comparisons
import framedata
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rc
import theorydata
import numpy as np

import hkmetric

PDE_NMESH = 1000

def prettify_mpl():
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    rc('text', usetex=True)
    rc('figure',figsize=(15,6))
    rc('axes',titlesize=12)
    rc('axes',labelsize=12)
    rc('xtick',labelsize=12)
    rc('ytick',labelsize=12)
    rc('legend',fontsize=12)

def heatmap(ax, func, nx = 0, ny = 0, square = True, trimfactor = None, vmin = None, vmax = None, cmap = "hot", extent = None):
    func = func.real
    if square:
        nx = len(func)
        ny = nx
    funcmat = func.reshape((nx, ny))
    if trimfactor is not None:
        bx = int(nx*trimfactor/2.0)
        by = int(ny*trimfactor/2.0)
        funcmat = funcmat[int(nx/2)-bx:int(nx/2)+bx,int(ny/2)-by:int(ny/2)+by]
    if vmax is None:
        plot = ax.imshow(funcmat, cmap = cmap, interpolation = "nearest", extent = extent)
    else:
        plot = ax.imshow(funcmat, cmap = cmap, interpolation = "nearest", vmin = vmin, vmax = vmax, extent = extent)

    return plot

def heatmap_triple_plot(pde_nmesh = PDE_NMESH, trimfactor = 0.2, rmax = 10.0, Lambda = 0, clist = [0.3,1,2.5], c = 1, useclist = False):
    prettify_mpl()
    fig = plt.figure(figsize=(10,8.6))

    if not useclist:
        clist = [c]

    nrows = len(clist)

    grid = plt.GridSpec(nrows=nrows,ncols=3)
    # grid.update(hspace=5)

    for n,c in enumerate(clist):
        ax0 = fig.add_subplot(grid[n,0])
        ax1 = fig.add_subplot(grid[n,1])
        ax2 = fig.add_subplot(grid[n,2])

        o = hkmetric.fdcomputeG(pde_nmesh = pde_nmesh, rmax = rmax, Lambda = Lambda, c = c)
        extent = [-rmax*trimfactor, rmax*trimfactor, -rmax*trimfactor, rmax*trimfactor]

        ax0.set_title("${\\mathcal I}$ for $c = %2.1f$" % c)
        plot0 = heatmap(ax0, o["integrand"], trimfactor = trimfactor, vmax = 10, vmin = 0, cmap = "binary", extent = extent)
        fig.colorbar(plot0, ax=ax0, fraction=0.046, pad=0.04)

        ax1.set_title("${\\mathcal I}^{\\mathrm{sf}}$ for $c = %2.1f$" % c)
        plot1 = heatmap(ax1, o["integrandsf"], trimfactor = trimfactor, vmax = 10, vmin = 0, cmap = "binary", extent = extent)
        fig.colorbar(plot1, ax=ax1, fraction=0.046, pad=0.04)

        ax2.set_title("${\\mathcal I}-{\\mathcal I}^{\\mathrm{sf}}$ for $c = %2.1f$" % c)
        plot2 = heatmap(ax2, o["integrand"]-o["integrandsf"], trimfactor = trimfactor, vmax = 4, vmin = -4, cmap = "seismic", extent = extent)
        fig.colorbar(plot2, ax=ax2, fraction=0.046, pad=0.04)

    grid.tight_layout(fig)

    return fig

def texname(s,theta=None):
    if len(s) == 2:
        return '$%s_%s$' % (s[0],s[1])
    elif len(s) == 4:
        out='$(%s_%s,%s_%s)$' % (s[0],s[1],s[2],s[3])
        if theta!=None:
            out += r' at $\vartheta=%g$' % theta
        return out
    else:
        raise NotImplementedError

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('--pde-nmesh',default=None,type=int)
    parser.add_argument('--trimfactor',default=None,type=float)
    parser.add_argument('--rmax',default=None,type=float)
    parser.add_argument('--Lambda',default=None,type=float)
    parser.add_argument('--c',default=None,type=float)
    parser.add_argument('-m','--multiplerows',action='store_true')
    parser.add_argument('-v','--verbose',action='store_true')
    parser.add_argument('-o','--output',help='output filename')

    args = parser.parse_args()

    if args.verbose:
        logconfig.loglevel(logging.INFO)

    if not args.output:
        args.output = 'paperfigheatmaps.pdf'

    extra_kwargs = dict()
    for k in ['pde_nmesh', 'trimfactor', 'rmax', 'Lambda', 'c']:
        if getattr(args,k) != None:
            extra_kwargs[k]=getattr(args,k)
    f = heatmap_triple_plot(useclist = args.multiplerows, **extra_kwargs)
    f.savefig(args.output)
