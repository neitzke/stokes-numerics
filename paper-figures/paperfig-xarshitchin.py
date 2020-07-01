'''Make matplotlib figures showing xars'''

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

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rc

import numpy as np
import integralequations
import integralequationsplotting


def prettify_mpl():
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    rc('text', usetex=True)
    rc('figure',figsize=(15,6))
    rc('axes',titlesize=15)
    rc('axes',labelsize=15)
    rc('xtick',labelsize=21)
    rc('ytick',labelsize=21)
    rc('legend',fontsize=13)

from math import log,exp
def xardata(Rlist=None):
    if Rlist is None:
        Rlist = [exp(-t/2.0) for t in range(20)]
    xarlist = []
    datalist = []
    for R in Rlist:
        xar = integralequations.computeXar(theoryname = "A1A2", R = R, failonmaxiter = True)
        xinstlist = xar.getRaydatum(0).xinstlists[0]
        tlist = xar.getRaydatum(0).tlist
        datalist.append([tlist,xinstlist])
        xarlist.append(xar)
    return datalist,xarlist,Rlist

def xarplothitchin():
    prettify_mpl()

    datalist,xarlist,Rlist = xardata()

    fig = integralequationsplotting.plotData(datalist, progressive=False, color="rainbow", tcutoff=20)

    return fig

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)

    # parser.add_argument('--pde-nmesh',default=None,type=int)
    # parser.add_argument('--trimfactor',default=None,type=float)
    # parser.add_argument('--rmax',default=None,type=float)
    # parser.add_argument('--Lambda',default=None,type=float)
    # parser.add_argument('--c',default=None,type=float)
    # parser.add_argument('-m','--multiplerows',action='store_true')
    parser.add_argument('-v','--verbose',action='store_true')
    parser.add_argument('-o','--output',help='output filename')

    args = parser.parse_args()

    if args.verbose:
        logconfig.loglevel(logging.INFO)

    if not args.output:
        args.output = 'paperfigxarshitchin.pdf'

    extra_kwargs = dict()
    # for k in ['pde_nmesh', 'trimfactor', 'rmax', 'Lambda', 'c']:
    #     if getattr(args,k) != None:
    #         extra_kwargs[k]=getattr(args,k)
    f = xarplothitchin()
    f.savefig(args.output)
