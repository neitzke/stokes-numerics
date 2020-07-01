'''Make matplotlib figures showing hyperkahler metric data'''

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
import jsonpickle

import comparisons
import framedata
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rc
import theorydata
import numpy as np

import namegen
import hkmetric

from math import exp

def prettify_mpl():
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    rc('text', usetex=True)
    rc('figure',figsize=(15,6))
    rc('axes',titlesize=20)
    rc('axes',labelsize=20)
    rc('xtick',labelsize=20)
    rc('ytick',labelsize=20)
    rc('legend',fontsize=20)

def loadmetricdata(Lambda = 0.0):
    clist = [0.03*exp(0.2*n) for n in range(25)]
    data = []
    for c in clist:
        ieqfilename = namegen.ieq_metric_filename(fullpath=True, c=c, Lambda=Lambda)
        ieqmetric = hkmetric.ieqmetricdata.load(ieqfilename)
        fdfilename = namegen.fd_metric_filename(fullpath=True, c=c, Lambda=Lambda)
        fdmetric = hkmetric.fdmetricdata.load(fdfilename)
        # Gsf is fast: just compute it on the fly
        Gsf = hkmetric.ieqcomputeG(c = c, R = 1, Lambda = Lambda, eps = 1e-6, leadingapprox = True, steps = 128)["G"]
        data.append( {"c": c, "Gfd": fdmetric.G, "Gieq": ieqmetric.G, "Gsf": Gsf} )
    return data

def hkmetric_plot(data):
    prettify_mpl()
    fig = plt.figure()
    grid = plt.GridSpec(1,2)
    grid.update(hspace=5)

    ax0 = fig.add_subplot(grid[0,0])
    ax1 = fig.add_subplot(grid[0,1])

    xtable = [i["c"] for i in data]
    ytable1 = [i["Gfd"] for i in data]
    ax0.plot(xtable, ytable1, linestyle = "none", color = "blue", marker = "+", label = "$g^{\mathrm{DE}}$")
    ytable2 = [i["Gieq"] for i in data]
    ax0.plot(xtable, ytable2, linestyle = "none", color = "blue", marker = "x", label = "$g^{\mathrm{IEQ}}$")
    ytable3 = [i["Gsf"] for i in data]
    ax0.plot(xtable, ytable3, linestyle = "dashed", color = "orange", marker = " ", label = "$g^{\mathrm{sf}}$")
    ax0.legend(loc="best", fancybox=True)
    ax0.set_xlabel("$c$")
    ax0.set_xlim(left=0)
    ax0.set_ylim(bottom=0)

    xtable = [i["c"] for i in data]
    ytable = [i["Gieq"]-i["Gfd"] for i in data]
    ax1.plot(xtable, ytable, linestyle = "none", color = "red", marker = "o", label = "$g^{\mathrm{DE}}-g^{\mathrm{IEQ}}$")
    ax1.legend(loc="best", fancybox=True)
    ax1.set_xlabel("$c$")
    ax1.set_xlim(left=0)
    ax1.set_ylim(bottom=0)

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

    # parser.add_argument('--pde-nmesh',default=None,type=int)
    # parser.add_argument('--trimfactor',default=None,type=float)
    parser.add_argument('--Lambda',default=0.0,type=float)
    parser.add_argument('-v','--verbose',action='store_true')
    parser.add_argument('-o','--output',help='output filename')

    args = parser.parse_args()

    if args.verbose:
        logconfig.loglevel(logging.INFO)

    if not args.output:
        args.output = 'paperfighkmetric.pdf'

    extra_kwargs = dict()
    for k in []:
        if getattr(args,k) != None:
            extra_kwargs[k]=getattr(args,k)

    f = hkmetric_plot(data=loadmetricdata(Lambda = args.Lambda), **extra_kwargs)

    f.savefig(args.output)
