'''Make dual-pane matplotlib figure showing spectral coordinates and differences (opers)'''

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


PALETTE = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']

def prettify_mpl():
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    rc('text', usetex=True)
    rc('figure',figsize=(15,6))
    rc('axes',titlesize=15)
    rc('axes',labelsize=15)
    rc('xtick',labelsize=15)
    rc('ytick',labelsize=15)
    rc('legend',fontsize=13)

def texname(s,theta=None,func=None):
    if len(s) == 2:
        n = '%s_%s' % (s[0],s[1])
        if func == 'abs':
            return '$|'+n+'|$'
        elif func == 'arg':
            return '$\\mathrm{arg}('+n+')$'
        else:
            return '$'+n+'$'
    elif len(s) == 4:
        out='$(%s_%s,%s_%s)$ opers' % (s[0],s[1],s[2],s[3])
        if theta!=None:
            out += ' at $\\vartheta=%g$' % theta
        return out
    else:
        raise NotImplementedError

def sign_convert(T):
    '''Convert a cluster list from the code convention (all Hitchin values negative) to the paper convention (all Hitchin values positive)'''
    return [-x for x in T]

def get_comparison_data(th,theta,hlims=(None,None)):
    logger.info('Finding comparison data for theory {} opers, theta={}'.format(th['name'],theta))
    comparison_data = []
    for h in th['abshlist']:
        if hlims[0]!=None and h < hlims[0]:
            continue
        if hlims[1]!=None and h > hlims[1]:
            continue
        
        try:
            d = comparisons.compareClusters(th['name'],oper=True,absh=h,theta=theta)
        except Exception as e:
            logger.warning('No data found for h={} (error was: {})'.format(h,str(e)))
            continue

        for k in d:
            if k.endswith('cluster'):
                d[k] = sign_convert(d[k])

        d['hinv'] = 1.0/h
        zeta = h*np.exp(theta*1j)
        d['logsfcluster'] = [ Z/zeta for Z in th['centralcharges'] ]
        comparison_data.append(d)
        
    return comparison_data

# Color lightening function by Ian Hincks (https://gist.github.com/ihincks/6a420b599f43fcd7dbd79d56798c4e5a)
def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])

def single_or_double(x,i):
    if type(x)==str:
        return x
    else:
        return x[i]

def guided_arg(z,theta):
    '''Among choices for arg(z), return the one closest to theta'''
    p = np.angle(z)
    delta = theta-p
    return p + (np.rint(delta/(2*np.pi)))*2*np.pi 
    
def spectral_coord_dual_plot(th,theta,hlims=(None,None),diffylims=(None,None),xlog=False,ylog=True,showwkb=[],exclude=[],legend=False,legend_ncol=None,title_suffix='',extra_left_yrange=None,tex_theoryname=None,diffmode='rel'):
    '''Make a two-pane plot: Spectral coordinates (left) and relative differences (right)'''

    cd = get_comparison_data(th,theta,hlims)

    if diffmode == 'rel':
        diffstr = 'Relative'
        diffkey = 'reldiff'
        errkey = 'relode'
    else:
        diffstr = 'Absolute'
        diffkey = 'absdiff'
        errkey = 'absode'

    wkbs = []
    for v in showwkb:
        if v == 'all':
            wkbs = th['Xclusternames']
        elif v == 'none':
            wkbs = []
        else:
            if v in th['Xclusternames']:
                wkbs.append(v)
            else:
                raise ValueError('Unknown cluster {}'.format(v))


    prettify_mpl()
    fig = plt.figure()
    grid = plt.GridSpec(2,2)
    grid.update(hspace=0)

    #--------------------------------------------------------------------------
    # arg(clusters) - Left bottom
    #--------------------------------------------------------------------------

    ax = fig.add_subplot(grid[1,0])
    if tex_theoryname == None:
        tex_theoryname = texname(th['name'],theta)

    if xlog:
        ax.set_xscale('log')
    ax.set_xlabel('$|\\hbar|^{-1}$')
    ax.set_ylabel('$\\mathrm{arg}(X_i)$')

    hinvs = [ x['hinv'] for x in cd ]

    for i,cname in enumerate(th['Xclusternames']):
        if cname in exclude:
            # Suppress this cluster
            continue
        
        c = PALETTE[i]
        prefix = texname(cname,func='arg')
        if cname in wkbs:
            ax.plot(hinvs, [ x['logsfcluster'][i].imag for x in cd ], color=lighten_color(c),label=prefix+' WKB approx',ls='--',zorder=10)

        ax.plot(hinvs, [ guided_arg(x['fdcluster'][i],x['logsfcluster'][i].imag) for x in cd ],marker='+',linestyle='none',color=c,label=prefix+' PDE',zorder=100,markersize=8)

        ax.plot(hinvs, [ guided_arg(x['xarcluster'][i],x['logsfcluster'][i].imag) for x in cd ],marker='x',linestyle='none',color=c,label=prefix+' IEQ',zorder=100,markersize=6)

    if extra_left_yrange != None:
        y0,y1 = ax.get_ylim()
        y1 = y0 + (1.0 + extra_left_yrange)*(y1-y0)
        ax.set_ylim(y0,y1)
        
    if legend:
        if legend_ncol == None:
            ax.legend()
        else:
            ax.legend(ncol=legend_ncol)
        
    x0,x1 = ax.get_xlim()


    #--------------------------------------------------------------------------
    # abs(clusters) - Left top
    #--------------------------------------------------------------------------

    ax = fig.add_subplot(grid[0,0],sharex=ax)
    if tex_theoryname == None:
        tex_theoryname = texname(th['name'],theta)

    # Suppress all tick labels on shared x axis
    plt.setp(ax.get_xticklabels(which='both'), visible=False)

    ax.set_title('{tn}: Results {suffix}'.format(tn=tex_theoryname,suffix=single_or_double(title_suffix,0)))
    if xlog:
        ax.set_xscale('log')
    if ylog:
        ax.set_yscale('log')

    ax.set_ylabel('$|X_i|$')
    
    for i,cname in enumerate(th['Xclusternames']):
        if cname in exclude:
            # Suppress this cluster
            continue
        
        c = PALETTE[i]
        prefix = texname(cname,func='abs')
        if cname in wkbs:
            ax.plot(hinvs, [ np.exp(x['logsfcluster'][i].real) for x in cd ], color=lighten_color(c),label=prefix+' WKB approx',ls='--',zorder=10)

        ax.plot(hinvs, [ abs(x['fdcluster'][i]) for x in cd ],marker='+',linestyle='none',color=c,label=prefix+' PDE',zorder=100,markersize=8)

        ax.plot(hinvs, [ abs(x['xarcluster'][i]) for x in cd ],marker='x',linestyle='none',color=c,label=prefix+' IEQ',zorder=100,markersize=6)

    if extra_left_yrange != None:
        y0,y1 = ax.get_ylim()
        y1 = np.exp( np.log(y0) + (1.0 + extra_left_yrange)*(np.log(y1)-np.log(y0)) )
        ax.set_ylim(y0,y1)
        
    if legend:
        if legend_ncol == None:
            ax.legend()
        else:
            ax.legend(ncol=legend_ncol)

    x0,x1 = ax.get_xlim()


    #--------------------------------------------------------------------------
    # differences - Right
    #--------------------------------------------------------------------------

    ax = fig.add_subplot(grid[:,1])
    ax.set_title('{tn}: {diffstr} difference{suffix}'.format(tn=tex_theoryname,diffstr=diffstr,suffix=single_or_double(title_suffix,1)))
    if xlog:
        ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(*diffylims)
    
    ax.set_xlim(x0,x1)

    ax.set_xlabel('$|\\hbar|^{-1}$')

    for i in range(len(th['Xclusternames'])):
        if th['Xclusternames'][i] in exclude:
            # Suppress this cluster
            continue
        c = PALETTE[i]
        ax.plot(hinvs, [ x[diffkey][i] for x in cd ],marker='o',linestyle='none',color=c,label=texname(th['Xclusternames'][i])+' rel. err.',markersize=4)
        ax.plot(hinvs, [ x['errest'][errkey][i] for x in cd ],marker='3',linestyle='none',color=c,label=texname(th['Xclusternames'][i])+' ODE rel err estimate',markersize=8) #4
    
    if legend:
        ax.legend()

    fig.tight_layout()

    return fig
                          

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('theoryname')
    parser.add_argument('theta',type=float)
    
    parser.add_argument('--exclude',default=[],nargs='+',help='Names of cluster variables to omit (space separated)')

    parser.add_argument('--xlog',
                    action='store_const',
                    default = None,
                    const = True,
                    dest='xlog')
    parser.add_argument('--xlinear',
                    action='store_const',
                    const = False,
                    dest='xlog')
    parser.add_argument('--ylog',
                    action='store_const',
                    default = None,
                    const = True,
                    dest='ylog')
    parser.add_argument('--ylinear',
                    action='store_const',
                    const = False,
                    dest='ylog')
    parser.add_argument('--absdiff',action='store_const',const='abs',dest='diffmode',help='Show absolute differences in right pane')
    parser.add_argument('--reldiff',action='store_const',const='rel',dest='diffmode',help='Show relative differences in right pane')
    parser.add_argument('--show-wkb',default=[],nargs='+',dest='showwkb',help='Which clusters to show WKB approximation ("all", "none", or space-separated cluster names)')
    parser.add_argument('--hmin',default=None,type=float)
    parser.add_argument('--hmax',default=None,type=float)
    parser.add_argument('--diffymin',default=None,type=float,help='y axis lower limit for difference plot (right pane)')
    parser.add_argument('--diffymax',default=None,type=float,help='y axis upper limit for difference plot (right pane)')
    parser.add_argument('--title-suffix',default=None,help='string to append to titles')
    parser.add_argument('--legend',action='store_true')
    parser.add_argument('--legend-ncol',default=None,type=int)
    parser.add_argument('--extra-left-yrange',type=float,default=None,help='fraction of original y range in left plot (coordinates) to add to top, allowing more space for the legend')
    parser.add_argument('--tex-theoryname',default=None,help='TeX version of theory name for display in title')
    parser.add_argument('-v','--verbose',action='store_true')
    parser.add_argument('-o','--output',help='output filename')

    parser.set_defaults(diffmode='rel')

    args = parser.parse_args()

    if args.verbose:
        logconfig.loglevel(logging.INFO)

    if not args.output:
        args.output = 'paperfig' + args.theoryname + '-oper.pdf'


    th = theorydata.getdata(args.theoryname)
    extra_kwargs = dict()
    for k in ['xlog','ylog','title_suffix','legend','legend_ncol','extra_left_yrange','tex_theoryname','diffmode']:
        if getattr(args,k) != None:
            extra_kwargs[k]=getattr(args,k)
    f = spectral_coord_dual_plot(th,args.theta,hlims=(args.hmin,args.hmax),diffylims=(args.diffymin,args.diffymax),showwkb=args.showwkb,exclude=args.exclude,**extra_kwargs)
    f.savefig(args.output)
