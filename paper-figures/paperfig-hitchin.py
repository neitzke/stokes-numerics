'''Make dual-pane matplotlib figure showing spectral coordinates and differences'''

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
import richardson

PALETTE = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']

def listlike(L):
    try:
        n = len(L)
        return True
    except TypeError:
        pass
    return False

def prettify_mpl():
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    rc('text', usetex=True)
    rc('figure',figsize=(15,6))
    rc('axes',titlesize=15)
    rc('axes',labelsize=15)
    rc('xtick',labelsize=15)
    rc('ytick',labelsize=15)
    rc('legend',fontsize=13)

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

def sign_convert(T):
    '''Convert a cluster list from the code convention (all Hitchin values negative) to the paper convention (all Hitchin values positive)'''
    return [-x for x in T]

def get_comparison_data(th,theta,pde_nmesh=framedata.PDE_NMESH,rlims=(None,None)):
    logger.info('Finding comparison data for theory {}, theta={}, pde_nmesh={}'.format(th['name'],theta,pde_nmesh))
    comparison_data = []
    for R in th['Rlist']:
        if rlims[0]!=None and R < rlims[0]:
            continue
        if rlims[1]!=None and R > rlims[1]:
            continue
        
        try:
            d = comparisons.compareClusters(th['name'],R=R,theta=theta,pde_nmesh=pde_nmesh,scratch=False)
        except Exception as e:
            logger.warning('No data found for R={} (error was: {})'.format(R,str(e)))
            continue

        for k in d:
            if k.endswith('cluster'):
                d[k] = sign_convert(d[k])
        d['R'] = R
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

def spectral_coord_dual_plot(th,theta,pde_nmesh=framedata.PDE_NMESH,rlims=(None,None),diffylims=(None,None),xlog=False,ylog=True,showsf=True,exclude=[],legend=False,legend_ncol=None,title_suffix='',extra_left_yrange=None,tex_theoryname=None,richardson_exponent_delta=None,diffmode='rel'):
    '''Make a two-pane plot: Spectral coordinates (left) and relative differences (right)'''

    cd = get_comparison_data(th,theta,pde_nmesh,rlims)
    if listlike(pde_nmesh):    
        errorbars = True
    else:        
        errorbars = False

    if diffmode == 'rel':
        diffstr = 'Relative'
        diffkey = 'reldiff'
        errkey = 'relpde'
    else:
        diffstr = 'Absolute'
        diffkey = 'absdiff'
        errkey = 'abspde'

    semiflats = []
    for v in showsf:
        if v == 'all':
            semiflats = th['Xclusternames']
        elif v == 'none':
            semiflats = []
        else:
            if v in th['Xclusternames']:
                semiflats.append(v)
            else:
                raise ValueError('Unknown cluster {}'.format(v))

    prettify_mpl()
    fig = plt.figure()

    #--------------------------------------------------------------------------
    # clusters - Left
    #--------------------------------------------------------------------------

    ax = fig.add_subplot(1,2,1)
    if tex_theoryname == None:
        tex_theoryname = texname(th['name'],theta)
        
    ax.set_title('{tn}: Results {suffix}'.format(tn=tex_theoryname,suffix=single_or_double(title_suffix,0)))
    if xlog:
        ax.set_xscale('log')
    if ylog:
        ax.set_yscale('log')
    ax.set_xlabel('$R$')
    
    Rs = [ x['R'] for x in cd ]

    for i,cname in enumerate(th['Xclusternames']):
        if cname in exclude:
            # Suppress this cluster
            continue
        
        # Indicate in label whether the sign is flipped in the figure
        # (we only display positive values)
        if len(cd) == 0 or cd[-1]['sfcluster'][i]>0:
            prefix=''
        else:
            prefix='-'
        c = PALETTE[i]
        prefix += texname(cname)
        if cname in semiflats:
            ax.plot(Rs, [ abs(x['sfcluster'][i]) for x in cd ], color=lighten_color(c),label=prefix+' sf approx',ls='--',zorder=10)

        ax.plot(Rs, [ abs(x['fdcluster'][i]) for x in cd ],marker='+',linestyle='none',color=c,label=prefix+' PDE',zorder=100,markersize=8)

        ax.plot(Rs, [ abs(x['xarcluster'][i]) for x in cd ],marker='x',linestyle='none',color=c,label=prefix+' IEQ',zorder=100,markersize=6)

    if extra_left_yrange != None:
        y0,y1 = ax.get_ylim()
        if ylog:
            y1 = np.exp( np.log(y0) + (1.0 + extra_left_yrange)*(np.log(y1)-np.log(y0)) )
        else:
            y1 = y0 + (1.0 + extra_left_range)*(y1-y0)
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

    ax = fig.add_subplot(1,2,2)
    ax.set_title('{tn}: {diffstr} difference{suffix}'.format(tn=tex_theoryname,diffstr=diffstr,suffix=single_or_double(title_suffix,1)))

    if xlog:
        ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(*diffylims)
    
    ax.set_xlim(x0,x1)

    ax.set_xlabel('$R$')

    for i,cname in enumerate(th['Xclusternames']):
        if cname in exclude:
            # Suppress this cluster
            continue
        c = PALETTE[i]
        ax.plot(Rs, [ x[diffkey][i] for x in cd ],marker='o',linestyle='none',color=c,label=texname(th['Xclusternames'][i])+' rel. diff.',markersize=4)
        if errorbars:
            if richardson_exponent_delta == None:
                cd_with_error = cd
            else:
                cd_with_error = [ x for x in cd if abs(x['fdrichardson']['observed_exponent'][i] - x['fdrichardson']['model_exponent'][i])<richardson_exponent_delta ]
            ax.plot([ x['R'] for x in cd_with_error ], [ x['errest'][errkey][i] for x in cd_with_error ],marker='_',linestyle='none',color=c,label=texname(cname)+' rel err estimate',markersize=8) #4
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
                    const = False,   # or False
                    dest='xlog')
    parser.add_argument('--ylog',
                    action='store_const',
                    default = None,
                    const = True,
                    dest='ylog')
    parser.add_argument('--ylinear',
                    action='store_const',
                    const = False,   # or False
                    dest='ylog')
    parser.add_argument('--absdiff',action='store_const',const='abs',dest='diffmode',help='Show absolute differences in right pane')
    parser.add_argument('--reldiff',action='store_const',const='rel',dest='diffmode',help='Show relative differences in right pane')
    parser.add_argument('--show-semiflat',default=[],nargs='+',dest='showsf',help='Which clusters to show semiflat approximation ("all", "none", or space-separated cluster names)')
    parser.add_argument('--pde-nmesh',default=None,help='Integer (to use single datum) or comma-separated list of integers (to use Richardson extrapolation)')
    parser.add_argument('--rmin',default=None,type=float)
    parser.add_argument('--rmax',default=None,type=float)
    parser.add_argument('--diffymin',default=None,type=float,help='y axis lower limit for difference plot (right pane)')
    parser.add_argument('--diffymax',default=None,type=float,help='y axis upper limit for difference plot (right pane)')
    parser.add_argument('--richardson-exponent-delta',default=0.4,type=float,help='Only show error bars when observed and expected exponents differ by less than this threshold')
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
        args.output = 'paperfig' + args.theoryname + '.pdf'

    if args.pde_nmesh != None:
        try:
            nm = int(args.pde_nmesh)
        except ValueError:
            # Assumption: Decreasing PDE_NMESH means increasing h
            nm = [ int(x) for x in args.pde_nmesh.split(',')]
            nm.sort(reverse=True)
        args.pde_nmesh = nm

    th = theorydata.getdata(args.theoryname)
    extra_kwargs = dict()
    for k in ['pde_nmesh','xlog','ylog','title_suffix','legend','legend_ncol','extra_left_yrange','tex_theoryname','richardson_exponent_delta','diffmode']:
        if getattr(args,k) != None:
            extra_kwargs[k]=getattr(args,k)
    f = spectral_coord_dual_plot(th,args.theta,rlims=(args.rmin,args.rmax),diffylims=(args.diffymin,args.diffymax),showsf=args.showsf,exclude=args.exclude,**extra_kwargs)
    f.savefig(args.output)
