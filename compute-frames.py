'''Compute frame data for given theory and parameters'''
# Note: This standalone executable will save framedata files outside the 
# usual path configured in paths.conf
from __future__ import absolute_import
from __future__ import print_function
import sys
import argparse
import theorydata
import framedata
import namegen
import os

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('theory',help='Theory name')
parser.add_argument('R',type=float)
parser.add_argument('-v','--verbose',action='store_true')
parser.add_argument('-m','--method',default='fourier',help='Computation method (fourier or euler)')
parser.add_argument('-t','--theta',default=0.0,type=float)
parser.add_argument('--pde-nmesh',type=int,default=framedata.PDE_NMESH)
parser.add_argument('--pde-thresh',type=float)
parser.add_argument('--pde-maxiter',type=int)
parser.add_argument('--ode-thresh',type=float)
parser.add_argument('--ode-rstep',type=float)
parser.add_argument('--ode-nsteps',type=int)
parser.add_argument('-s','--store-metric',action='store_true',help='Save metric to output file')
parser.add_argument('-o','--output',help='Output filename (default: do not save)')
parser.add_argument('--outdir',help='Choose a filename automatically and save in this directory')
parser.add_argument('-j','--print-json',action='store_true',help='Print vertices and clusters as a JSON object (default: print as a table)')
parser.add_argument('--no-clobber',action='store_true',help='Exit if output file already exists')

args = parser.parse_args()

if args.verbose:
    import logging
    from logconfig import logconfig
    logconfig(filename=None)


if not args.theory in theorydata.theorydata:
    sys.stderr.write('Unknown theory: %s\n' % args.theory)
    sys.exit(1)

if args.outdir:
    args.output = os.path.join(args.outdir,namegen.fd_filename(args.theory,R=args.R,theta=args.theta,oper=False,fullpath=False,pde_nmesh=args.pde_nmesh))

if args.output and args.no_clobber:
    if os.path.exists(args.output):
        sys.stderr.write('ERROR: Output file \"{}\" exists.\n'.format(args.output))
        sys.exit(1)

if args.method == 'fourier' and args.pde_nmesh == None:
    args.pde_nmesh = 511

# handle optional parameters that we want to fall back to framedata defaults for
extra_args = dict()
for k in ['pde_nmesh', 'pde_thresh', 'pde_maxiter', 'ode_thresh', 'ode_rstep', 'ode_nsteps']:
    v = getattr(args,k)
    if v != None:
        extra_args[k] = v

th = theorydata.getdata(args.theory)
fd = framedata.framedata(theory=th,
                         R=args.R,
                         theta=args.theta,
                         storemetric=args.store_metric,
                         method=args.method,
                         **extra_args)

outobj = {
    'clusters': fd.getCluster(),
    'vertices': fd.vertexlist
}

if args.print_json:
    import json
    # Convert numpy.complex128 to tuples
    outobj['vertices'] = [ [ (x.real, x.imag) for x in vec ] for vec in outobj['vertices'] ]
    print(json.dumps(outobj))
else:
    for k in outobj:
        print('{}: {}'.format(k,outobj[k]))


if args.output:
    fd.save(args.output)
