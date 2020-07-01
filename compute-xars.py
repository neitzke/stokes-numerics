'''Use integral equations to compute clusters for given theory and parameters'''
from __future__ import absolute_import
from __future__ import print_function
import sys
import argparse

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('theory',help='Theory name')
parser.add_argument('R',type=float)
parser.add_argument('-v','--verbose',action='store_true')
parser.add_argument('-m','--method',default='fourier',help='Computation method (fourier or simps)')
parser.add_argument('-t','--theta',default=0.0,type=float,help='Theta (only used if clusters are requested)')
parser.add_argument('--L',type=float)
parser.add_argument('--tolerance',type=float)
parser.add_argument('--steps',type=int)
parser.add_argument('-o','--output',help='Output filename for xar (default: do not save)')
parser.add_argument('-c','--cluster-output',help='Output filename for clusters (default: do not save)')
parser.add_argument('--outdir',help='Choose a filename automatically and save in this directory')
parser.add_argument('-j','--print-json',action='store_true',help='Print clusters to stdout as JSON')
parser.add_argument('--no-clobber',action='store_true',help='Exit if output file already exists')

args = parser.parse_args()

if args.verbose:
    import logging
    from logconfig import logconfig
    logconfig(filename=None)

import theorydata
import integralequations as ieq
import namegen
import os

if not args.theory in theorydata.theorydata:
    sys.stderr.write('Unknown theory: %s\n' % args.theory)
    sys.exit(1)

if args.outdir:
    args.output = os.path.join(args.outdir,namegen.ie_filename(args.theory,R=args.R,oper=False,fullpath=False))
    args.cluster_output = os.path.join(args.outdir,namegen.ie_filename(args.theory,R=args.R,oper=False,fullpath=False,clusteronly=True))

if args.output and args.no_clobber:
    if os.path.exists(args.output):
        sys.stderr.write('ERROR: Output file \"{}\" exists.\n'.format(args.output))
        sys.exit(1)
if args.cluster_output and args.no_clobber:
    if os.path.exists(args.cluster_output):
        sys.stderr.write('ERROR: Output file \"{}\" exists.\n'.format(args.cluster_output))
        sys.exit(1)

# handle optional parameters that we want to fall back to IEQ defaults for
extra_args = dict()
for k in [ 'method', 'L', 'tolerance', 'steps' ]:
    v = getattr(args,k)
    if v != None:
        extra_args[k] = v

th = theorydata.getdata(args.theory)
x = ieq.computeXar(args.theory, R=args.R, **extra_args)


outobj = {
    'clusters': x.getCluster(theta=args.theta),
}

if args.print_json:
    import json
    print(json.dumps(outobj))
else:
    for k in outobj:
        print('{}: {}'.format(k,outobj[k]))

if args.output:
    x.save(args.output)

if args.cluster_output:
    x.saveclusters(args.cluster_output, args.theta)
