"""Make a LaTeX table of hitchin section comparisons for one theory"""
import numpy as np
import math
import sys
import os
import argparse
import logging

logger = logging.getLogger(__name__)

# Allow imports from parent directory
sys.path.insert(1, os.path.join(sys.path[0], '..'))

import theorydata
import comparisons
import richardson
from latexutil import *

parser = argparse.ArgumentParser(description=__doc__)

parser.add_argument('--scratch',action='store_true',help='Compute from scratch; do not load from files')
parser.add_argument('--richardson-exponent-delta',default=0.4,type=float,help='Only show pde error estimates when observed and expected exponents differ by less than this threshold')
parser.add_argument('--pde-nmesh',help='List of PDE_NMESH values to apply Richardson interpolation to',default='2047,4095,8191')
parser.add_argument('-v','--verbose',action='store_true')

args = parser.parse_args()

if args.verbose:
    logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

pde_nmesh = [ int(x) for x in args.pde_nmesh.split(',') ]
pde_nmesh.sort(reverse=True)
if len(pde_nmesh) != 3:
    logger.error('Only Richardson interpolation from 3 points is supported now')
    sys.exit(1)

tn = 'A1A2'
th = theorydata.getdata(tn)

table_log_rs = [ -6.0, -3.0, 0.0, 1.5, 1.875 ]

def approx_in(x,L):
    return min(abs(x-y) for y in L) < 0.001

Rlist  = [ r for r in th['Rlist'] if approx_in(math.log(r),table_log_rs) ]

if args.scratch:
    verb = 'Computing'
else:
    verb = 'Loading'

L = []
for i,R in enumerate(Rlist):
    sys.stderr.write('{} for R={} ({} of {})\n'.format(verb,R,i+1,len(Rlist)))
    d = comparisons.compareClusters(th['name'],R=R,theta=0.1,oper=False,scratch=args.scratch,pde_nmesh=pde_nmesh)
    d['R'] = R
    d['logR'] = math.log(R)
    L.append(d)

print(r'''\begin{tabular}{@{}lllll@{}}
\toprule
& & & & Rel.~PDE\\
$R$ & $X_1^{\DE}$ & $X_1^{\IEQ}$ & $\mathrm{reldiff}(X_1^{\DE},X_1^{\IEQ})$ & error est. \\
\midrule''')

for r in sorted(L,key=lambda x:x['R']):
    x1de, x2de = [-x for x in r['fdcluster']]    # convert from code to paper sign convention
    x1ieq, x2ieq = [-x for x in r['xarcluster']] # convert from code to paper sign convention
    x1reldiff, x2reldiff = r['reldiff']
    x1relerr, x2relerr = r['errest']['relpde']
    if abs(r['fdrichardson']['observed_exponent'][0] - r['fdrichardson']['model_exponent'][0]) < args.richardson_exponent_delta:
        errstr = latex_real(x1relerr,prec=1)
    else:
        errstr = '--'
    
    row = [ latex_exp(r['logR']), 
            latex_real(x1de,prec=9), 
            latex_real(x1ieq,prec=9), 
            latex_real(x1reldiff,prec=1),
            errstr ]
    print(' & '.join(row) + r'\\')
print(r'''\midrule
& & & & Rel.~PDE \\
$R$ & $X_2^{\DE}$ & $X_2^{\IEQ}$ & $\mathrm{reldiff}(X_2^{\DE},X_2^{\IEQ})$ & error est.\\
\midrule''')
for r in sorted(L,key=lambda x:x['R']):
    x1de, x2de = [-x for x in r['fdcluster']]    # convert from code to paper sign convention
    x1ieq, x2ieq = [-x for x in r['xarcluster']] # convert from code to paper sign convention
    x1reldiff, x2reldiff = r['reldiff']
    x1relerr, x2relerr = r['errest']['relpde']
    if abs(r['fdrichardson']['observed_exponent'][1] - r['fdrichardson']['model_exponent'][1]) < args.richardson_exponent_delta:
        errstr = latex_real(x2relerr,prec=1)
    else:
        errstr = '--'

    row = [ latex_exp(r['logR']), 
            latex_real(x2de,prec=9), 
            latex_real(x2ieq,prec=9), 
            latex_real(x2reldiff,prec=1),
            errstr ]
    print(' & '.join(row) + r'\\')
print(r'\bottomrule')
print(r'\end{tabular}')
