"""Make a LaTeX table of oper comparisons for one theory"""
import math
import sys
import os
import argparse
import logging

# Allow imports from parent directory
sys.path.insert(1, os.path.join(sys.path[0], '..'))

import theorydata
import comparisons
from latexutil import *

parser = argparse.ArgumentParser(description=__doc__)

parser.add_argument('--scratch',action='store_true',help='Compute from scratch; do not load from files')
parser.add_argument('-v','--verbose',action='store_true')

args = parser.parse_args()

if args.verbose:
    logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

tn = 'A1A2_Lambda=0.8j_cr=1_ci=0'
th = theorydata.getdata(tn)

table_log_abshinvs = [ -6.0, -3.0, 0.0, 1.5, 2.25, 3.0]

def approx_in(x,L):
    return min(abs(x-y) for y in L) < 0.001

abshlist = [ h for h in th['abshlist'] if approx_in(-math.log(h),table_log_abshinvs) ]

if args.scratch:
    verb = 'Computing'
else:
    verb = 'Loading'

L = []
for i,h in enumerate(abshlist):
    sys.stderr.write('{} for h={} ({} of {})\n'.format(verb,h,i+1,len(abshlist)))
    d=comparisons.compareClusters(th['name'],absh=h,theta=0.0,oper=True,scratch=args.scratch)
    d['abshinv'] = 1.0/h
    d['logabshinv'] = -math.log(h)
    L.append(d)

print(r'''\begin{tabular}{@{}lllll@{}}
\toprule
& & & & Rel.~ODE\\
$|\hbar|^{-1}$ & $X_1^{\DE}$ & $X_1^{\IEQ}$ & $\mathrm{reldiff}(X_1^{\DE},X_1^{\IEQ})$ & err.~est. \\
\midrule''')
for r in sorted(L,key=lambda x:x['logabshinv']):
    x1de, x2de = [-x for x in r['fdcluster']]    # convert from code to paper sign convention
    x1ieq, x2ieq = [-x for x in r['xarcluster']] # convert from code to paper sign convention
    x1reldiff, x2reldiff = r['reldiff']
    x1relerr, x2relerr = r['errest']['relode']
    row = [ latex_exp(r['logabshinv']), 
            latex_real(x1de.real,prec=12), 
            latex_real(x1ieq.real,prec=12), 
            latex_real(x1reldiff,prec=1),
            latex_real(x1relerr,prec=1) ]
    print(' & '.join(row) + r'\\')
    row = [ "", 
            latex_real(x1de.imag,prec=12)+"i", 
            latex_real(x1ieq.imag,prec=12)+"i", 
            "",
            "" ]
    print(' & '.join(row) + r'\\')
print(r'''\bottomrule
\end{tabular}''')
print(r'''\begin{tabular}{@{}lllll@{}}
\toprule
& & & & Rel.~ODE\\
$|\hbar|^{-1}$ & $X_2^{\DE}$ & $X_2^{\IEQ}$ & $\mathrm{reldiff}(X_2^{\DE},X_2^{\IEQ})$ & err.~est. \\
\midrule''')
for r in sorted(L,key=lambda x:x['logabshinv']):
    x1de, x2de = [-x for x in r['fdcluster']]    # convert from code to paper sign convention
    x1ieq, x2ieq = [-x for x in r['xarcluster']] # convert from code to paper sign convention
    x1reldiff, x2reldiff = r['reldiff']
    x1relerr, x2relerr = r['errest']['relode']
    row = [ latex_exp(r['logabshinv']), 
            latex_real(x2de.real,prec=12), 
            latex_real(x2ieq.real,prec=12), 
            latex_real(x2reldiff,prec=1),
            latex_real(x2relerr,prec=1) ]
    print(' & '.join(row) + r'\\')
    row = [ "", 
            latex_real(x2de.imag,prec=12)+"i", 
            latex_real(x2ieq.imag,prec=12)+"i", 
            "",
            "" ]
    print(' & '.join(row) + r'\\')
print(r'\bottomrule')
print(r'\end{tabular}')
