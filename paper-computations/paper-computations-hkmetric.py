'''Compute hyperkahler metric data for the Dumas-Neitzke experimental paper.'''

import sys
import subprocess
import json
import os

from math import exp

# Presumed root of the repo = parent of script directory
pkgroot = os.path.realpath(os.path.join(sys.path[0], '..'))
# Allow imports from repo root
sys.path.insert(1, pkgroot)

def computemetricdata(**kwargs):
    clist = [0.03*exp(0.2*n) for n in range(25)]
    for c in clist:
        yield 'import hkmetric; extra_args={extra_args}; hkmetric.computeAndSaveGs(c={c},**extra_args)'.format(
            extra_args=repr(kwargs),
            c=c
            )

if __name__ == '__main__':
    import argparse
    import logging

    parser = argparse.ArgumentParser(description=__doc__,epilog='Parameters not specified revert to the default values from framedata.py.')
    
    parser.add_argument('--pde-nmesh', type=int)
    parser.add_argument('--eps', type=float)
    parser.add_argument('--Lambda', default=0.0, type=float)
    parser.add_argument('-v','--verbose',action='store_true')

    parser.add_argument('--run',dest='mode',action='store_const',const='run',help='Run the calculations (one at a time)')
    parser.add_argument('--print',dest='mode',action='store_const',const='print',help='Print shell commands to run the calculations (for easy parallelization)')
    parser.add_argument('--json',dest='mode',action='store_const',const='json',help='Print one json object for each calculation, with a "script" attribute containing python code to perform the work')

    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)

    if args.mode == None:
        logging.info('Using default mode "print" since no mode specified.')
        args.mode = 'print'
    
    extra_args = dict()
    for k in ['pde_nmesh','eps','Lambda']:
        if getattr(args,k) != None:
            extra_args[k] = getattr(args,k)
    logger.info('extra_args = {}'.format(extra_args))

    for script in computemetricdata(**extra_args):
        if args.mode == 'run':
            # Run in subinterpreter to make run and print options truly equivalent
            subprocess.call([sys.executable,'-c',script],cwd=pkgroot)
        elif args.mode == 'print':
            print(sys.executable + ' -c "' + script + '"')
        elif args.mode == 'json':
            d = dict(script=script)
            print(json.dumps(d))
