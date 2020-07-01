'''Compute IEQ values for Hitchin Stokes data of examples in the Dumas-Neitzke experimental paper.'''
import sys
import itertools
import subprocess
import json
import os

# Presumed root of the repo = parent of script directory
pkgroot = os.path.realpath(os.path.join(sys.path[0], '..'))
# Allow imports from repo root
sys.path.insert(1, pkgroot)

import namegen
import integralequations as ieq
import theorydata

def computetheory(tn,thetalist,**kwargs):
    th = theorydata.getdata(tn)
    for theta in thetalist:
        for R in th['Rlist']:
            yield 'import integralequations as ieq; _,fn = ieq.computeAndSaveXar(\'{theoryname}\',oper=False,savexar=False,saveclusters=True,thetalist=[{theta}],R={R}); print(fn)'.format(
            theoryname=tn,
            R=R,
            theta=theta
            )

def inventorytheory(tn,thetalist,**kwargs):
    '''Generator yielding a sequence of filenames that contain the calculations for a theory'''
    th = theorydata.getdata(tn)
    for theta in thetalist:
        for R in th['Rlist']:
            yield namegen.ie_filename(tn,oper=False,R=R,theta=theta,clusteronly=True,**kwargs)

if __name__ == '__main__':
    import argparse
    import logging

    parser = argparse.ArgumentParser(description=__doc__)
    
    parser.add_argument('-v','--verbose',action='store_true')

    parser.add_argument('--run',dest='mode',action='store_const',const='run',help='Run the calculations (one at a time)')
    parser.add_argument('--print',dest='mode',action='store_const',const='print',help='Print shell commands to run the calculations (for easy parallelization)')
    parser.add_argument('--json',dest='mode',action='store_const',const='json',help='Print one json object for each calculation, with a "script" attribute containing python code to perform the work')
    parser.add_argument('--inventory',dest='mode',action='store_const',const='inventory',help='Print a list of files that would be created by these calculations')
    parser.add_argument('--fullpath',action='store_true',help='In inventory mode, print the full path')

    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)

    if args.mode == None:
        logging.info('Using default mode "print" since no mode specified.')
        args.mode = 'print'

    extra_args = dict()
    if args.mode == 'inventory':
        action = inventorytheory
        if args.fullpath:
            extra_args['fullpath'] = True
    else:
        action = computetheory
    
    gens = [
            action('A1A2',thetalist=[0.0,0.1],**extra_args),
            action('A1A2_Lambda=0.8j_cr=1_ci=0',thetalist=[0.0],**extra_args),
            action('A1A3',thetalist=[0.1],**extra_args),
            action('A2A1',thetalist=[0.0],**extra_args),
            action('A2A2',thetalist=[0.1],**extra_args),
           ]

    for script in itertools.chain(*gens):
        if args.mode == 'run':
            # Run in subinterpreter to make run and print options truly equivalent
            subprocess.call([sys.executable,'-c',script],cwd=pkgroot)
        elif args.mode == 'print':
            print(sys.executable + ' -c "' + script + '"')
        elif args.mode == 'json':
            d = dict(script=script)
            print(json.dumps(d))
        elif args.mode == 'inventory':
            print(script)
