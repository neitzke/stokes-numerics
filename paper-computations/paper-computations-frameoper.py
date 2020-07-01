'''Compute opers Stokes data for the examples in the Dumas-Neitzke experimental paper.'''

import itertools
import sys
import subprocess
import json
import os

# Presumed root of the repo = parent of script directory
pkgroot = os.path.realpath(os.path.join(sys.path[0], '..'))
# Allow imports from repo root
sys.path.insert(1, pkgroot)

import framedata as fd
import theorydata
import namegen

def computetheory(tn,thetalist,**kwargs):
    th = theorydata.getdata(tn)
    for theta in thetalist:
        for absh in th['abshlist']:
            yield 'import framedata as fd; extra_args={extra_args}; _,fn = fd.computeAndSaveFrames(\'{theoryname}\',oper=True,theta={theta},absh={absh},**extra_args); print(fn)'.format(
                extra_args=repr(kwargs),
                theoryname=tn,
                absh=absh,
                theta=theta
                )

def inventorytheory(tn,thetalist,**kwargs):
    '''Generator yielding a sequence of filenames that contain the calculations for a theory'''
    kwargs.pop('ode_thresh',None) # fd_filename does not accept this param
    kwargs.pop('rmax',None) # fd_filename does not accept this param
    th = theorydata.getdata(tn)
    for theta in thetalist:
        for absh in th['abshlist']:
            yield namegen.fd_filename(tn,oper=True,absh=absh,theta=theta,**kwargs)

if __name__ == '__main__':
    import argparse
    import logging

    parser = argparse.ArgumentParser(description=__doc__,epilog='Parameters not specified revert to the default values from framedata.py.')
    
    parser.add_argument('--ode-thresh',
                    type=float,
                    default=fd.ODE_THRESH)
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
    for k in ['ode_thresh']:
        if getattr(args,k) != None:
            extra_args[k] = getattr(args,k)
    logger.info('extra_args = {}'.format(extra_args))

    if args.mode == 'inventory':
        action = inventorytheory
        if args.fullpath:
            extra_args['fullpath'] = True
    else:
        action = computetheory

    gens = [
            action('A1A2',thetalist=[0.0,0.1],**extra_args),
            action('A1A2_Lambda=0.8j_cr=1_ci=0',thetalist=[0.0],**extra_args),
            action('A1A3',thetalist=[0.1],rmax=5,**extra_args), # poor accuracy with default rmax
            action('A2A1',thetalist=[0.0],**extra_args),
            action('A2A1_c=0.5',thetalist=[0.0],**extra_args),
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
