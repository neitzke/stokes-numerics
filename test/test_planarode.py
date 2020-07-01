import logging
import logconfig
logconfig.logconfig(filename=None)
logconfig.loglevel(logging.INFO)

from planarode import *

def test_planarode():
    '''Integrate oper connection along 200 rays and print results'''
    global ODE
    ODE = planarODE.buildoperODEfrompolys(Pcoefs = [[0.5,-1,0,1]], K = 2, r = 5.0)
    vlist = ODE.integrate_rays(200)
    for v in enumerate(vlist):
        print(v)
