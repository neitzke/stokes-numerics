import logging
import logconfig
logconfig.logconfig(filename=None)
logconfig.loglevel(logging.INFO)

from extsde import *

def test_extsde(K = 2, method = "euler"):
    '''Solve the self duality equation, integrate along many rays'''
    coefs = [0.5, -1, 0, 1]

    H = ExtendedSelfDualityEquation(K=K, coefs=coefs, pde_nmesh = 100, rmax = 5.0, pde_thresh = 1e-7, pde_maxiter = 5000, method = method)
    ODE = planarODE.buildhiggsODE(H)
    vlist = ODE.integrate_rays(200)
    for v in enumerate(vlist):
        print(v)
