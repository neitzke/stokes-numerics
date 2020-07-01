import logging
import logconfig
logconfig.logconfig(filename=None)
logconfig.loglevel(logging.INFO)

import squaregrid
from sde import *

def test_sde():
    def c(z):
        return 2.0*z*z*z - 1j*z + 0.2

    gr = squaregrid.SquareGrid(3.0,255)

    def report(Q):
        """Print data about solution in SelfDualityEquation object Q"""
        j = int(gr.ny / 2)
        for i in range(0,gr.nx,gr.nx // 10):
            z = Q.grid.zm[j,i]
            u = Q.u[j,i]
            u0 = Q.u0[j,i]
            print('u(%g%+gi) = \t%f  (diff from uzero is %f)' % (z.real,z.imag,u,u-u0))
            
    print("----------------------------------------------------------------------")
    print("                          FOURIER METHOD")
    print("----------------------------------------------------------------------")
    global QF
    QF = SelfDualityEquation(3,c,gr,method='fourier')
    report(QF)

    print("----------------------------------------------------------------------")
    print("                          EULER METHOD")
    print("----------------------------------------------------------------------")
    global QE
    QE = SelfDualityEquation(3,c,gr,method='euler')
    report(QE)
