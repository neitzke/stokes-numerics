import logging
import logconfig
logconfig.logconfig(filename=None)
logconfig.loglevel(logging.INFO)

from nlp_euler import *
import squaregrid

def test_nlp_euler():
    def u0(z):
        return 5.0
    
    def f(u,z):
        return np.exp(u) - np.abs(z**3 - 3.0*z + 1) * np.exp(-u)

    def d1f(u,z):
        return np.exp(u) + np.abs(z**3 - 3.0*z + 1) * np.exp(-u) 

    g = squaregrid.SquareGrid(4.0,127);
    S = NLPEuler(f,d1f,np.vectorize(u0),g)
    print('Solution u has shape: ',S.u.shape)
