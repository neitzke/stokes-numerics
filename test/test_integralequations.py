import logging
import logconfig
logconfig.logconfig(filename=None)
logconfig.loglevel(logging.INFO)

from integralequations import *

def test_integralequations(theoryname = "A1A2", R = 0.4, oper = True):
    steps = 256
    tolerance = 1e-11

    def report(xar):
        """Print a few quantities derived from given xar"""
        print("cluster X-variables at theta = 0.00: %s" % xar.getCluster())
        print("cluster X-variables at theta = 0.10: %s" % xar.getCluster(theta = 0.10))
        print("cluster X-variables at theta = 0.20: %s" % xar.getCluster(theta = 0.20))

    print("Computing xars in theory %s with R = %0.8f, steps = %d, tolerance = %s" % (theoryname,R,steps,tolerance))

    global xarfourier,xarsfourier
    xarfourier = computeXar(theoryname = theoryname, R = R, oper = oper, tolerance = tolerance, steps = steps, method = "fourier")
    global xarsimps,xarssimps
    xarsimps = computeXar(theoryname = theoryname, R = R, oper = oper, tolerance = tolerance, steps = steps, method = "simps")

    print("----------------------------------------------------------------------")
    print("                          FOURIER METHOD")
    print("----------------------------------------------------------------------")
    report(xarfourier)

    print("----------------------------------------------------------------------")
    print("                          SIMPS METHOD")
    print("----------------------------------------------------------------------")
    report(xarsimps)

    print("----------------------------------------------------------------------")
    print("                          APPROXIMATE CLUSTER                         ")
    print("----------------------------------------------------------------------")
    for nterms in range(2):
        print("cluster X-variables at theta = 0.00, nterms = %d: %s" % (nterms,getApproxCluster(theoryname = theoryname, R = R, theta = 0, nterms = nterms, oper = oper)))
        if not oper:
            print("cluster X-variables at theta = 0.10, nterms = %d: %s" % (nterms,getApproxCluster(theoryname = theoryname, R = R, theta = 0.10, nterms = nterms, oper = oper)))
            print("cluster X-variables at theta = 0.20, nterms = %d: %s" % (nterms,getApproxCluster(theoryname = theoryname, R = R, theta = 0.20, nterms = nterms, oper = oper)))
