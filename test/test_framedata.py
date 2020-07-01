import logging
import logconfig
logconfig.logconfig(filename=None)
logconfig.loglevel(logging.INFO)

from framedata import *

def _moduletest(theoryname = "A1A2", R = 0.4, absh = 0.4, theta = 0.0, oper = True):

    def report(fd):
        print("cluster X-variables: %s" % list(zip(theorydata.getdata(theoryname)["Xclusternames"], fd.getCluster())))

    if not oper:
        nmesh = 255
        thresh = 1e-10

        print("Computing frames in theory %s with R = %0.8f, theta = %0.4f, pde_nmesh = %d, pde_thresh = %0.4e" % (theoryname, R, theta, nmesh, thresh))

        global fdfourier
        fdfourier = computeFrames(theoryname = theoryname, R = R, theta = theta, pde_nmesh = nmesh, pde_thresh = thresh, method = "fourier")

        global fdeuler
        fdeuler = computeFrames(theoryname = theoryname, R = R, theta = theta, pde_nmesh = nmesh, pde_thresh = thresh, method = "euler")

        print("----------------------------------------------------------------------")
        print("                          FOURIER METHOD")
        print("----------------------------------------------------------------------")
        report(fdfourier)

        print("----------------------------------------------------------------------")
        print("                          EULER METHOD")
        print("----------------------------------------------------------------------")
        report(fdeuler)
        print("----------------------------------------------------------------------")
        print("This would be saved to: {}".format(fd_filename(theoryname = theoryname, theta = theta, oper = oper, absh = absh, pde_nmesh = nmesh)))
        print("")

    if oper:
        print("Computing frames for oper in theory %s with abs(h) = %0.8f, theta = %0.4f" % (theoryname, absh, theta))
        fd = computeFrames(theoryname = theoryname, theta = theta, oper = True, absh = absh)
        print("----------------------------------------------------------------------")
        report(fd)
        print("----------------------------------------------------------------------")
        print("This would be saved to: {}".format(fd_filename(theoryname = theoryname, theta = theta, oper = oper, absh = absh)))
        print("")
        # Xrank = theorydata.getdata(theoryname)["Xrank"]
        # Z = theorydata.getdata(theoryname)["centralcharges"]
        # Xestimates = [exp(R*Z[n]) for n in range(Xrank)]
        # print("----------------------------------------------------------------------")
        # print("X estimates: %s" % Xestimates)

def test_framedata_oper():
    '''Compute frames for A1A2 oper example and print clusters'''
    _moduletest(oper=True)

def test_framedata_hitchin():
    '''Compute frames for A1A2 example by Fourier and Euler methods and print clusters'''
    _moduletest(oper=False)
