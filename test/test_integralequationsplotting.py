import logging
import logconfig
logconfig.logconfig(filename=None)
logconfig.loglevel(logging.INFO)

import os
from integralequationsplotting import *

def test_integralequations_rayplots(tmpdir):
    R = 3
    steps = 1024
    tolerance = 1e-12
    theoryname = "A1A2"

    global xar
    xar = integralequations.computeXar("A1A2", R = R, tolerance = tolerance, steps = steps, method = "fourier", oper = True)

    plt = xar.rayplots(rayid = 0)[0]
    outpath = tmpdir.join('test_integralequations_rayplots.pdf')
    print('Writing {}'.format(outpath))
    plt.savefig(str(outpath))
    assert outpath.check(exists=1), "Expected output file is missing: {}".format(outpath)
