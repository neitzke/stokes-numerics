import logging
import logconfig
logconfig.logconfig(filename=None)
logconfig.loglevel(logging.INFO)

import hkmetric
import numpy as np

def test_compare_metrics():
	data = hkmetric.comparemetrics(pde_nmesh = 300)
	assert abs(data["fd"]-data["ieq"]) <  0.001, "difference between PDE and IEQ computations unexpectedly big"
