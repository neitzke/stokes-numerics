import logging
import logconfig
logconfig.logconfig(filename=None)
logconfig.loglevel(logging.INFO)

import squaregrid
from discdiff import *
import pytest

@pytest.mark.filterwarnings("ignore:the matrix subclass")
def test_discdiff():
    '''Create and print finite difference derivative matrices for a small square grid'''
    # These tests just make sure the code runs without an exception, no real check!
    gr = squaregrid.SquareGrid(3.0,5)
    print(gr)
    print(discdiff(gr,"ddx").todense())
    print(discdiff(gr,"ddy").todense())
    print(discdiff(gr,"laplacian").todense())
