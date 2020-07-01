from __future__ import absolute_import
from __future__ import print_function
from scipy import sparse
import sys
import os.path
import fileinfo
import numpy as np

import logging
logger = logging.getLogger(__name__)

moduleversion = 2

def discdiff(g, kind = "laplacian"):
    '''Make a matrix representing a discrete differential operator for a square mesh

    Args:
      g = SquareGrid object representing the mesh
    '''

#   sys.stderr.write('PDE: Generating discrete laplacian matrix\n')
    if g.dx != g.dy or g.nx != g.ny:
        raise NotImplementedError("Discrete differential operators for non-uniform grid spacing not supported")

    n = g.nx
    N = g.N

    filename = "grid%d.v%d.disc%s.cache" % (n,moduleversion,kind)
    pathname = os.path.join(fileinfo.CACHEPATH, filename)

    if os.path.isfile(pathname + ".npz"):
        logger.info("Loading discrete %s from %s.npz" % (kind,filename))
        m = sparse.load_npz(pathname + ".npz")
    else:
        logger.info("Generating discrete %s with size %d" % (kind,n))
        m = sparse.dok_matrix( (N,N) )

        if kind == "laplacian":
            # Interior points -- centered finite difference 2nd deriv
            for i in range(1,n-1):
                for j in range(1,n-1):
                    # d/dx^2
                    m[g.toidx(i,j),g.toidx(i-1,j)] += 1
                    m[g.toidx(i,j),g.toidx(i+1,j)] += 1
                    m[g.toidx(i,j),g.toidx(i,j)] += -2
                    # d/dy^2
                    m[g.toidx(i,j),g.toidx(i,j+1)] += 1
                    m[g.toidx(i,j),g.toidx(i,j-1)] += 1
                    m[g.toidx(i,j),g.toidx(i,j)] += -2
        elif kind == "ddx":
            # non-centered version
            # TODO: see if symmetrizing makes a difference
            for i in range(1,n-1):
                for j in range(1,n-1):
                    m[g.toidx(i,j),g.toidx(i,j)] += -1
                    m[g.toidx(i,j),g.toidx(i+1,j)] += 1
        elif kind == "ddy":
            # non-centered version
            # TODO: see if symmetrizing makes a difference
            for i in range(1,n-1):
                for j in range(1,n-1):
                    m[g.toidx(i,j),g.toidx(i,j)] += -1
                    m[g.toidx(i,j),g.toidx(i,j+1)] += 1

        m = m.tocsr()

        logger.info("Saving discrete %s to %s.npz" % (kind,filename))
        sparse.save_npz(pathname, m)

    idx = 1.0 / (g.dx*g.dx)
    mscaled = m*idx

    return mscaled
