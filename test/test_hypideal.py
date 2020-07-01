import sys
from hypideal import *
    
def test_h2line(tmpdir,t0=-3.2,t1=1.5):
    L = H2Line(t0,t1)
    for k in dir(L):
        if k.startswith('__'):
            continue
        print('{}: {}'.format(k,getattr(L,k)))
    
    import matplotlib.pyplot as plt
    from matplotlib import rc

    rc('figure',figsize=(9,9))
    fig = plt.figure()
    ax = plt.gca()
    ax.axis('off')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim(-1.05,1.05)
    ax.set_ylim(-1.05,1.05)
    ax.set_aspect(1)
    ax.add_patch(patches.Circle((0,0),radius=1,fc='none',color='gray',lw=3))
    ax.add_patch(L.patch(fc='none',color='red',lw=2))
    plt.tight_layout()
    outpath = tmpdir.join('test_h2line.pdf')
    print('Writing {}'.format(outpath))
    plt.savefig(str(outpath))
    assert outpath.check(exists=1), "Expected output file is missing: {}".format(outpath)

def test_h2polygon(tmpdir,verts = [0.1, 0.5, 3.9]):
    P = H2IdealPolygon(verts)

    import matplotlib.pyplot as plt
    from matplotlib import rc

    rc('figure',figsize=(9,9))
    fig = plt.figure()
    ax = plt.gca()
    ax.axis('off')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim(-1.05,1.05)
    ax.set_ylim(-1.05,1.05)
    ax.set_aspect(1)
    ax.add_patch(patches.Circle((0,0),radius=1,fc='none',color='gray',lw=3,zorder=1))
    ax.add_patch(P.patch(fc='red',lw=0,zorder=5))
    pc = PatchCollection(P.vertex_patches(radius=0.02),facecolor='orange',zorder=100)
    ax.add_collection(pc)
    plt.tight_layout()
    
    outpath=tmpdir.join('test_h2polygon.pdf')
    print('Writing {}'.format(outpath))
    plt.savefig(str(outpath))
    assert outpath.check(exists=1), "Expected output file is missing: {}".format(outpath)
