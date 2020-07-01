from squaregrid import *

def test_squaregrid():
    gr = SquareGrid(3.0,6)
    expected_indices = [ (0,0), (1,0), (2,0), (3,0), (4,0), (5,0),
                         (0,1), (1,1), (2,1), (3,1), (4,1), (5,1),
                         (0,2), (1,2), (2,2), (3,2), (4,2), (5,2),
                         (0,3), (1,3), (2,3), (3,3), (4,3), (5,3),
                         (0,4), (1,4), (2,4), (3,4), (4,4), (5,4),
                         (0,5), (1,5), (2,5), (3,5), (4,5), (5,5), ]

    for k,er in enumerate(expected_indices):
        assert gr.fromidx(k) == er

    # Index 7 is the second element of row 2, whicfh should be -1.8-1.8j
    # because the grid spacing is 1.2 = (2*3.0)/(6-1)
    assert gr.zv[7]==-1.8-1.8j

