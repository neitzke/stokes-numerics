import logging
import logconfig
logconfig.logconfig(filename=None)
logconfig.loglevel(logging.INFO)

from comparisons import *

PDE_NMESH = 255

def one_cluster_test(theoryname, oper, delta_thresh, **kwargs):
    ccl = compareClusters(theoryname=theoryname, scratch=True, oper = oper, pde_nmesh=PDE_NMESH, **kwargs)
    print("--------------------")
    print(ccl)
    print("--------------------")
    for k,delta in enumerate(ccl['reldiff']):
        assert delta < delta_thresh, "Large difference in {tn} cluster {k}: rel_diff={delta} > {delta_thresh}".format(
            tn=theoryname,
            k=k,
            delta=delta,
            delta_thresh=delta_thresh
        )

def test_compare_clusters_A1A2_oper():
    one_cluster_test(theoryname="A1A2", absh=1, oper=True, delta_thresh = 1e-7)

def test_compare_clusters_A1A3_oper():
    one_cluster_test(theoryname="A1A3", absh=1, oper=True, delta_thresh = 1e-7, theta = 0.1, rmax = 5)

def test_compare_clusters_A2A1_oper():
    one_cluster_test(theoryname="A2A1", absh=1, oper=True, delta_thresh = 1e-7)

def test_compare_clusters_A2A2_oper():
    one_cluster_test(theoryname="A2A2", absh=1, oper=True, delta_thresh = 1e-7, tolerance = 2e-15)

def test_compare_clusters_A1A2_hitchin():
    one_cluster_test(theoryname="A1A2", R=1, oper=False, delta_thresh = 0.001)

def test_compare_clusters_A1A3_hitchin():
    one_cluster_test(theoryname="A1A3", R=1, oper=False, delta_thresh = 0.001)

def test_compare_clusters_A2A1_hitchin():
    one_cluster_test(theoryname="A2A1", R=1, oper=False, delta_thresh = 0.0015)

def test_compare_clusters_A2A2_hitchin():
    one_cluster_test(theoryname="A2A2", R=1, oper=False, delta_thresh = 0.04)
