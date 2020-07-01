from __future__ import print_function
import sys
import subprocess
import os

sys.path.insert(1, os.path.join(sys.path[0], 'paper-figures'))
paperfig_oper = __import__("paperfig-oper")

def test_paperfig_oper(tmpdir):
    scriptfn = paperfig_oper.__file__
    interpfn = sys.executable

    # Modify environment so the script doesn't find any frame or xar data
    # to make the test independent of the filesystem outside of the repo
    env = os.environ.copy()
    env["XARPATH"] = str(tmpdir.join('xars'))
    env["FRAMEPATH"] = str(tmpdir.join('frames'))
    
    outpath=tmpdir.join('paperfig-oper-test.pdf')
    subprocess.call([interpfn,scriptfn,'A1A2','0.0','-o',str(outpath)],env=env)

    assert outpath.check(exists=1), "Expected output file is missing: {}".format(outpath)
