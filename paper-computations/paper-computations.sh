#!/bin/bash

# Reproduce all of the computations necessary to generate the tables and
# figures in the associated paper "Opers and nonabelian Hodge: numerical studies"

# These commands will generate the data files and save them to paths
# configured in ../paths.conf (or subdirectories of /tmp if paths.conf does
# not exist).  The actual figures and tables are generated by the scripts in
# ../paper-figures

usage() {
    echo "Usage: $0 MODE"
    echo "where MODE is one of:"
    echo "  --run   : to run all computations serially"
    echo "  --print : to print a few thousand one-line shell scripts that can be run in parallel"
    exit 1
}

if [ "$#" -ne 1 ]; then
    usage
fi

MODE=$1

[ "$MODE" == "--run" ] || [ "$MODE" == "--print" ] || ( echo "Unknown mode \"$MODE\""; usage )

# We assume the paper-computations-*.py scripts are in the same directory
# as this script, but that may not be CWD.  Here we determine the script
# directory:
SCRIPTDIR=$(readlink -f "$(dirname "${BASH_SOURCE[0]}")")

/usr/bin/env python3 $SCRIPTDIR/paper-computations-xarhitchin.py $MODE
/usr/bin/env python3 $SCRIPTDIR/paper-computations-xaroper.py $MODE
/usr/bin/env python3 $SCRIPTDIR/paper-computations-frameoper.py $MODE
/usr/bin/env python3 $SCRIPTDIR/paper-computations-framehitchin.py --pde-nmesh 2047 $MODE
/usr/bin/env python3 $SCRIPTDIR/paper-computations-framehitchin.py --pde-nmesh 4095 $MODE
/usr/bin/env python3 $SCRIPTDIR/paper-computations-framehitchin.py --pde-nmesh 8191 $MODE
/usr/bin/env python3 $SCRIPTDIR/paper-computations-hkmetric.py --pde-nmesh 1400 $MODE