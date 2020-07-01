# Scripts to create figures and tables

This directory holds scripts that read data files created by the scripts in `../papercomputations/` and reproduce the tables and plots in the paper _Opers and nonabelian Hodge: numerical studies_.
These scripts require matplotlib.

The `Makefile` contains rules to generate the individual tables and plots by calling the scripts `papertable-*.py` and `paperfig-*.py`.  Thus, to generate all of them, one need only run

```
make
```

in this directory.  This takes about 30 minutes on a fast CPU (as of mid-2020), and can be sped up significantly with a parallel make, e.g. `make -j8`.