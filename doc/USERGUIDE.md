# stokes-numerics:  User Guide

## Basic use

Much of the mathematical background is desrcibed in the paper _Opers and nonabelian Hodge: numerical studies_ by David Dumas and Andrew Neitzke (henceforth called "the paper").  This guide focuses on details of using the code itself.

The programs and modules are meant to be used from the root of a local clone of the repository.  Several parts of the code embed the hash of the current commit in their output files, so it is necessary to run them in a git repository.

The main functions are:

1. Compute Stokes data and flat connections associated to polynomial opers and polynomial Higgs bundles
2. Solve integral equations to find predicted Stokes data according to the Twistorial Riemann-Hilbert (TRH) conjecture of Gaiotto-Moore-Neitzke and Gaiotto.
3. Compute Hitchin's hyperkahler metric tensor restricted to the Hitchin section
4. Solve integral equations to find predicted values for the metric tensor from (3) according to TRH.

Of these, (1) supports all polynomial examples in ranks 2 and 3, while (2)-(4) work only in certain hard-coded families.  For (2) the families are enumerated in `theorydata.py`; for (3) and (4) only a single family (z^3 - Lambda z - c)dz^2 in rank 2 is supported.

## Configuration

This package was designed for experiments to compare TRH conjecture predictions to the results of direct calculations based on solving the relevant PDE or ODE.  Many of the calculations are long-running, and so the design assumes most work will be split into two steps:

1. Computing Stokes data (or metric data) through long-running integral equation, PDE, or ODE calculations, saving the results to JSON files.
2. Interactive or script-based work to collect and analyze the results after loading from files.

Data files are automatically named based on calculation parameters, and are stored in directories specified in an INI-style configuration file `paths.conf`.  An example `paths.conf-example` is included in the repository.  If `paths.conf` is missing, subdirectories of `/tmp` will be used.

A typical `paths.conf` looks like

```
[paths]
XARPATH = /home/username/data/xars
FRAMEPATH = /home/username/data/frames
CACHEPATH = /home/username/data/cache
HKMETRICPATH = /home/username/hk
```

with the variable meanings as follows:

* `XARPATH` - Results of integral equation calculations (called "xars" or "xarrays" in various parts of the code)
* `FRAMEPATH` - Results of PDE/ODE calculations (which consists of horizontal frames and parallel transport data for certain flat connections)
* `CACHEPATH` - Cache area for discrete differentiation matrices
* `HKMETRICPATH` - Hyperkaehler metric data (both direct and TRH conjectural)

Modules that deal with data files include logic to select filenames automatically based on the calculation parameters.

Any of the paths in `paths.conf` can be overridden for a single calculation by setting the corresponding environment variable.

## Sample calculations

A jupyter/ipython notebook of sample calculations with explanation is included in [sample-computations.ipynb](../sample-computations.ipynb).

## Logical structure

The core computational functionality is concentrated in the following modules:

* `theorydata.py` contains a list of "theories" (parameterized families of Higgs bundles and opers) for which both integral equation and direct PDE/ODE computations are possible. Its main function `getdata(theoryname)` takes a string `theoryname` (e.g. `"A1A2"`, `"A1A3"`, `"A2A1"`, `"A2A2"`, or `"A1A2_Lambda=0.8j_cr=1_ci=0"`) and returns a dict that contains homological data about the spectral curve, periods, and explicit forms of the spectral coordinates in terms of Stokes matrices.  In most cases the theory data is computed by hand and simply tabulated here; a few cases have free parameters and compute periods numerically on request.  Many functions in other modules (e.g. most in `comparisons`, and the constructors of the main classes of `framedata` and `integralequations`) expect either a theory name or a `theorydata`

* `comparisons` loads (or computes) Stokes data for one of the named families by both methods (direct and integral equations) and compares them.  Plots in the paper are produced by scripts that use this module to collect comparison data.  The actual calculations (or loading of data) is handled by the modules `framedata` and `integralequations`.

* `framedata` is the main interface to PDE/ODE Stokes data calculations and saving/loading such data.  Its main class `framedata.framedata` is the abstraction of a polynomial oper or Higgs bundle paired with computation of its Stokes data, parallel transports, etc..  This class provides a higher-level interface to the PDE solver in module `extsde` (solving the self-duality equation) and the ODE solver in module `planarode` (computing parallel transports).  It also contains the logic for turning connection data into spectral coordinates by computing certain rational functions of the parallel transport matrices, in the method `getCluster()`.

* `integralequations` handles all integral equation calculations, implementing the iteration procedure described in the paper to compute predicted values of spectral coordinates.  Its main class `integralequations.xarray` represents a polynomial oper or Higgs bundle paired with arrays of values of its spectral coordinates on certain rays in the complex plane.  Some of the logic for setting up the iteration is handled by the `theory` module, whose main class `theory.theory` is a sort of enhanced version of `theorydata.getdata()` providing extra data specific to the integral equation calculations.
The method `integralequations.xarray.getCluster()` returns conjectural predictions for the same quantities computed by `framedata.getCluster()`.

* `hkmetric` handles all hyperkaehler metric calculations, both direct and conjectural, and resembles a hybrid of the `framedata` and `integralequations` modules.

* `extsde` is a thin wrapper adding extra features to the module `sde`, which in turn translates a high-level description of the Higgs bundle into a low-level explicit form for its self-duality equation, represented as a non-linear Poisson (NLP) equation $\Delta(u) = F(u,x,y)$ on a rectangular grid.  Thus for example `sde` contains the explicit forms of the Higgs field matrices.  A NLP solver is then called to compute the self-dual metric; two modules implement such equations, and `sde` can use either one:  `nlp_euler` and `nlp_fourier`

* `nlp_euler` implements the `euler` solver described in the paper, which is based on solving the finite difference linearization of the NLP equation.  This method has high memory consumption but converges in relatively few iterations.

* `nlp_fourier` implements the `fourier` spectral method described in the paper.  This method has much lower memory consumption, and if the grid size has the form 2^n-1, each iteration is very fast, but many iterations are required for convergence.  (For grid sizes exceeding 511x511, this method is usually best.)

* `planarode` contains the description of the parallel transport ODE for the flat connections we consider, and calls ODE solvers from `scipy.integrate` to solve them.

## Reproducing experiments from the paper

### Generating the data

The subdirectory `paper-computations` contains scripts that automate computing integral equation and frame data for many theories and parameter values (and storing the results in data files). 

The script functions are as follows:

* `paper-computations-xaroper.py` - oper integral equation computation series (~400 calculations, each taking a few minutes)
* `paper-computations-xarhitchin.py` - Hitchin section integral equation computation series (~400 calculations, each teaking a few minutes)
* `paper-computations-frameoper.py` -  direct (ODE-based) oper computation series (~400 calculations, each taking a few minutes)
* `paper-computations-framehitchin.py` - direct (PDE-based) Higgs bundle computation series for a given grid size `--pde-nmesh` (~400 calculations, each taking a few minutes and a few hundred MiB RAM for `pde_nmesh=1023` or several hours and 14GiB RAM for `pde_nmesh=8191`)
* `paper-computations-hkmetric.py` - direct and integral equation hyperkaehler metric series for a given grid size (25 calculations, running time and memory use highly dependent on `pde_nmesh`)
* `paper-computation.sh` - Run each of the scripts above, and run `paper-computations-framehitchin.py` for each `pde_nmesh` in `[2047,4095,8191]`.

All of the scripts support two modes of operation, controlled by a command-line option:
* `--run` runs the calculations, one at a time
* `--print` prints a list of one-line shell scripts that can be run independently (in parallel on one machine, or across several machines) to perform the calculations.  These scripts should run with the repository root as the working directory.

Thus running the single command
```
./paper-compuations/paper-computations.sh --run
```
in the repository root will reproduce the full dataset associated with the paper (written to files in the directories configured in `paths.conf`), however running the calculations serially will take several months on a fast CPU.  In practice one is more likely to use the `--print` mode to distribute the calculations across several machines, or at least run several in parallel on a single machine.  For example, if GNU parallel is available, the command
```
./paper-computations/paper-computation.sh --print | nice -19 parallel --line-buffer
```
will run as many calculations as there are CPU cores, in parallel, until all are completed.  Some calculations use up to 14 GiB of RAM, so it may also be necessary to limit the number of parallel processes with GNU parallel's `-j` option.

## Producing tables and plots

The subdirectory `paper-figures` contains scripts that read the data files created by `paper-computations.sh` (or the dataset provided with the paper) and produce the plots and tables from the paper using matplotlib.

The script functions are as follows:

* `paperfig-hitchin.py` - generates side-by-side plots of spectral coordinates and relative differences on the Hitchin section for a given theory and range of R values

* `paperfig-oper.py` - generates side-by-side plots of spectral coordinates and relative differences for opers for a given theory and range of hbar values

* `paperfig-hkmetric.py` - generates side-by-side plots of the Hitchin metric tensor and relative differences for the A1A2 theory

* `paperfig-heatmaps.py` - generates a grid of heatmap-style plots of the Hitchin metric integrand

* `paperfig-xaroper.py` - generates a plot of the integrand in the TRH integral equations for A1A2 opers

* `paperfig-xarshitchin.py` - generates a plot of the integrands in the TRH integral equations for A1A2 Hitchin section for several values of R

* `papertable-hitchin.py` - generates a table of spectral coordinates and relative differences for the A1A2 Hitchin section

* `papertable-hitchin.py` - generates a table of spectral coordinates and relative differences for A1A2 opers

The `Makefile` in this directory contains rules to generate all of the figures and tables in the paper.
