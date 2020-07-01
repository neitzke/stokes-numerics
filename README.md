# stokes-numerics

![Alt text](doc/framevis-example.png?raw=true "Spectral network and polygon")

Python programs and modules for computing Stokes data for polynomial opers and the polynomial Hitchin section of the complex plane, and the predictions of the twistorial Riemann-Hilbert conjecture of Gaiotto-Moore-Neitzke and Gaiotto.

Experiments conducted with these programs are reported in the paper _Opers and nonabelian Hodge: numerical studies_ by David Dumas and Andrew Neitzke.

## Getting Started

### Prerequisites

Required:
* git
* python 2.7 or 3.4+
* jsonpickle >= 0.9.2
* numpy >= 1.8.2
* scipy >= 0.19

Recommended:
* GitPython >= 1.0.1 (packaged as python-git or python3-git in Debian and Ubuntu), needed to tag output files with commit hash of the code that generated them

Some optional features require:
* jupyter / ipython notebook (to view [sample-computations.ipynb](./sample-computations.ipynb))
* matplotlib >= 2.2.5 (necessary to generate plots)
* latex (necessary to generate plots)
* dvipng (necessary to save plots as png)
* pytest (to run the test suite)
* Mathematica 11+ (for generating spectral network pictures)
* GNU make (to generate all the plots using the `Makefile`)

### Installing

Ensure the prerequisite python packages are available and clone this repository, e.g.
```
git clone https://github.com/neitzke/stokes-numerics
```

Detailed installation instructions can be found in [doc/INSTALL.md](doc/INSTALL.md).

## Sample usage

### Stokes data of a polynomial Higgs bundle in rank two

For any real numbers R and theta, the harmonic map from the complex plane to H^2 with Hopf differential -4(R exp(-i\*theta))^2(z^3-1)dz^2 has image an ideal pentagon, determined up to isometry by two cross ratio invariants (*spectral coordinates*).

The predicted cross ratio invariants of the pentagon for R=1 and theta=0.1 by the TRH conjecture can be calculated in the python interpreter as follows:
```
>>> import integralequations
>>> xar = integralequations.computeXar(theoryname="A1A2", R=1)
>>> xar.getCluster(theta=0.1)
[-0.00880077224294364, -0.558911091886089]

```
These values appear in Table 6 of the paper as IEQ spectral coordinates at R=1.

The actual cross ratio invariants for the same parameters in the interpreter:
```
>>> import framedata
>>> fd = framedata.computeFrames(theoryname="A1A2", R=1, theta=0.1, pde_nmesh=2047)
>>> fd.getCluster()
[-0.008800831391373946, -0.5589152011265417]
```
Here the parameter `pde_nmesh` controls the grid size, which can be increased for greater accuracy at the expense of longer computation time and higher memory consumption.  (It should be of the form 2^n-1 for maximum efficiency.)  If `pde_nmesh` is increased to 8191 in the command above, it will exactly reproduce the values reported in the DE column of Table 6 of the paper for R=1, though that calculation takes several hours and requires approximately 14GiB of RAM.

The calculations described above can be done outside the interpreter using the scripts:
```
python compute-frames.py A1A2 1.0 --theta 0.1 --pde-nmesh 2047
python compute-xars.py A1A2 1.0 --theta 0.1
```

In these examples, the theory name `"A1A2"` refers to the family of Hopf differentials -4(R exp(-i\*theta))^2(z^3-1)dz^2.

The TRH predictions can only be calculated for certain named families of holomorphic differentials (in `theorydata.py`), but the direct method can be applied to any polynomial.  The script `compute-polygon.py` takes coefficients of a polynomial P(z) and an angle theta and computes the ideal vertices of the harmonic map image corresponding to Hopf differential -4(exp(-i\*theta))^2 P(z).  Thus the examples above compute cross ratios of the vertices given by
```
python ./compute-polygon.py -1 0 0 1 --theta 0.1
```
with output
```
0.9979225495492065,-0.06442503473969946
-0.19469284421110275,0.9808642599325307
0.04884362929768648,0.998806437642965
0.18665367992663442,0.9824257751961954
0.9995683857444064,-0.02937758023257562
```
which are five points on the unit circle.

### Stokes data of polynomial SL(3,C) opers

Compute predicted spectral coordinates associated to cyclic opers of rank 3 with cubic differential (1/2)(z^2-1)dz^2 by TRH conjecture:

```
>>> import integralequations
>>> xar = integralequations.computeXar(theoryname="A2A1", R=1, oper=True)
>>> xar.getCluster()
[(-0.05368179194231135-0.10524332799475357j), (0.5870816874374098+0.8095276970404678j)]
```

Compute numerical approximations of the actual spectral coordinates by ODE methods:
```
>>> import framedata
>>> fd = framedata.computeFrames(theoryname="A2A1", R=1, oper=True)
>>> fd.getCluster()
[(-0.053681791942411926-0.10524332799481517j), (0.5870816874365472+0.8095276970410847j)]
```

### Spectral network and polygon visualization

This example requires the Wolfram Mathematica kernel `math` to be in the search path.

Draw the spectral network and polygon for Hopf differential -4 R^2 P(z)dz^2, R=0.2, P(z) = z^3-1:
```
python framevis.py -1 0 0 1  -R 0.2 -o example.pdf
```
The image at the top of the README shows the output file `example.pdf`.

### More sample calculations

Many more sample calculations can be found in the jupyter/ipython notebook [sample-computations.ipynb](./sample-computations.ipynb).

## Running the tests

The subdirectory `tests` contains `pytest`-compatible test scripts.  To run all of the tests, in the repository root:

```
python -m pytest
```

## Documentation

Mathematical background and conventions are explained in detail in the paper _Opers and nonabelian Hodge: numerical studies_.

A brief (and incomplete) user guide is available at [doc/USERGUIDE.md](doc/USERGUIDE.md).

## Built with

The PDE solver is derived from [blaschke](https://github.com/daviddumas/blaschke) by David Dumas and Michael Wolf.

The spectral network visualization is derived from the Mathematica notebook [swn-plotter.nb](https://arxiv.org/src/1704.01522v1/anc/swn-plotter.nb) included with the paper [_Integral Iterations for Harmonic Maps_ by Andrew Neitzke](https://arxiv.org/abs/1704.01522).

## Version

This is stokes-numerics version 1.0.0.

Version history:
* 1.0.0: June 29, 2020 (this version)

## Authors

* **David Dumas** - [github.com/daviddumas](https://github.com/daviddumas)
* **Andrew Neitzke** - [github.com/neitzke](https://github.com/neitzke)

## Acknowledgement

The authors were supported in part by the US National Science Foundation, through grants NSF DMS 1709877 (DD) and DMS 1711692 (AN).

This material is based upon work supported by the National Science Foundation. Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author and do not necessarily reflect the views of the National Science Foundation.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
