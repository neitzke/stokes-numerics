# stokes-numerics: Installation instructions

These instructions use distro packages where possible, and when packages are installed with pip, they are installed system-wide.  They also install stokes-numerics to work with Python 3.

## Ubuntu 18.04+ or recent Debian

All prerequisites can be satisfied by distribution packages.

```
sudo apt update
sudo apt install -y --no-install-recommends python3-scipy python3-numpy \
    python3-pytest python3-matplotlib python3-jsonpickle \
    texlive texlive-latex-extra
# ---- next line optional but recommended so output has version tagging
sudo apt install -y python3-git
# ----
git clone https://github.com/neitzke/stokes-numerics
```

## Ubuntu 16.04

The distribution-packaged matplotlib and scipy are too old to work, but newer versions can be installed with pip.  (The numpy version installed by pip as a dependency of scipy must also be controlled so as to get a version compatible with Python 3.5.2.)

```
sudo apt update
sudo apt install -y --no-install-recommends python3-pytest \
    python3-jsonpickle texlive texlive-latex-extra
# ---- next line optional but recommended so output has version tagging
sudo apt install -y python3-git
# ----
sudo apt build-dep -y python3-matplotlib
sudo apt build-dep -y python3-scipy
sudo apt install -y python3-pip
sudo -H python3 -m pip install numpy==1.13.3 scipy==1.4 matplotlib==2.2.5
git clone https://github.com/neitzke/stokes-numerics
```

## Ubuntu 14.04

Some prerequisites are not packaged for this release, and others are too old.  Several workarounds are needed to install versions old enough to be compatible with python 3.4 but new enough to work with this program.

```
sudo apt-get update
sudo apt-get install -y --no-install-recommends python3-pytest \
    texlive texlive-latex-extra \
    git build-essential
sudo apt-get install -y python3-pip
sudo apt-get build-dep -y python3-matplotlib
sudo apt-get build-dep -y python3-scipy
sudo -H python3 -m pip install --upgrade pip==18.0 setuptools==41.6.0
sudo -H python3 -m pip install jsonpickle matplotlib==2.2.5 \
    scipy==0.19 --ignore-installed six
# ---- next line optional but recommended so output has version tagging
sudo -H python3 -m pip install GitPython
# ----
git clone https://github.com/neitzke/stokes-numerics    
```

## Testing the installation

Run a simple calculation

```
cd stokes-numerics
python3 -c 'import comparisons; print(comparisons.compareClusters("A1A2",R=1,theta=0,scratch=True))'
```

Running the test suite

```
cd stokes-numerics
python3 -m pytest
```
