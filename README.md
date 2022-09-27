# chspec: Local models in XSPEC for X-ray spectra from Astrophysical plasma with CHIANTI database

## Introduction

The chspec package includes multiple models for calculating photon spectrum for
different temperature distributions that can be loaded in
[XSPEC](https://heasarc.gsfc.nasa.gov/xanadu/xspec/) spectral fitting package. At
present it includes two models:

* `chisoth`   : Emission from isothermal plasma
* `chgausdem` : Emission from Gaussian distribution of Differential Emission Measure

These models calculate the synthetic photon spectrum from a spectral library, generated
by using the atomic database package [CHIANTI](https://www.chiantidatabase.org/). The
pre-calculated spectra over grid of temperatures is stored in a FITS file.

The isothermal model takes logarithm of the temperature (LogT), abundances (Ax) of the
elements with Z=2 to Z=30, and the volume emission measure as input parameters. The
volume emission measure is implemented as a normalization factor, which is in the units
of 10e46/cm3. For the gaussian DEM model, peak temperature and width of Gaussian (in
log) are input parameters in addition to the abundances. Normalization gives peak
emission measure in units of 10e46/cm3.

Details of chisoth model and the tabulated spectra used in the model is given in the
appendix of [Mondal, Biswajit et al. 2021, ApJ.](https://doi.org/10.3847/1538-4357/ac14c1)

## Pre-requisites

[HEASOFT](https://heasarc.gsfc.nasa.gov/docs/software/heasoft/) package along with
[XSPEC](https://heasarc.gsfc.nasa.gov/xanadu/xspec/) is required for the installation of
chspec package.

## Installation

- Download the chspec package from this repository by selecting _"Code"_ > _"Download
  ZIP"_.

- Unzip and place the resulting chspec directory in desired location.

- Go to this chspec directory:
``` shell
cd chspec
```
Note: If you have not executed the init script for heasoft in your environment, it needs
to be done (usually by executing `heainit` if that alias is set).

- To install chspec package, run the `install.sh` script:
``` shell
./install.sh
```
This will first download the FITS table file using `wget` and save it in the directory
`chspec/tables`. Then, the model will be compiled and installed.

## Usage

Open XSPEC command prompt and run the following:
``` shell
xspec > lmod chspec /path to/chspec/
```
You can also include this in your .xspec/xspec.rc file to load the model on XSPEC
startup automatically.

Now the `chisoth` and `chgausdem` models will be available for use as as an [XSPEC
model](https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/node27.html). They may be
defined like:
``` shell
xspec > model chisoth
xspec > model chgausdem
```
For an example of the use of `chisoth` model in PyXSPEC, see [this Jupyter
notebook](https://github.com/xastprl/xsm-analysis/blob/master/notebooks/ch2xsm_spectralFit_isothermal_demo.ipynb).
