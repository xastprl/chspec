# chisoth: A local model in XSPEC to fit the X-ray spectrum of Astrophysical plasma.

The  “chisoth”  model  calculates  the synthetic photon spectrum from a spectral library, gen-erated by using the atomic database package, CHIANTI v10. The model takes logarithm of the temperature (LogT), abundances (Ax) ofthe elements with Z=2 to Z=30, and the volume emission measure as input parameters. The volume emissionmeasure is implemented as a normalization factor, which is in the units of 10e46/cm3. 

Details of chisoth can be found in the appendix of:
* [Mondal, Biswajit et al. 2021, ApJ.](https://doi.org/10.3847/1538-4357/ac14c1)

## Pre-requisites

[HEASOFT](https://heasarc.gsfc.nasa.gov/docs/software/heasoft/) package along with the [XSPEC](https://heasarc.gsfc.nasa.gov/xanadu/xspec/) is required for the installation of chisoth.

## Installation

The chisoth package can be downloaded from the Github link-
https://github.com/biswajitmb/CHISOTH.git

After downloading it, go to the CHISOTH directory:

```
cd CHISOTH
```
Initialized Heasoft in your environment if it is not globally set in
your environment.

```
heanit
```
Make the setup.sh file executable

```
chmod 777 setup.sh
```

To install chisoth, run the setup file:

```
./setup.sh
```

## Usage

Open XSPEC command prompt and run the following

```
xspec > lmod isoth model_directory
```
Here, `isoth' is the pre-defined package name and the `model_directory' (e.g., home/CHISOTH/chisoth_gen or home/CHISOTH/chisoth_xsm) is
the full path of model directory in your system.

Now the `chisoth' model will be loaded in your xspec environment and you can call/use it as an [XSPEC model](https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/node27.html).
