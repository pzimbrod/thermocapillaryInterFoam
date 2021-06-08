[![DOI](https://zenodo.org/badge/368456876.svg)](https://zenodo.org/badge/latestdoi/368456876)
![GitHub](https://img.shields.io/github/license/pzimbrod/lpbfFoam)
![Travis (.org)](https://img.shields.io/travis/pzimbrod/lpbfFoam)
![foam](https://img.shields.io/badge/made%20for-OpenFOAM-blue)
![2006](https://img.shields.io/badge/Version-2006-green)
![GitHub issues](https://img.shields.io/github/issues/pzimbrod/lpbfFoam)
![GitHub last commit](https://img.shields.io/github/last-commit/pzimbrod/lpbfFoam)

# OpenFOAM custom solver lpbfFoam

Source Code for a custom solver, based von OpenFOAM v2006 solver `icoReactingMultiphaseInterFoam`. Its main intended purpose is to model Powder Bed Fusion additive manufacturing (PBF), hence the name.
This solver was created due to the lack of a universal formulation for surface tension in the standard solvers. It has been shown however that marangoni stresses influence melt pools during laser melting considerably. This aims to provide a formulation that resolves those stresses in an explicit manner.

## Features

* Multiphase, VoF-based solver for n incompressible phases
* Multi-species phase formulations
* Melting, solidification and evaporation models
* Heat input using `laserDTRM`
* Fluid surface tension forces
* Marangoni force

## Theory

The numerical implementation of capillary physics uses the capillary stress tensor formulation:

<img src="https://render.githubusercontent.com/render/math?math=\frac{\partial T_{i,j}}{\partial x_j} = \frac{\partial}{\partial x_j} \left[ \sigma \delta_s \left(\delta_{i,j} - n_i n_j \right)\right]">

Which in a finite volume framework translates to:

<img src="https://render.githubusercontent.com/render/math?math=\delta_s = \left\lvert \frac{\partial \alpha}{\partial x_i} \right\rvert">

<img src="https://render.githubusercontent.com/render/math?math=n_i = \frac{1}{\left\lvert \frac{\partial \alpha}{\partial x_i} \right\rvert} \frac{\partial \alpha}{\partial x_i}">

The gradient of the surface tension coefficient is computed explicitly in the code. Therefore, the model can account for Marangoni forces due to temperature gradients (as long as the coefficient itself is being modelled as a function of temperature using `temperatureDependent`) and concentration gradients, using the species modelling capabilities.

The thory is based on the following references:

> [1] B. Lafaurie, C. Nardone, R. Scardovelli, S. Zaleski, und G. Zanetti, „Modelling Merging and Fragmentation in Multiphase Flows with SURFER“, Journal of Computational Physics, Bd. 113, Nr. 1, S. 134–147, Juli 1994, doi: [10.1006/jcph.1994.1123](http://doi.org/10.1006/jcph.1994.1123).

> [2] D. Gueyffier, J. Li, A. Nadim, R. Scardovelli, und S. Zaleski, „Volume-of-Fluid Interface Tracking with Smoothed Surface Stress Methods for Three-Dimensional Flows“, Journal of Computational Physics, Bd. 152, Nr. 2, S. 423–456, Juli 1999, doi: [10.1006/jcph.1998.6168](http://doi.org/10.1006/jcph.1998.6168).

## Usage

This solver uses a set of proprietary models that come with it. For the original compilation, OpenFOAM ESI in Version 2006 was used. Upwards and downwards compatibility can not be ensured, but should likely be possible with few modifications to the Makefiles.

Compile the source code using the `Allwmake` script. Requires the standard set of compilers used by standard OpenFOAM installations.

## Citing

Please consider citing this code when using it in your projects.
