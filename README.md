# OpenFOAM custom solver `lpbfFoam`

Source Code for a custom solver, based von OpenFOAM v2006 solver `icoReactingMultiphaseInterFoam`.

Features:

* Multiphase, VoF-based solver for n incompressible phases
* Melting, solidification and evaporation
* Heat input using `laserDTRM`
* Fluid surface tension forces
* Marangoni force

Some models require additional libraries. Use the separately hosted, custom libraries and save them in:
```$WM_PROJECT_USER_DIR/src```

Compile the source code using the `Allwmake` script. Requires the standard set of compilers used by standard OpenFOAM installations.