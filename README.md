# DeepFlame v0.1.0
DeepFlame is a computational fluid dynamics suite for single or multiphase, laminar or turbulent, reacting flows in all speeds with machine learning capabilities. 

## Dependencies
OpenFOAM-7, Cantera C++ lib 2.6.0, Torch C++ lib 1.11.0

## Features
- Native Cantera reader for chemical mechanisms in `.cti`, `.xml` or `.ymal` formats
- Full compatiblity with Cantera's `UnityLewis`, `Mix` and `Multi` transport models
- Zero-dimensional constant pressure or constant volume reactor solver `df0DFoam`
- Pressued-based low-Mach number reacting flow solver `dfLowMachFoam`
- Density-based high-speed reacting flow solver `dfHighSpeedFoam`
- Two-phase Lagrangian/Euler spray reacting flow solver `dfSprayFoam`
- Cantera's native SUNDIALS CVODE solver for chemical reaction rate evaluation
- Torch's tensor operation functionality for neutral network I/O and calculation
- Interface for DNN model to obtain chemical reaction rates 
- Multiple example and tutorial cases with `Allrun` and `Allclean` scripts
  - 0D Perfectly Stirred Reactor
  - 1D Freely Propagating Premixed Flame
  - 2D Lifted Partially Premixed Triple Flame
  - 3D Taylor-Green Vortex with Flame
  - 1D Detotation Wave in Homogeneous Premixed Mixture
  - 3D Aachen Bomb Spray Flame

## How to install
```shell
# source your OpenFOAM

# Note: libcantera does not yet support Arch (i.e. Apple M1 Chip). You can set your libcantera path manually in deepflame/bashrc
conda create -n libcantera
conda activate libcantera
conda install -c cantera libcantera-devel

. install.sh
```

Some compiling issues may happen, try to consider compile your own Cantera and torch C++ libraries, instead of using conda installed Cantera C++ lib and the downloaded torch C++ lib.

## How to use

```shell
# source your OpenFOAM
source deepflame/bashrc
```
cd to examples, execute Allrun

for torch model, please contact ? or read this paper ?

## Citing DeepFlame
If you use DeepFlame for a publication, please use the citation: 

DeepFlame: An OpenFOAM-based CFD suite for multiphase turbulent reacting flows at all speeds with machine learning. URL:https://github.com/deepmodeling/deepflame-dev, 2022.
