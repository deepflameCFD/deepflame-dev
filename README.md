# DeepFlame v0.1.0
DeepFlame is a computational fluid dynamics package for single or multiphase, laminar or turbulent, reacting flows in all speeds with machine learning capabilities.  

## Dependencies
OpenFOAM-7, Cantera c++ lib 2.6.0, torch C++ lib 1.11.0

## Features
You can use Cantera or ANN(torch) to calculate multi-step chemitry.
This repository offers solvers for zero-D, low mach flow, high speed flow and spray.
The original OpenFOAM multi-component species are totally removed, replaced by the newly designed CanteraMixture class.

## How to install
```shell
# source your OpenFOAM

# Note: libcantera does not yet support Arch system. You can set your libcantera path manually in deepflame/bashrc
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
