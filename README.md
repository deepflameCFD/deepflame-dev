# DeepFlame v0.1.0
DeepFlame is a computational fluid dynamics suite for single or multiphase, laminar or turbulent, reacting flows in all speeds with machine learning capabilities. It aims to provide an open-source platform bringing together the individual strengths of OpenFOAM, Cantera and Torch libraries for machine learning assisted reacting flow simulations. It is also has the scope to incorporate next-generation heterogenous supercomputing and AI acceleration infrustructures such as GPU and FPGAs.  

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
The installation of DeepFlame is simple and requires OpenFOAM-7, libcantera and libtorch. 

1. Install OpenFOAM-7 (if not already installed)
```
Refer to https://openfoam.org/version/7 for detailed installation instructions.
```
2. Source your OpenFOAM at the default path below (or your own path for OpenFOAM bashrc)
```
source $HOME/OpenFOAM/OpenFOAM-7/etc/bashrc 
```
3. Install precompiled libcantera via conda (Note: libcantera does not yet support Arch (e.g. Apple M1 Chip). You can set your libcantera path manually in deepflame-dev/bashrc)
```
conda create -n libcantera
conda activate libcantera
conda install -c cantera libcantera-devel
```
4. Clone the DeepFlame repository
```
git clone https://github.com/deepmodeling/deepflame-dev.git
cd deepflame-dev
```
5. Install precompiled libtorch 
```
wget https://download.pytorch.org/libtorch/cpu/libtorch-cxx11-abi-shared-with-deps-1.11.0%2Bcpu.zip
unzip libtorch-cxx11-abi-shared-with-deps-1.11.0+cpu.zip -d thirdParty
```
6. Install DeepFlame
```
. install.sh
```
Some compiling issues may happen due to system compatability. Instead of using conda installed Cantera C++ lib and the downloaded Torch C++ lib, try to compile your own Cantera and Torch C++ libraries.

## Running DeepFlame examples
1. Source your OpenFOAM, for example:
```
source $HOME/OpenFOAM/OpenFOAM-7/etc/bashrc 
```
2. Source deepflame-dev/bashrc, for example:
```
source $HOME/deepflame-dev/bashrc
```
3. Go to an example case directory, for example:
```
cd $HOME/deepflame-dev
cd examples/zeroD_cubicReactor/H2/cvodeSolver/
./Allrun
```

Note: For the example cases with torchSover, an additional DNN model file in the `.pt` format is required. Please contact the developers if you would like a test run. 


## Citing DeepFlame
If you use DeepFlame for a publication, please use the citation: 

DeepFlame: An OpenFOAM-based CFD suite for multiphase turbulent reacting flows at all speeds with machine learning. URL:https://github.com/deepmodeling/deepflame-dev, 2022.
