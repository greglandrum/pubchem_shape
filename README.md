# Wrapper to make using 

## Compilation and installation

Building this requires you to have an builds of the RDKit and Boost::Python available; the easiest way to get these is to set up a conda environment:

```
% conda create -n py312_shape python=3.12 cmake rdkit libboost-devel libboost-python-devel boost-cpp
% conda activate py312_shape
% wget -O $CONDA_PREFIX/include/rdkit/RDBoost/pyint_api.h https://raw.githubusercontent.com/rdkit/rdkit/Release_2024_03/Code/RDBoost/pyint_api.h
% mkdir build
% cd build
% cmake ..
% make -j4
% ctest --output-on-failure
```
**Note** that the horrible bit above with the `wget` should not be necessary after the 2024.03.6 release of the RDKit is available in conda-forge.

