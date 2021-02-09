# miniVite + umap

This page describes how to build [miniVite](https://github.com/Exa-Graph/miniVite) with umap.


## Overview

miniVite has a mode that uses [Metall](https://github.com/LLNL/metall) to store a graph into persistent memory.
Metall has a mode that uses umap instead of system mmap() internally.
Therefore, to build miniVite with umap, all one has to do is have Metall use umap.

The miniVite version that works with Metall comes with a CMake file.
On this instruction page, we use the CMake file to build miniVite with umap. 


## Required
We assume that the following items are already available (installed) on the system:
- GCC 8.1 or more.
- MPI  
- CMake 3.10 or more.
- umap


## Build with Spack

```bash
# Install and load umap
spack install umap
spack load umap

# Install and load Metall
spack install metall
spack load metall

# Build miniVite
git clone git@github.com:Exa-Graph/miniVite.git
cd minivite
git checkout metallds2
mkdir build
cd build
cmake ../ \
 -DUSE_METALL=ON \
 -DUSE_UMAP=ON \
make
```

Use `CMAKE_CXX_COMPILER=/path/to/g++` and `MPI_CXX_COMPILER=/path/to/mpic++` CMake options to specify a C++ compiler and a MPI compiler, respectively.


## Build without Spack

### Example

Here are the CMake variables to specify the locations of Boost C++ Libraries, Metall, and umap manually.
* `BOOST_ROOT=/path/to/boost`
* `METALL_ROOT=/path/to/metall`
* `UMAP_ROOT=/path/to/umap/install/dir/root`

miniVite uses header files of Boost C++ Libraries and Metall. One does not need to build them.
On the other hand, one has to build and install umap.


```bash
# Download Boost (Boost C++ Libraries 1.64 or more is required)
wget https://dl.bintray.com/boostorg/release/1.75.0/source/boost_1_75_0.tar.gz
tar xvf boost_1_75_0.tar.gz
export BOOST_ROOT=$PWD/boost_1_75_0

# Download Metall
git clone git@github.com:LLNL/metall.git
export METALL_ROOT=${PWD}/metall

# Build miniVite
git clone git@github.com:Exa-Graph/miniVite.git
cd minivite
git checkout metallds2
mkdir build
cd build
cmake ../ \
 -DBOOST_ROOT=${BOOST_ROOT} \
 -DUSE_METALL=ON \
 -DMETALL_ROOT=${METALL_ROOT} \
 -DUSE_UMAP=ON \
 -DUMAP_ROOT=/path/to/umap/install
make
```


## Run miniVite (example)

```bash
# Run community detection and store graph:
mpiexec -n 2 ./miniVite -n 64 -s "/tmp"

# Load graph and run community detection
# (requires the same number of processes as above):
mpiexec -n 2 ./miniVite -c "/tmp"
```