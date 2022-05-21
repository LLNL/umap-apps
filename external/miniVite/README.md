# Introduction

MiniVite is benchmark in the ECP ExaGraph suite that implements a single phase of the Louvain method for community detection.
This page describes how to build [miniVite](https://github.com/Exa-Graph/miniVite) with UMap.


miniVite has a mode that uses [Metall](https://github.com/LLNL/metall) to store a graph into persistent memory.
Metall has a mode that uses UMap instead of system mmap() internally.
Therefore, to build miniVite with UMap, all one has to do is have Metall use UMap.

The miniVite version that works with Metall comes with a CMake file.
On this instruction page, we use the CMake file to build miniVite with UMap.


# Build
## Required

We assume that the following items are already available (installed) on the system:
- GCC 8.3 or more.
- MPI
- CMake 3.14 or more.

The write-protect feature of userfaultfd() must be available on the system's Linux kernel.

## Build Manually

miniVite uses header files of Boost C++ Libraries and Metall. One does not needs to build them.
On the other hand, one has to build and install UMap.
To install UMap, see [here](https://github.com/LLNL/umap).

Here are the CMake variables to specify the locations of Boost C++ Libraries, Metall, and UMap manually.
* `BOOST_ROOT=/path/to/boost`
* `METALL_ROOT=/path/to/metall`
* `UMAP_ROOT=/path/to/umap/install/dir/root`


```bash
# Download Boost (Boost C++ Libraries 1.64 or more is required)
wget https://dl.bintray.com/boostorg/release/1.75.0/source/boost_1_75_0.tar.gz
tar xvf boost_1_75_0.tar.gz
export BOOST_ROOT=$PWD/boost_1_75_0

# Download Metall
git clone https://github.com/LLNL/metall.git
export METALL_ROOT=${PWD}/metall

# Build miniVite
git clone https://github.com/Exa-Graph/miniVite.git
cd minivite
git checkout metallds2
mkdir build
cd build
cmake ../ \
 -DBOOST_ROOT=${BOOST_ROOT} \
 -DUSE_METALL=ON \
 -DMETALL_ROOT=${METALL_ROOT} \
 -DUSE_UMAP=ON \
 -DUMAP_ROOT=${UMAP_ROOT}
make
```


## Build with Spack

One can use [Spack](https://spack.io/) to install UMap and Metall.
Spack also can set environmental variables properly,
and the CMake file in miniVite uses the variables to locate UMap and Metall.
On the other hand, one still needs to build miniVite manually using CMake.

```bash
# Install and load UMap
spack install umap
spack load umap

# Install and load Metall
spack install metall
spack load metall

# Build miniVite
git clone https://github.com/Exa-Graph/miniVite.git
cd minivite
git checkout metallds2
mkdir build
cd build
cmake ../ \
 -DUSE_METALL=ON \
 -DUSE_UMAP=ON
make
```


# Run miniVite (example)

```bash
# Run community detection and store graph.
# -n : The total number of vertices of the generated graph.
# -g : A path to store the graph.
mpiexec -n 2 ./miniVite -n 64 -g "/mnt/ssd/graph"

# Load graph and run community detection
# (requires the same number of processes as above).
# -i : A path to the stored graph.
mpiexec -n 2 ./miniVite -i "/mnt/ssd/graph"
```

See details [miniVite README](https://github.com/Exa-Graph/miniVite/tree/metallds2).
