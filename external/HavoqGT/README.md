
# Required to Build HavoqGT with UMap

We assume that the following items are already available (installed) on the system:
- GCC 8.1 or more.
- MPI
- CMake 2.6 or more.

Since HavoqGT is not a read-only application,
the write-protect support for userfaultfd() must be supported on the system.

# Build

HavoqGT depends on Boost C++ Libraries and [Metall](https://github.com/LLNL/metall).
Metall also depends on Boost C++ Libraries.
One can install them and UMap using Spack.

An example of building HavoqGT using Spack is:
```bash
spack install umap
spack load umap

# A proper version of Boost C++ Libraries will be installed along with Metall.
spack install metall
spack load metall

git clone https://github.com/LLNL/havoqgt.git
cd havoqgt
git checkout develop
mkdir build_dir
cd build_dir
cmake ../ \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_FLAGS="-std=c++17 -lrt -lstdc++fs -lpthread" \
  -DHAVOQGT_BUILD_TEST=TRUE \
  -DMPIEXEC_NUMPROC_FLAG="-n" \
  -DUSE_UMAP=on
make
make test # option
make install # option
```

Use `CMAKE_CXX_COMPILER=/path/to/g++` and `MPI_CXX_COMPILER=/path/to/mpic++` CMake options to specify a C++ compiler and a MPI compiler, respectively.
To change the install directory, one can use `CMAKE_INSTALL_PREFIX` CMake option.


## Build without Spack

HavoqGT uses header files of Boost C++ Libraries and Metall. One does not need to build them.
On the other hand, one has to build and install UMap.

Here are the CMake variables to tell HavoqGT the locations of Boost C++ Libraries, Metall, and UMap manually.
* `BOOST_ROOT=/path/to/boost`
* `METALL_ROOT=/path/to/metall`
* `UMAP_ROOT=/path/to/umap/install/dir/root`


```bash
# Assume that umap is already available on this system
# To install umap, see https://github.com/LLNL/umap

# Download Boost (Boost C++ Libraries 1.64 or more is required)
wget https://dl.bintray.com/boostorg/release/1.75.0/source/boost_1_75_0.tar.gz
tar xvf boost_1_75_0.tar.gz
export BOOST_ROOT=${PWD}/boost_1_75_0

# Download Metall
git clone git@github.com:LLNL/metall.git
export METALL_ROOT=${PWD}/metall

git clone https://github.com/LLNL/havoqgt.git
cd havoqgt
git checkout develop
mkdir build_dir
cd build_dir
cmake ../ \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_FLAGS="-std=c++17 -lrt -lstdc++fs -lpthread" \
  -DHAVOQGT_BUILD_TEST=TRUE \
  -DMPIEXEC_NUMPROC_FLAG="-n" \
  -DBOOST_ROOT=${BOOST_ROOT} \
  -DMETALL_ROOT=${METALL_ROOT} \
  -DUSE_UMAP=on \
  -DUMAP_ROOT=/path/to/umap/install
make
make test # option
make install # option
```
