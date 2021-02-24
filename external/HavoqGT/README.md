
# Required to Build HavoqGT with UMap

We assume that the following items are already available (installed) on the system:
- GCC 8.1 or more.
- MPI
- CMake 2.6 or more.

[HavoqGT](https://github.com/LLNL/havoqgt) creates, analyzes, modifies, and saves graphs. To support this usage, the write-protect feature of userfaultfd() must be available on the system's Linux kernel.

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
mkdir build_dir
cd build_dir
cmake ../ \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_FLAGS="-std=c++17 -lrt -lstdc++fs -lpthread" \
  -DHAVOQGT_BUILD_TEST=TRUE \
  -DMPIEXEC_NUMPROC_FLAG="-n" \
  -DUSE_UMAP=on
make

export HAVOQGT_BUILD_ROOT=${PWD}

# Option.
# runs two graph algorithms as correctness testing: breadth-first search (BFS) and connected components (CC).
make test
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

export HAVOQGT_BUILD_ROOT=${PWD}

# Option.
# runs two graph algorithms as correctness testing: breadth-first search (BFS) and connected components (CC).
make test
```


# Run

Here is how to run two example graph programs (BFS and CC) in HavoqGT.
HavoqGT uses MPI for distributed-memory communication.
Please note that the actual MPI launch command depends on your environment.

```bash
export GRAPH_PATH=/dev/shm/graph
# Graph is constructed at '/dev/shm/graph'
# -o <path> : Output graph base filename
mpiexec -n 2 ${HAVOQGT_BUILD_ROOT}/src/ingest_edge_list -o ${GRAPH_PATH} /path/to/edge_list/file1 /path/to/edge_list/file2 # List edge list files at the end

# -i <path> : Input graph base filename
# -s <int>  : Source vertex of BFS (Default is 0)
mpiexec -n 2 ${HAVOQGT_BUILD_ROOT}/src/run_bfs -i ${GRAPH_PATH} -s 0
```
