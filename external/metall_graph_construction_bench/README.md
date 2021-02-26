
# Graph Construction Benchmark

[Metall](https://github.com/LLNL/metall) contains a shared-memory dynamic graph construction benchmark.
In this page, we describe how to run the benchmark with UMap (Metall can use UMap instead of system mmap() to allocate memory into persistent memory).

In this benchmark, edges are added to a graph data structure to create a graph. More details about the benchmark are written in [this paper (see Section III)](https://www.osti.gov/servlets/purl/1576900).


# Required

We assume that the following items are already available (installed) on the system:
- GCC 8.1 or more.
- CMake 3.10 or more.
- UMap

The write-protect feature of userfaultfd() must be available on the system's Linux kernel.


*The instructions below were tested with Umap (commit ID e784957) and Metall v0.10.*


# Build

```bash
# Download Boost (Boost C++ Libraries 1.64 or more is required)
wget https://dl.bintray.com/boostorg/release/1.75.0/source/boost_1_75_0.tar.gz
tar xvf boost_1_75_0.tar.gz
export BOOST_ROOT=${PWD}/boost_1_75_0


git clone https://github.com/LLNL/metall.git
cd metall
mkdir build
cd build
cmake .. -DBOOST_ROOT=${BOOST_ROOT} -DUMAP_ROOT=/path/to/umap/install/ -DBUILD_BENCH=ON
make

```



# Run

```
# In 'build/bench/adjacency_list' directory
# The benchmark creates a directory '/mnt/ssd/graph' and stores graph data (files) there
bash ../../../bench/adjacency_list/run_bench.sh -d /mnt/ssd/graph

# For more command-line options,
# see the script: bench/adjacency_list/run_bench.sh
```

Initial evaluation results are available [here](http://sc20.supercomputing.org/proceedings/tech_poster/tech_poster_pages/rpost156.html).
