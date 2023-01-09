# Memory-mapped Ligra Graph processing


# UMap [https://github.com/LLNL/umap]

# Ligra [https://github.com/jshun/ligra]

##This page ports Ligra to use UMap for file-backed memory mapping..
###1. support for capture graph metadata and input graph
###2. support memory-mapping routines to replace graph I/O routines
###5. built-in tests of graph applications using synthetic and real graphs
###6. support performance tuning on key UMap configurable parameters.

*The instructions below were tested with Ligra (commit ID 7755d95), and UMap (commit ID e784957).

# Build

```bash
# Download Ligra

git clone https://github.com/jshun/ligra.git
cd ligra
git checkout 7755d95

cp ligraumap.patch .

git apply ligra-umap.patch

```

