# Memory-mapped Ligra Graph processing


### UMap [https://github.com/LLNL/umap]

### Ligra [https://github.com/jshun/ligra]

#### This page ports Ligra to use UMap for file-backed memory mapping..
#### 1. support for capture graph metadata and input graph
#### 2. support memory-mapping routines to replace graph I/O routines
#### 3. built-in tests of graph applications using synthetic and real graphs
#### 4. support performance tuning on key UMap configurable parameters.

*The instructions below were tested with Ligra (commit ID 7755d95), and UMap (develop branch).

# Build
#### an example script for applying the UMap patch to Ligra is included in ./build.sh
#### Or, following the steps:

```bash
# Download Ligra

git clone https://github.com/jshun/ligra.git
cd ligra
git checkout 7755d95

cp ../ligra-umap.patch ligra-umap.patch
git apply ligra-umap.patch

cd apps
export UMAP_INSTALL_PATH=<path to Umap installation>
make BFS

```

# New Options: 
```bash
# generate datastore files for memory-map through mmap or umap later
  -g 
# memory map with UMap 
  -umap
# memory map with mmap
  -mmap 
 ```
  
# Run

```bash
# Example generating edge and vertex datastores from graph file
./BFS -s -g input_uncompressed_graph_filename

# Example running BFS using memory-mapped option with mmap
./BFS -s -mmap input_uncompressed_graph_filename

# Example running BFS using memory-mapped option with umap
./BFS -s -umap input_uncompressed_graph_filename

```
