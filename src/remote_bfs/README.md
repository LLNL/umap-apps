## Build:

(1) build margo library from https://xgitlab.cels.anl.gov/sds/margo.git

(2) build network-based UMap library from https://github.com/LLNL/umap/tree/remote_region

(3) in umap-app root directory, mkdir build && cd build && cmake3 -DUMAP_INSTALL_PATH=<where UMap is installed> -DMARGO_ROOT=<where Margo is installed> ..

(4) make install


## Run BFS:

(1) start the servers

```bash
./bin/remote_bfs -n [#of vertices] -m [#of edges] -g [/path/to/graph_file] -s &
```

(2) run the bfs clients

```bash
[UMap environment variables] ./bin/remote_bfs -n [#of vertices] -m [#of edges] -g [/path/to/graph_file]
```

## (Optional) Generate Graph 500 Inputs using an R-MAT Generator

(1) generate edge lists and save them to a file system

```bash
./bin/rmat_edge_generator/generate_edge_list -o [/path/to/edge_lists]/edge_list -v [#scale of the graph] -e [#of edges]
````

(2) Ingest Edge List (construct CSR graph)

```bash
./bin/ingest_edge_list -g [/path/to/graph_file]/csr_graph_file [list of the generated edge lists]
```
