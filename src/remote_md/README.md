CoMD
====

Classical molecular dynamics proxy application.

This is CoMD version 1.1
------------------------

CoMD is a reference implementation of typical classical molecular
dynamics algorithms and workloads.  It is created and maintained by
ExMatEx: Exascale Co-Design Center for Materials in Extreme Environments
(<a href="http://exmatex.org">exmatex.org</a>).  The
code is intended to serve as a vehicle for co-design by allowing
others to extend and/or reimplement it as needed to test performance of 
new architectures, programming models, etc.

To view the generated Doxygen documentation for CoMD, please visit
<a href="http://exmatex.github.io/CoMD/doxygen-mpi/index.html">exmatex.github.io/CoMD/doxygen-mpi/index.html</a>.

To contact the developers of CoMD send email to exmatex-comd@llnl.gov.

-----------------------

Link to UMap network-based datastore:

(1) build margo library from https://xgitlab.cels.anl.gov/sds/margo.git

(2) build network-based UMap library from https://github.com/LLNL/umap/tree/remote_region

(3) cd ./src-mpi && make MARGO_DIR=<where Margo is installed> UMAP_DIR=<where Umap is installed>

-----------------------

Run MPI version:

(1) start the simulation

    env LD_LIBRARY_PATH=<margo lib>:<umap lib>:$LD_LIBRARY_PATH srun -n 8 ./bin/CoMD-mpi -e -i 2 -j 2 -k 2 -x 40 -y 40 -z 40 -N 600 > output 2>&1 &

(2) start the monitoring client

    env LD_LIBRARY_PATH=<margo lib>:<umap lib>:$LD_LIBRARY_PATH srun -n 1 ./src-mpi/remote_diagnostic 8 50