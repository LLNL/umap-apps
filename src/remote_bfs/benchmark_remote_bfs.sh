#!/bin/bash

# -------------------------------------------------------- #
# Usage
# -------------------------------------------------------- #
# cd path/to/build
# sh path/to/umap-apps/src/remote_bfs/benchmark_remote_bfs.sh

# -------------------------------------------------------- #
# Input Configuration
# -------------------------------------------------------- #
EXE="./bin/remote_bfs"
dir="/p/lscratchh/peng8"  #"/l/ssd"

# -------------------------------------------------------- #
# Machine Configuration
# -------------------------------------------------------- #
serverNode=flash1
clientNode=flash1
#"flash3,flash4,flash5,flash6,flash7,flash8,flash9,flash10,flash11,flash12,flash13"

margo_lib_dir="/g/g90/peng8/flash/umap/build/lib"
lpath="LD_LIBRARY_PATH=${margo_lib_dir}:/usr/tce/packages/openmpi/openmpi-4.0.0-gcc-8.1.0/lib:/usr/tce/packages/python/python-3.7.2/lib:/g/g90/peng8/cuda/lib64:/usr/tce/packages/cuda/cuda-10.0.130/lib64:/g/g90/peng8/cuda/lib64:/usr/tce/packages/cuda/cuda-10.0.130/lib64"

# -------------------------------------------------------- #
# constants
# -------------------------------------------------------- #

K=1024
M=$((K*K))
G=$((K*K*K))

# -------------------------------------------------------- #

# -------------------------------------------------------- #
# Functions
# -------------------------------------------------------- #
prepare_input(){
    if test ! -f $graph_file_path ;then
	numThreads=8
	cmd="env $lpath OMP_NUM_THREADS=$numThreads ./bin/generate_edge_list -o $dir/edge_list -v $scale -e $((2 ** scale * 16))"
	echo $cmd
	eval $cmd

	cmd="env $lpath OMP_NUM_THREADS=$numThreads ./bin/ingest_edge_list -g $graph_file_path "
	t=$(( numThreads-1 ))
	for i in $(seq 0 $t)
	do
	    cmd=$cmd" "$dir/edge_list_$i
	done
	echo $cmd
	eval $cmd
    fi
}

print_gcc_version() {
    echo ""
    ret=$(strings $1 | grep "GCC")
    echo ${ret}
    echo ""
}


# -------------------------------------------------------- #
# Run benchmark varying configuration
# -------------------------------------------------------- #
main() {
    # ---- Print some system information ---- #
    print_gcc_version $EXE

    for scale in 27;do
	
	graph_file_path=$dir"/csr_graph_s${scale}"
	out_prefix="bfs_s${scale}"
    
	num_vertices=$(( 2**scale ))
	num_edges=$(( num_vertices* 16 * 2))

	prepare_input

	BASE_OPTIONS="-n${num_vertices} -m${num_edges} -g${graph_file_path}"

	# ---- Start the server  ---- #
	rm -rf serverfile
	cmd="env $lpath srun --nodelist=${serverNode} --ntasks-per-node=1 -N 1 $EXE ${BASE_OPTIONS} -s &"
	echo $cmd
	eval $cmd

	# -----Start the client ----- #
	while [ ! -f serverfile ]; do
            sleep 1
	done

	for numNodes in 1 #{1..12..1}
	do
	    for numProcPerNode in 1 2 4
	    do
		for numOMPThreads in 1 4 8
		do
		    for pSize in $((1*M)) $((256*K)) #$((64*K)) $((16*K)) # $((4*K)) $((1*M)) $((4*M)) #umap_page_size
		    do
		    cmd="env $lpath OMP_NUM_THREADS=$numOMPThreads OMP_SCHEDULE=static UMAP_PAGESIZE=$pSize  srun --nodelist=${clientNode} --ntasks-per-node=$numProcPerNode -N $numNodes  $EXE ${BASE_OPTIONS}"
		    echo $cmd
		    time eval $cmd
		    done
		done
	    done
	done
    done

    pkill remote_bfs
    sleep 3
    echo "Done"
}

main "$@"
