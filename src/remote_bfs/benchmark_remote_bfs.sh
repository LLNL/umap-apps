#!/bin/bash

EXE="./bin/remote_bfs"
dir="/p/lscratchh/peng8"

margo_lib_dir="/g/g90/peng8/flash/deps/lib"
umap_lib_dir="/g/g90/peng8/flash/forked_umap_remote_region/build4/lib"
lpath="LD_LIBRARY_PATH=${margo_lib_dir}:${umap_lib_dir}:$LD_LIBRARY_PATH"

K=1024
M=$((K*K))
G=$((K*K*K))

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

    for scale in 27 28 29;do
	
	graph_file_path=$dir"/csr_graph_s${scale}"
	out_prefix="bfs_s${scale}"
    
	num_vertices=$(( 2**scale ))
	num_edges=$(( num_vertices* 16 * 2))

	prepare_input

	BASE_OPTIONS="-n${num_vertices} -m${num_edges} -g${graph_file_path}"

	# ---- Start the server  ---- #
	for numServerNode in 4 3 2 1;do
	    for numServerProc in 1 2 4;do
		rm -rf serverfile
		serverThreads=$(( 24/numServerProc ))
		cmd="env $lpath OMP_NUM_THREADS=$serverThreads UMAP_PAGESIZE=4194304 srun --ntasks-per-node=${numServerProc} -N ${numServerNode} $EXE ${BASE_OPTIONS} -s &"
		echo $cmd
		eval $cmd

		# -----Start the client ----- #
		while [ ! -f serverfile ]; do
		    sleep 3
		done

		for numClientNodes in 4 3 2 1;do
		    for numClientProcPerNode in 1 2 4;do
			numOMPThreads=$(( 24/numClientProcPerNode ))
			for pSize in $((1*M)) $((256*K)) #$((64*K)) $((16*K)) # $((4*K)) $((1*M)) $((4*M)) #umap_page_size
			do
			    cmd="env $lpath OMP_NUM_THREADS=$numOMPThreads OMP_SCHEDULE=static UMAP_PAGESIZE=$pSize srun --ntasks-per-node=$numClientProcPerNode -N $numClientNodes $EXE ${BASE_OPTIONS}"
			    echo $cmd
			    time eval $cmd
			done
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
