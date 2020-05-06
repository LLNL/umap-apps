#!/bin/bash

EXE="./bin/remote_bfs"

if test -z "$input_dir";then
   input_dir="$(pwd)"
fi

if test -z "$margo_lib_dir";then
    echo "$margo_lib_dir not set!"
    exit
fi

lpath="LD_LIBRARY_PATH=${margo_lib_dir}:$LD_LIBRARY_PATH"

prepare_input(){
    if test ! -f $graph_file_path ;then
	numThreads=8
	cmd="env $lpath OMP_NUM_THREADS=$numThreads ./bin/generate_edge_list -o ${input_dir}/edge_list -v $scale -e $((2 ** scale * 16))"
	echo $cmd
	eval $cmd

	cmd="env $lpath OMP_NUM_THREADS=$numThreads ./bin/ingest_edge_list -g $graph_file_path "
	t=$(( numThreads-1 ))
	for i in $(seq 0 $t)
	do
	    cmd=$cmd" "${input_dir}/edge_list_$i
	done
	echo $cmd
	eval $cmd
    fi
}

main() {

    for scale in 27 28 29 30;do 
	
	graph_file_path=${input_dir}"/csr_graph_s${scale}"
	out_prefix="bfs_s${scale}"
    
	num_vertices=$(( 2**scale ))
	num_edges=$(( num_vertices* 16 * 2))

	prepare_input

	BASE_OPTIONS="-n${num_vertices} -m${num_edges} -g${graph_file_path}"

	# ---- Start the server  ---- #
	for numServerNode in 2;do
	    for numServerProc in 4;do
		rm -rf serverfile
		serverThreads=$(( 24/numServerProc ))
		cmd="env $lpath OMP_NUM_THREADS=$serverThreads UMAP_PAGESIZE=1048576 srun --ntasks-per-node=${numServerProc} -N ${numServerNode} $EXE ${BASE_OPTIONS} -s &"
		echo $cmd
		eval $cmd

		# -----Start the client ----- #
		while [ ! -f serverfile ]; do
		    sleep 10
		done
		sleep 30

		for numClientNodes in 4;do
		    for numClientProcPerNode in 1 ;do
			for k in 128 16;do
			    numOMPThreads=$(( 24/numClientProcPerNode ))
			    pSize=$(( k * 1024 ))
			    cmd="env $lpath OMP_NUM_THREADS=$numOMPThreads OMP_SCHEDULE=static UMAP_PAGESIZE=$pSize srun --ntasks-per-node=$numClientProcPerNode -N $numClientNodes $EXE ${BASE_OPTIONS}"
			    echo $cmd
			    time eval $cmd
			    sleep 120
			done
		    done
		done
		pkill remote_bfs
		sleep 120
	    done
	done
    done
    
    echo "Done"
}

main "$@"
