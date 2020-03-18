#!/bin/bash

# -------------------------------------------------------- #
# Usage
# -------------------------------------------------------- #
# cd path/to/build
# sh path/to/umap-apps/src/remote_bfs/benchmark_remote_bfs.sh

# -------------------------------------------------------- #
# Configuration
# -------------------------------------------------------- #
scale=18
EXE="./bin/remote_bfs"
dir="/p/lscratchh/peng8"  #"/l/ssd"
graph_file_path=$dir"/csr_graph_s${scale}"
out_prefix="bfs_s${scale}"
margo_lib_dir="/g/g90/peng8/flash/umap/build/lib"
lpath="LD_LIBRARY_PATH=${margo_lib_dir}:$LD_LIBRARY_PATH"
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

	exit
    fi
}

execute() {
    echo "$@" |& tee -a ${out_file}
    time "$@" |& tee -a ${out_file}
}

used_gcc_version() {
    ret=$(strings $1 | grep "GCC")
    echo ${ret}
}

run() {
    num_vertices=$((2**${scale}))
    num_edges=$((${num_vertices}*32))

    BASE_OPTIONS="-n${num_vertices} -m${num_edges} -g${graph_file_path}"

    # ---- Generate output file name ---- #
    if [ $usemmap -eq 1 ]; then
	MAP_OPTION=" -s " # use system mmap
	out_file="${out_prefix}_m${usemmap}_t${numOMPThreads}_ra${read_ahead_size}.log"
    else
	MAP_OPTION=""
	out_file="${out_prefix}_m${usemmap}_t${numOMPThreads}_f${umap_page_fillers}_e${umap_page_evictors}_h${umap_high_evict}_l${umap_low_evict}_p${umap_page_size}_r${umap_read_ahead}.log"
    fi
    date | tee ${out_file}

    
    # ---- Print some system information ---- #
    local gcc_version=$(used_gcc_version $EXE )
    echo "USED GCC is " ${gcc_version} >> ${out_file}
    echo "" | tee -a ${out_file}

    
    # ---- Set some environmental variables ---- #
    if [ $usemmap -eq 0 ]; then
      export UMAP_PAGESIZE=$umap_page_size
      env | grep "UMAP" |& tee -a ${out_file}
    fi    

    # ---- Run the benchmark ---- #
    execute env $lpath OMP_NUM_THREADS=$numOMPThreads OMP_SCHEDULE=static $EXE ${BASE_OPTIONS} ${MAP_OPTION} &
    #echo "" |& tee -a ${out_file}

    #date | tee -a ${out_file}
}


# -------------------------------------------------------- #
# Run benchmark varying configuration
# -------------------------------------------------------- #
main() {

    prepare_input
    
    # ---- Run benchmark with mmap ---- #
    usemmap=1
    for numOMPThreads in 48 #96
    do
        run
    done

    # ---- Run benchmark with umap ---- #
    K=1024
    M=$((K*K))
    G=$((K*K*K))

    sleep 3
    usemmap=0
    for numOMPThreads in 48 #96
    do
	for umap_page_size in $((256*K)) #$((64*K)) $((16*K)) # $((4*K)) $((1*M)) $((4*M)) #umap_page_size
	do
	    umap_read_ahead=0
            run
	done
    done
}

main "$@"
