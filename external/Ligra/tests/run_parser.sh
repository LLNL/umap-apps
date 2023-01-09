#!/bin/bash

for input in ligra_rMat_s_n26 ligra_rMat_s_n27 ligra_rMat_s_n28 ligra_rMat_s_n29 com-friendster com-orkut soc-LiveJournal1 soc-pokec cit-Patents
do
    for app in BFS BC Radii Components PageRank
    do
	exe=${app}"_mmap"
	output=numax2_${exe}_${input}.log
	echo $output
	grep --max-count=1 "Running time " $output | awk '{print $4}'

	exe=${app}"_umap"
	for psize in 16384 32768 65536 131072 262144 524288 1048576
	do
	    output=numax2_${exe}_${input}_psize${psize}.log
	    grep --max-count=1 "Running time " $output | awk '{print $4}'
	done
	echo ""
    done
done

echo "Done"
exit
    

for s in 24 25 26 27 28 29
do
    exe="bfs_mmap"
    #output=numa1_pmem1_cpu1_${exe}_s${s}.log
    output=numax2_${exe}_s${s}.log
    grep --max-count=1 "Running time " $output | awk '{print $4}'

    exe="bfs_umap"
    for psize in 16384 32768 65536 131072 262144 524288 1048576
    do
	#output=numa1_pmem1_cpu1_${exe}_psize${psize}_s${s}.log
	output=numax2_${exe}_psize${psize}_s${s}.log
	grep --max-count=1 "Running time " $output | awk '{print $4}'
    done
    echo ""
done

echo "Done"
    
