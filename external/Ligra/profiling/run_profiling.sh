#!/bin/bash

umap_lib=/home/peng8/skylake/ligra/umap/build_profiler/lib
ld_path="${umap_lib}:${LD_LIBRARY_PATH}"
input_pref=/mnt/pmem/pm1
input=com-friendster

for app in BFS #BC Radii Components PageRank 
do
    exe=${app}_umap
    output=profile_${exe}_umapbuf5GB_24omp_${input}.log
    cmd="env LD_LIBRARY_PATH=${ld_path} OMP_NUM_THREADS=24 UMAP_PAGESIZE=4096 UMAP_BUFSIZE=1310720 ./${exe} -maxiters 1 -s -m -r 101  ${input_pref}/${input} > ${output}"
    echo $cmd
    #eval $cmd
done

echo "Done"
    
