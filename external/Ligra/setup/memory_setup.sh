#!/bin/bash

function free_mem {
    m=`grep MemFree /proc/meminfo | awk -v N=2 '{print $N}'`
    fm=$(((${m}/1024)/1024))
    echo ${fm}" GB Free on the system"
}

function avail_mem {
    m=`grep MemAvail /proc/meminfo | awk -v N=2 '{print $N}'`
    fm=$(((${m}/1024)/1024))
    echo ${fm}" GB Free on the system"
}

function drop_page_cache {
    echo "Dropping page cache"
    sudo sh -c "/usr/bin/echo 3 > /proc/sys/vm/drop_caches"
    #/home/perma/drop_buffer_cache
    #srun --drop-caches=pagecache hostname
}

function system_info {
  uname -a
}

function test_setup {
  drop_page_cache
  avail_mem
  numactl -H
  system_info
}

function waste_memory_1 {
    free=`numactl -H | grep "node 1 free:" | awk '{print $4}'`
    echo ${free}" MB free on NUMA node 1"
    waste=$(( free / 1024 - 10 ))
    echo "Wasting $waste GB of memory"
    if [ $waste > 0 ]
    then
	echo dd if=/dev/zero of=/dev/shm/${waste}GB bs=4096 count=$((${waste}*256*1024))
	numactl -C 1 -N 1 dd if=/dev/zero of=/dev/shm/${waste}GB bs=4096 count=$((${waste}*256*1024))
    fi
}
function waste_memory_0 {
    free=`numactl -H | grep "node 0 free:" | awk '{print $4}'`
    echo ${free}" MB free on NUMA node 0"
    waste=$(( free / 1024 - 10 ))
    echo "Wasting $waste GB of memory"
    if [ $waste > 0 ]
    then
	echo dd if=/dev/zero of=/dev/shm/${waste}GB bs=4096 count=$((${waste}*256*1024))
	numactl -C 0 -N 0 dd if=/dev/zero of=/dev/shm/${waste}GB bs=4096 count=$((${waste}*256*1024))
    fi
}

function waste_memory_old {
    m=`grep MemFree /proc/meminfo | awk -v N=2 '{print $N}'`
    fm=$(((${m}/1024)))
    echo ${fm}" MB Free on the system"
    waste=$(( fm / 1024 - 50 ))
    echo "Wasting $waste GB of memory"
    if [ $waste > 0 ]
    then
	echo dd if=/dev/zero of=/dev/shm/${waste}GB bs=4096 count=$((${waste}*256*1024))
	dd if=/dev/zero of=/dev/shm/${waste}GB bs=4096 count=$((${waste}*256*1024))
    fi
}

function amounttowaste {
  m=`grep MemFree /proc/meminfo | awk -v N=2 '{print $N}'`
  fm=$(((${m}/1024)/1024))
  waste=$((${fm}-${memtoleave}))
  echo $fm GB Available, Wasting $waste GB
}

test_setup
#drop_page_cache
#waste_memory_0
#waste_memory_1
exit

umap_lib=/home/peng8/umap/build/lib
ld_path="${umap_lib}:/home/peng8/skylake/Software/perfmon2-libpfm4/lib:/home/peng8/umap/ext/profiling/pebs/lib:${LD_LIBRARY_PATH}"
input_dir=/mnt/pmem/pm0/ligra

declare -a graph_arr=(
    "ligra_rMat_s_n20"
    "ligra_rMat_s_n26"
    "ligra_rMat_s_n27"
    "ligra_rMat_s_n28"
    "ligra_rMat_s_n29"
    "com-friendster"
    "com-orkut"
    "soc-LiveJournal1"
    "soc-pokec"
    "cit-Patents"
)

declare -a source_arr=(
    "0"
    "0"
    "0"
    "0"
    "0"
    "101"
    "1"
    "0"
    "1"
    "3858241"
)

len=5 #${#graph_arr[@]}
#numa_setup=" numactl --membind=1 --cpunodebind=1 "

for i in $(seq 4 $((len-1)))
do
    input="${graph_arr[$i]}"
    root="${source_arr[$i]}"

    for app in BFS #BC Radii Components PageRank
    do
        exe=${app}
        output=${exe}_${input}.log
        test_setup > ${output}
        cmd="env LD_LIBRARY_PATH=${ld_path} ${numa_setup} ./${exe} -maxiters 1 -s -m -r $root ${input_dir}/${input}" #>> ${output}"
        echo $cmd #>> ${output}
        #eval $cmd
        continue

        exe=${app}"_umap"
        for psize in 16384 32768 65536 131072 262144 524288 1048576
        do
            output=numax2_${exe}_${input}_psize${psize}.log
            test_setup > ${output}
            #cmd="env UMAP_PAGESIZE=${psize} LD_LIBRARY_PATH=${ld_path} numactl --membind=1 --cpunodebind=1 ./${exe} -s -m ${input_pref}${s} >> ${output}"
            cmd="env UMAP_PAGESIZE=${psize} LD_LIBRARY_PATH=${ld_path} ./${exe} -maxiters 1 -s -m -r $root ${input_dir}/${input} >> ${output}"
            echo $cmd >> ${output}
            eval $cmd
        done
    done
done

echo "Done"
    
