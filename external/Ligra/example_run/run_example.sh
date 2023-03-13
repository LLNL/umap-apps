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


cd ./apps

export UMAP_INSTALL_PATH=<UMAP_INSTALL_PATH>

#This require sudo 
drop_page_cache

#Option 1: run the original version with file I/O
./BFS -s /mnt/ssd/ip/soc-LiveJournal1

#Option 2: run the memory-mapped version with UMap
env UMAP_PAGESIZE=131072 ./BFS -s -umap /mnt/ssd/ip/soc-LiveJournal1

#Option 2: run the memory-mapped version with UMap using different page sizes
env UMAP_PAGESIZE=262144 ./BFS -s -umap /mnt/ssd/ip/soc-LiveJournal1
