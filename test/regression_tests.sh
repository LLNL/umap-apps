#!/bin/bash
set -e

##################################################
# This test script does not require sudo privilege 
# This test download and compile Umap and Umap-app
# if they are not in the currrent directory
# It then run BFS, UMapsort, and Churn tests
# with different parameters
##################################################

SSD_MNT_PATH="/mnt/ssd" 

##############################################
# Download and compile Umap and Umap-app
# if they are not specified nor exit
##############################################
if [ -z "$UMAP_ROOT" ];
then
    UMAP_ROOT=$(pwd)/umap
    if [ ! -d "$UMAP_ROOT" ];
    then
	git clone --depth 1 -b develop https://github.com/LLNL/umap.git 
    fi
    cd $UMAP_ROOT
    cmake3 -DCMAKE_INSTALL_PREFIX=. .
    make -j
    make install
    cd ..
fi

if [ -z "$UMAP_APP_ROOT" ];
then
    UMAP_APP_ROOT=$(pwd)/umap-apps
    if [ ! -d "$UMAP_APP_ROOT" ];
    then
	git clone --depth 1 -b develop https://github.com/LLNL/umap-apps.git 
    fi
    cd $UMAP_APP_ROOT
    cmake3 -DCMAKE_INSTALL_PREFIX=. -DUMAP_INSTALL_PATH=$UMAP_ROOT .
    make -j
    make install
    cd ..
fi


BIN_PATH=$UMAP_APP_ROOT/bin


##############################################
# BFS
##############################################
INPUT_GRAPH=$SSD_MNT_PATH/test_graph
$BIN_PATH/ingest_edge_list -g $INPUT_GRAPH $UMAP_APP_ROOT/src/bfs/data/edge_list_rmat_s10_?_of_4
$BIN_PATH/test_bfs -n 1017 -m 32768 -g $INPUT_GRAPH -l $UMAP_APP_ROOT/src/bfs/data/bfs_level_reference
/bin/rm -f $INPUT_GRAPH


##############################################
# UMapSort default
##############################################
DATA_SIZE=$(( 16*1024*1024*1024 ))

for t in 144 96 48; do
    for umap_psize in 65536 16384 4096;do
	rm -f ${SSD_MNT_PATH}/sort_perf_data
	cmd="env UMAP_PAGESIZE=$umap_psize $BIN_PATH/umapsort -f ${SSD_MNT_PATH}/sort_perf_data -p $((DATA_SIZE/umap_psize)) -N 1 -t $t"
	date
	echo $cmd
	time sh -c "$cmd"
    done
done


##############################################
# UMapSort out-of-core
##############################################
DATA_SIZE=$(( 16*1024*1024*1024 ))
BUF_SIZE=$((  12*1024*1024*1024 )) 

for t in 144 96 48; do
    for umap_psize in 65536 16384 4096;do
	rm -f ${SSD_MNT_PATH}/sort_perf_data
	cmd="env UMAP_PAGESIZE=$umap_psize UMAP_BUFSIZE=$((BUF_SIZE/umap_psize)) $BIN_PATH/umapsort -f ${SSD_MNT_PATH}/sort_perf_data -p $((DATA_SIZE/umap_psize)) -N 1 -t $t"
	date
	echo $cmd
	time sh -c "$cmd"
    done
done

exit


# Test for median calculation using fits files
tar -xvf ../median_calculation/data/test_fits_files.tar.gz -C /tmp/
test_median_calculation -f /tmp/test_fits_files/asteroid_sim_epoch
/bin/rm -f /tmp/test_fits_files/*
/bin/rmdir  /tmp/test_fits_files


churn --directio -f /tmp/regression_test_churn.dat -b 10000 -c 20000 -l 1000 -d 60 
/bin/rm -f /tmp/regression_test_churn.dat /tmp/regression_test_sort.dat

