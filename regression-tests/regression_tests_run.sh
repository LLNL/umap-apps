#!/bin/bash

##################################################
# This test script does not require sudo privilege 
# The tests run BFS, UMapsort, and Churn tests
# with different parameters
##################################################

if [ -z "$UMAP_ROOT" ];
then
    echo "UMAP_ROOT is not set. Try to use the current directory."
    UMAP_ROOT=$(pwd)/umap
    if [ ! -d "$UMAP_ROOT" ];
    then
	echo "Cannnot find umap. Please set UMAP_ROOT."
	exit
    fi
fi


if [ -z "$UMAP_APP_ROOT" ];
then
    echo "UMAP_APP_ROOT is not set. Try to use the current directory."
    UMAP_APP_ROOT=$(pwd)/umap-apps
    if [ ! -d "$UMAP_APP_ROOT" ];
    then
	echo "Cannnot find umap-apps. Please set UMAP_APP_ROOT."
	exit
    fi
fi



BIN_PATH=$UMAP_APP_ROOT/build/bin

SSD_MNT_PATH="/mnt/ssd" 
if [ ! -d $SSD_MNT_PATH ];
then
    echo "$SSD_MNT_PATH does not exit!"
    exit
fi

export LD_LIBRARY_PATH=$UMAP_ROOT/build/lib:$LD_LIBRARY_PATH

echo "##############################################"
echo "# Churn Test "
echo "##############################################"
CHURN_FILE=$SSD_MNT_PATH/regression_test_churn.dat
$UMAP_ROOT/build/bin/churn --directio -f $CHURN_FILE -b 10000 -c 20000 -l 1000 -d 60
/bin/rm -f $CHURN_FILE
echo ""


echo "##############################################"
echo "# BFS "
echo "##############################################"
INPUT_GRAPH=$SSD_MNT_PATH/test_graph
$BIN_PATH/ingest_edge_list -g $INPUT_GRAPH $UMAP_APP_ROOT/src/bfs/data/edge_list_rmat_s10_?_of_4
env OMP_SCHEDULE=static $BIN_PATH/run_bfs -n 1017 -m 32768 -g $INPUT_GRAPH -l $UMAP_APP_ROOT/src/bfs/data/bfs_level_reference
/bin/rm -f $INPUT_GRAPH
echo ""


echo "##############################################"
echo "# UMapSort default"
echo "##############################################"
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
echo ""


echo "##############################################"
echo "# UMapSort out-of-core"
echo "##############################################"
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
echo ""

echo "Regression Tests complete successfully!"
 
exit


