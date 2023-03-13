#!/bin/bash


OUTPUT_DIR=/mnt/ssd

#Convert real-world graphs into symmetric graphs
#SNAPtoAdj converts a graph in SNAP format and converts it to Ligra's adjacency graph format

# LJ
wget http://snap.stanford.edu/data/soc-LiveJournal1.txt.gz
gunzip soc-LiveJournal1.txt.gz
./SNAPtoAdj -s soc-LiveJournal1.txt ${OUTPUT_DIR}/soc-LiveJournal1

# CP
wget http://snap.stanford.edu/data/cit-Patents.txt.gz
gunzip cit-Patents.txt.gz
./SNAPtoAdj -s cit-Patents.txt ${OUTPUT_DIR}/cit-Patents

# FT
wget http://snap.stanford.edu/data/bigdata/communities/com-friendster.ungraph.txt.gz
gunzip com-friendster.ungraph.txt.gz
./SNAPtoAdj -s com-friendster.ungraph.txt ${OUTPUT_DIR}/com-friendster

# TW
wget https://suitesparse-collection-website.herokuapp.com/MM/SNAP/twitter7.tar.gz
#http://sparse-files.engr.tamu.edu/MM/SNAP/twitter7.tar.gz
tar -xzvf twitter7.tar.gz
./SNAPtoAdj -s twitter7.txt ${OUTPUT_DIR}/twitter7

