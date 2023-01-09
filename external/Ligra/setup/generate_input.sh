#!/bin/bash


#Creates symmetric rMAT graphs

./rMatGraph -s 1048576 /l/ssd/ligra_rMat_s_n20
./rMatGraph -s 67108864 /l/ssd/ligra_rMat_s_n26
./rMatGraph -s 134217728 /l/ssd/ligra_rMat_s_n27
./rMatGraph -s 268435456 /l/ssd/ligra_rMat_s_n28
./rMatGraph -s 536870912 /l/ssd/ligra_rMat_s_n29
./rMatGraph -s 1073741824 /l/ssd/ligra_rMat_s_n30
./rMatGraph -s 2147483648 /l/ssd/ligra_rMat_s_n31



#Convert real-world graphs into symmetric graphs
#SNAPtoAdj converts a graph in SNAP format and converts it to Ligra's adjacency graph format

# LJ
wget http://snap.stanford.edu/data/soc-LiveJournal1.txt.gz
gunzip soc-LiveJournal1.txt.gz
./SNAPtoAdj -s soc-LiveJournal1.txt /mnt/pmem/pm1/soc-LiveJournal1

# CP
wget http://snap.stanford.edu/data/cit-Patents.txt.gz
gunzip cit-Patents.txt.gz
./SNAPtoAdj -s cit-Patents.txt /mnt/pmem/pm1/cit-Patents

# FT
wget http://snap.stanford.edu/data/bigdata/communities/com-friendster.ungraph.txt.gz
gunzip com-friendster.ungraph.txt.gz
./SNAPtoAdj -s com-friendster.ungraph.txt /mnt/pmem/pm1/com-friendster

# TW
wget https://suitesparse-collection-website.herokuapp.com/MM/SNAP/twitter7.tar.gz
#http://sparse-files.engr.tamu.edu/MM/SNAP/twitter7.tar.gz
tar -xzvf /mnt/pmem/pm0/twitter7.tar.gz
./SNAPtoAdj -s /mnt/pmem/pm0/twitter7.txt /mnt/pmem/pm0/twitter7

