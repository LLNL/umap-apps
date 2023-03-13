#!/bin/bash

OUTPUT_DIR=/mnt/ssd

# Creates symmetric rMAT graphs                                                                                                                       
# Attention: ensure [n=number of vertices] located right before filename

./rMatGraph -s 1048576 ${OUTPUT_DIR}/ligra_rMat_s_n20
./rMatGraph -s -m 536870912  33554432  ${OUTPUT_DIR}/ligra_rMat_sym_s25_k16
./rMatGraph -s 536870912 ${OUTPUT_DIR}/ligra_rMat_s_n29
./rMatGraph -s 1073741824 ${OUTPUT_DIR}/ligra_rMat_s_n30
