#!/bin/bash

git clone https://github.com/jshun/ligra.git
cd ligra
git checkout 7755d95

cp ../ligra-umap.patch ligra-umap.patch

git apply ligra-umap.patch


