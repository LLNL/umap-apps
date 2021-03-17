#!/bin/bash

git clone https://github.com/LivGen/LMAT.git
cd LMAT
git checkout 5fca89e

cp ../CMakeLists.txt .
cp ../run_rl.sh ./bin/run_rl.sh
cp ../redoall.sh redoall.sh
cp ../umap-patch.patch umap-patch.patch

bash redoall.sh U

