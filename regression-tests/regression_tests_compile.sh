#!/bin/bash
set -e

##############################################
# Compile Umap and Umap-app
##############################################


############################################## 
# Check if UMAP_ROOT or UMAP_APP_ROOT are not set
##############################################
if [ -z "$UMAP_ROOT" ];
then
    echo "UMAP_ROOT is not set. Try to use the current directory."
    UMAP_ROOT=$(pwd)/umap
    if [ ! -d "$UMAP_ROOT" ];
    then
	echo "Cannnot find umap. Please set UMAP_ROOT."
	exit
    else
	echo "Start Compiling Umap."
	cd $UMAP_ROOT
	mkdir -p build && cd build
	rm -f CMakeCache.txt
	cmake3 -DCMAKE_INSTALL_PREFIX=. ..
	make -j && make install
	cd ../..
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
    else
	echo "Start Compiling Umap-apps."
	cd $UMAP_APP_ROOT
	mkdir -p build && cd build
	rm -f CMakeCache.txt
	cmake3 -DCMAKE_INSTALL_PREFIX=. -DUMAP_INSTALL_PATH=$UMAP_ROOT/build ..
	make -j && make install
	cd ../..
    fi
fi

echo "Complete!"
 
exit


