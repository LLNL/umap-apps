#!/bin/bash


##################################################
# Download Umap and UmapApp from github.com
# By default, the develop branch is downloaded
##################################################

if [ -z "$umap_branch" ];
then
    umap_branch="develop"
fi

if [ -z "$app_branch" ];
then
    app_branch="develop"
fi

echo "Start downloading Umap."
git clone --depth 1 -b $umap_branch https://github.com/LLNL/umap.git 

echo "Start downloading Umap."
git clone --depth 1 -b $app_branch https://github.com/LLNL/umap-apps.git 

echo "Complete downloading!"
 
exit


