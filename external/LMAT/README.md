# Livermore  Metagenomics  Analysis  Toolkit
**Taxonomic classification, content summarization
 and gene identification: all-in-1 metagenomic analysis toolkit**

[https://github.com/LivGen/LMAT]

# UMap [https://github.com/LLNL/umap]

# PERMA [https://github.com/khyox/perm-je]

This page ports PERMA to use UMap for file-backed memory mapping. Transparently, LMAT built atop PERMA will be able to leverage UMap for application-specific tuning of page management.


*The instructions below were tested with LMAT (commit ID 5fca89e), PERMA (commit ID 4c53d09), and UMap (commit ID e784957).


# Build

```bash
# Download LMAT

git clone https://github.com/LivGen/LMAT.git
cd LMAT
git checkout 5fca89e

cp ../CMakeLists.txt .
cp ../run_rl.sh ./bin/run_rl.sh
cp ../redoall.sh redoall.sh
cp ../umap-patch.patch umap-patch.patch

bash redoall.sh U

```

