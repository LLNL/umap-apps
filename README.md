# UMAP Applications

[![Documentation Status](https://readthedocs.org/projects/umap-apps/badge/?version=latest)](https://umap-apps.readthedocs.io/en/develop/?badge=develop)

This repository contains applications and prototypes that use the umap()
library.

## Link to a Network-based UMap handler

(1) build network-based UMap library from https://github.com/LLNL/umap/tree/remote_region

(2) in umap-app root directory, mkdir build && cd build && cmake3 -DUMAP_INSTALL_PATH=<where UMap is installed> -DMARGO_ROOT=<where Margo is installed> ..

(3) make install

## License

- The license is [LGPL](/LICENSE).
- [thirdparty_licenses.md](/thirdparty_licenses.md)

`LLNL-CODE-733797`

## Contact

Primary contact/Lead developer

## Documentation

Both user and code documentation is available [here](http://umap-apps.readthedocs.io/).

If you have build problems, we have comprehensive [build sytem documentation](https://umap-apps.readthedocs.io/en/develop/advanced_configuration.html) too!

- Maya Gokhale (gokhale2@llnl.gov)

Other developers

- Marty McFadden  (mcfadden8@llnl.gov)
