.. _umapsort:

===================
UMap Sort Benchmark
===================

The umapsort benchmark uses file-backed memory mapping to sort a large array (potentially out-of-core) in parallel.

This application supports using UMap for user space page management, as well as regular mmap.

Using UMap, umapsort provides the option of using either the default backing store interface, or a sparse multi-file backing store interface.
Please refer to UMap documentations for details about different backing store interfaces.

^^^^^^
Usage:
^^^^^^

* Specify the UMap page size by assigning the page size in bytes to the environment variable $UMAP_PAGESIZE, for example, to set a 16KB page size:

.. code-block:: bash

  $ export UMAP_PAGESIZE=16384

The default page size is 4096 bytes (4KB).

* To display various options:

.. code-block:: bash

  $ ./umapsort --help
  --help                      - This message
  --initonly                  - Initialize file, then stop
  --noinit                    - Use previously initialized file
  --usemmap                   - Use mmap instead of umap
  --use_sparse_store          - Use UMap with a sparse backing store object (Cannot be used with --usemmap)
  --shuffle                   - Shuffle memory accesses (instead of sequential access)
  -p # of pages               - default: 10000000
  -t # of app threads         - default: 2
  -a # pages to access        - default: 0 - access all pages
  -N # of files               - default: 1
  -f [file name]              - backing file name.  Or file basename if multiple files
  -d [directory name]         - backing directory name.  Or dir basename if multiple dirs
  Environment Variable Configuration (command line arguments obsolete):
  UMAP_PAGESIZE                   - currently: 4096 bytes
  UMAP_PAGE_FILLERS               - currently: 96 fillers
  UMAP_PAGE_EVICTORS              - currently: 96 evictors
  UMAP_READ_AHEAD                 - currently: 0 evictors
  UMAP_BUFSIZE                    - currently: 45124056 pages
  UMAP_EVICT_LOW_WATER_THRESHOLD  - currently: 70 percent full
  UMAP_EVICT_HIGH_WATER_THRESHOLD - currently: 90 percent full
  
* Example using UMap with a sparse multi-file store object:

.. code-block::	bash
  
  $ ./umapsort -d /mnt/ssd/umapsort_multi_file_sparsestore -t 24 -p 1024 -N 16 --use_sparse_store

This command will tell umapsort to generate an array os size $UMAP_PAGESIZE*p and map it to files under "/mnt/ssd/umapsort_multi_file_sparsestore".  
The mapped files will be generated dynamically and as needed, meaning that files will only be created when touched for the first time after a page fault.
The array is then initialized and sorted in parallel. 
