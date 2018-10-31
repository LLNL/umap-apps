.. _advanced_configuration:

======================
Advanced Configuration
======================

Listed below are the umap-specific options which may be used when configuring
your build directory with cmake.  Some CMake-specific options have also been
added to show how to make additional changes to the build configuration.

.. code-block:: bash

    cmake -DUMAP_LIBRARY_PATH="<path to installed umap libraries>"
    cmake -DUMAP_INCLUDE_PATH="<path to installed umap include files>"

Here is a summary of the configuration options, their default value, and meaning:

      ===========================  ======== ===============================================================================
      Variable                     Default  Meaning
      ===========================  ======== ===============================================================================
      ``UMAP_LIBRARY_PATH``        not set  Location of umap library
      ``UMAP_INCLUDE_PATH``        not set  Location of umap include files
      ``CFITS_LIBRARY_PATH``       not set  Location of cfitsio library
      ``CFITS_INCLUDE_PATH``       not set  Location of cfitsio include files
      ``CMAKE_CXX_COMPILER``       not set  C++ compiler to use
      ``DCMAKE_CC_COMPILER``       not set  C compiler to use
      ===========================  ======== ===============================================================================

These arguments are explained in more detail below:

* ``UMAP_INCLUDE_PATH`` and ``UMAP_LIBRARY_PATH``
  If these are specified, then the applications that use FITS files as the
  backing store for umap() will be built.

* ``CFITS_INCLUDE_PATH`` and ``CFITS_LIBRARY_PATH``
  If these are specified, then the applications that use FITS files as the
  backing store for umap() will be built.
