.. _advanced_configuration:

======================
Advanced Configuration
======================

Listed below are the umap-specific options which may be used when configuring
your build directory with cmake.  Some CMake-specific options have also been
added to show how to make additional changes to the build configuration.

.. code-block:: bash

    cmake -DUMAP_INSTALL_PATH="<path to where umap is installed>"

Here is a summary of the configuration options, their default value, and meaning:

      ===========================  ======== ===============================================================================
      Variable                     Default  Meaning
      ===========================  ======== ===============================================================================
      ``UMAP_INSTALL_PATH``        not set  Location of umap
      ``CFITS_LIBRARY_PATH``       not set  Location of cfitsio library
      ``CFITS_INCLUDE_PATH``       not set  Location of cfitsio include files
      ``CMAKE_CXX_COMPILER``       not set  C++ compiler to use
      ``DCMAKE_CC_COMPILER``       not set  C compiler to use
      ===========================  ======== ===============================================================================

These arguments are explained in more detail below:

* ``UMAP_INSTALL_PATH``
  Location of prerequisite umap installation.

* ``CFITS_INCLUDE_PATH`` and ``CFITS_LIBRARY_PATH``
  If these are specified, then the applications that use FITS files as the
  backing store for umap() will be built.
