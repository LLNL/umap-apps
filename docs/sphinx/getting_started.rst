.. _getting_started:

===============
Getting Started
===============

^^^^^^^^^^^^^
Dependencies
^^^^^^^^^^^^^
At a minimum, the programs under umap-apps depend cmake 3.5.1 or greater and
upon umap being installed.

Follow the instructions for installing umap from
the umap repository located
`here <https://llnl-umap.readthedocs.io/en/develop/getting_started.html>`_.

--------------------------------
Umap-apps Build and Installation
--------------------------------
Once the prerequisites are taken care of, the following lines should get you
up and running:

.. code-block:: bash

  $ git clone https://github.com/LLNL/umap-apps.git
  $ mkdir build && cd build
  $ cmake -DCMAKE_INSTALL_PREFIX=<path to install umap-apps> -DUMAP_INSTALL_PATH=<path to where umap is installed> ..
  $ make -j

By default, umap-apps will build a Release type build and will use the system
defined directories for installation.  To specify different build types or
specify alternate installation paths, see the :doc:`advanced_configuration`.

Should you wish to also install the applications in their respective installation
directories, type:

.. code-block:: bash

  $ make install

Umap-apps installs files to the ``lib``, ``include`` and ``bin`` directories
of the ``CMAKE_INSTALL_PREFIX``. 
