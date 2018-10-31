.. _getting_started:

===============
Getting Started
===============

This page provides information on how to quickly get up and running with
umap-apps.

------------
Installation
------------

Umap-apps is hosted on GitHub `here <https://github.com/LLNL/umap-apps>`_.
To clone the repo into your local working space, type:

.. code-block:: bash

  $ git clone https://github.com/LLNL/umap-apps.git

or

.. code-block:: bash

  $ git clone git@github.com:LLNL/umap-apps.git

^^^^^^^^^^^^^
Dependencies
^^^^^^^^^^^^^
At a minimum, the programs under umap-apps depend upon umap being installed.
Follow the instructions for installing umap from the umap repository located
here <https://github.com/LLNL/umap>`_.

^^^^^^^^^^^^^^^^^^
Building umap-apps
^^^^^^^^^^^^^^^^^^

Umap uses CMake to handle builds. Make sure that you have a modern
compiler loaded and the configuration is as simple as:

.. code-block:: bash

  $ mkdir build && cd build
  $ cmake -DUMAP_LIBRARY_PATH="<path to installed umap libraries>" -DUMAP_INCLUDE_PATH="<path to installed umap include files>" ../

By default, umap-apps will build a Release type build and will use the system
defined directories for installation.  You may specify a different build type
such as Debug and your own specific directory for installation by using the
two "-D" parameters above.  CMake will provide output about which compiler iscw
being used. Once CMake has completed, umap-apps can be built with Make as follows:

.. code-block:: bash

  $ make

For more advanced configuration, see :doc:`advanced_configuration`.

^^^^^^^^^^^^^^^^^^^^
Installing umap-apps
^^^^^^^^^^^^^^^^^^^^

To install umap-apps, just run:

.. code-block:: bash

  $ make install

Umap-apps install files to the ``lib``, ``include`` and ``bin`` directories of the
``CMAKE_INSTALL_PREFIX``. 

