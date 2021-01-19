************
Installation
************


Installing the C++ library
==========================

This is a CMake_ based project. The installation procedure is standard. First,
clone the repository. Then, ``cd`` into the root directory of the
vect2rast project. Let
``vect2rast_INSTALL_PREFIX`` be the path to the directory
where vect2rast should be installed::

  $ git clone https://github.com/sbrisard/vect2rast
  $ cd vect2rast
  $ mkdir build
  $ cd build
  $ cmake -DCMAKE_INSTALL_PREFIX=vect2rast_INSTALL_PREFIX ..
  $ cmake --build . --config Release
  $ cmake --install . --config Release

.. note:: The ``--config`` option might not be available, depending on the
   selected generator.

At this point, vect2rast should be installed. You can now
run the tests::

  $ ctest . -C Release

.. note:: Depending on the system, you might need to add
   ``vect2rast_INSTALL_PREFIX`` to your ``PATH`` environment
   variable.


Compiling your first vect2rast program
==========================================================

``cd`` into the ``example`` subdirectory. The provided example program should be
compiled and linked against vect2rast::

  $ mkdir build
  $ cd build
  $ cmake -Dvect2rast_DIR=vect2rast_INSTALL_PREFIX/lib/cmake/vect2rast ..
  $ cmake --build . --config Release

An executable called ``example_vect2rast`` should be present
in the ``build/Release`` subdirectory.


Building the documentation
==========================

The documentation of vect2rast requires Sphinx_. The C++ API
docs are built with Doxygen_ and the Breathe_ extension to Sphinx_.

To build the HTML version of the docs in the ``public`` subdirectory::

  $ cd docs
  $ sphinx-build -b html . ../public

To build the LaTeX version of the docs::

  $ cd docs
  $ make latex


Installing the Python bindings
==============================

To install the vect2rast module, ``cd`` into the
``python`` subdirectory and edit the ``setup.cfg`` file. Set the ``include_dir``
and ``library_dir`` to the appropriate paths. These should be::

  [vect2rast]
  include_dir = ${CMAKE_INSTALL_PREFIX}/include
  library_dir = ${CMAKE_INSTLAL_PREFIX}/lib

Then, issue the following command::

  $ python setup.py install --user

or (if you intend to edit the project)::

  $ python setup.py develop --user

To run the tests with Pytest_::

  $ python -m pytest tests

.. _Breathe: https://breathe.readthedocs.io/
.. _CMake: https://cmake.org/
.. _Doxygen: https://www.doxygen.nl/
.. _Pytest: https://docs.pytest.org/
.. _Sphinx: https://www.sphinx-doc.org/
