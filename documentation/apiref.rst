API Reference
=============

The API Reference aims to convey exhaustive and detailed description of the package functions, classes, methods etc.

useful info
-----------

Here are some thing you might find interesting.

.. toctree::
  arrays_and_numbers

core features
-------------

Core features are those essential for reaction diffusion system building and simulation.

directly exposed interface
^^^^^^^^^^^^^^^^^^^^^^^^^^

This is the part of the interface directly accessible using "import strengths"

.. toctree::
  rdnetwork
  rdspace
  rdsystem
  rdscript
  rdoutput
  simulation
  units
  typechecking

indirectly exposed interface
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is the part of the interface that is not directly accessible from "import strengths".
You'll need to import their respective modules in order to access their features

.. toctree::
  rdengine
  engine_collection
  constants

side features
-------------

Side features corver other aspeces such as plotting, and their corresponding submodules
must be imported in order to access them.

.. toctree::
  plot
  kinetics
