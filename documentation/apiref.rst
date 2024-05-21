API Reference
=============

The API Reference aims to convey an exhaustive and detailed description of the package functions, classes, methods etc.

Useful info
-----------

Here are some things you might find interesting.

.. toctree::
  arrays_and_numbers

Core features
-------------

Core features are those essential for reaction-diffusion system building and simulation.

Directly imported
^^^^^^^^^^^^^^^^^

This is the part of the interface that is (mostly) directly accessible using "import strengths":

.. toctree::
  rdnetwork
  rdspace
  rdsystem
  rdscript
  rdoutput
  simulation
  units
  engine_collection

Not directly imported
^^^^^^^^^^^^^^^^^^^^^

This is the interface that is not directly accessible from "import strengths".
You will need to import those specific modules to access their features:

.. toctree::
  rdengine
  constants
  kinetics

Side features
-------------

Side features cover other aspects such as plots, and their corresponding submodules
must be imported to access them:

.. toctree::
  plot
