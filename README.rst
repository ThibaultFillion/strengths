The Strengths Package
=====================

Introduction
------------

Strengths is a modeling and simulation tool for reaction diffusion systems, interfaced as a Python package.
It stands for "Simulation and modeling Tool for Reaction diffusion Networks in Graphs and Tridimensional Heterogenous Systems".
The design of reaction-diffusion systems is facilitated by the use of dictionary as a way to define most of the key objects.

Installing and using Strengths
------------------------------

For Strengths to work properly, a few dependencies have to be installed first, including `NumPy <https://numpy.org>`_ (required) and
`Matplotlib <https://matplotlib.org>`_ (optional). All dependencies can be installed at once with

.. code:: bash

  pip install -r "https://raw.githubusercontent.com/ThibaultFillion/strengths/main/requirements.txt"

The package can then be installed from the `Python Package Index <https://pypi.org>`_ with

.. code:: bash

  pip install strengths

Alternatively, you can build Strengths from source. More details on how to do so are given in the documentation, in the
`"Building Strengths from source" section <https://strengths.readthedocs.io/en/latest/building_strengths_from_source.html>`_.

Once it is successfully installed, you should be able to import it in Python.
Below is a short example, featuring the simulation of the trajectory for a simple
reaction system:

.. code:: python

  # importing strengths

  import strengths as strn
  import strengths.plot as strnplt

  # defining a reaction system:
  # this one have two species, A
  # and B. One can be converted to the other,
  # through the reversible reaction A <=> B.

  system = strn.rdsystem_from_dict({
      "network" : {
          "species" : [
              {"label" : "A", "density" : 150},
              {"label" : "B", "density" : 50}
              ],
          "reactions" : [
              {"eq" : "A -> B", "k+" : 0.01, "k-" : 0.012}
              ]
          }
      })

  # running a simulation
  trajectory = strn.simulate(system, t_sample=list(range(500)))

  # plotting the trajectories of A and B
  strnplt.plot_trajectory(trajectory, ["A", "B"])

More examples, using more advanced features - especially diffusion -, are available in the documentation, in the `"Using Strengths" section <https://strengths.readthedocs.io/en/latest/using_strengths.html>`_.

Documentation
-------------

Detailed information on the package and how to use it are given in the `documentation <https://strengths.readthedocs.io/en/latest/>`_.
Especially, to get started with the package, you may look at the `"Using Strengths" section <https://strengths.readthedocs.io/en/latest/using_strengths.html>`_,
which presents the key functionalities through examples.
For detailed and more exhaustive information on the accessible interface,
please refer to the `"API Reference" section <https://strengths.readthedocs.io/en/latest/apiref.html>`_, where all relevant functions, classes,
methods and attributes are covered in detail.

Source code and contribution
----------------------------

Strengths has a `repository <https://github.com/ThibaultFillion/strengths/tree/main>`_ hosted on GitHub.

If you wish to contribute to the package,
whether by giving your feedback, reporting bugs or errors,
improving the documentation or writing code,
please refer to the project's `community and contribution guidelines <https://github.com/ThibaultFillion/strengths/blob/main/community_and_contribution_guidelines.rst>`_.

Testing
-------

Tests rely on the `pytest <https://docs.pytest.org>`_ package.
To run the unit tests, run either :code:`run_all_tests.py` or :code:`pytest` in the 
`tests` directory.

Licence
-------

Strengths's source code and documentation are licensed under the terms of the `MIT Licence <https://github.com/ThibaultFillion/strengths/blob/main/LICENCE>`_.
You'll find the licence text in your strengths installation root file, or in the root file of the
project's GitHub repository.

Initial authors
---------------

* Thibault Fillion
  * University of Orléans, Avenue du Parc Floral, BP 6749, 45067 Orléans, France
  * Center for Molecular Biophysics (CBM), Rue Charles Sadron CS 80054, 45071 Orléans, France
  * University of Florence, Department of Experimental and Clinical Medicine, Viale Giovanni Batista Morgagni 50, 50314 Firenze, Italy

* Francesco Piazza,
  * University of Florence, Department of Physics & Astronomy, Via Giovanni Sansone 1, 50019 Sesto Fiorentino, Italy
  * Instituto Nazionale di Fisica Nucleare (INFN) sezione di Firenze, Via Giovanni Sansone 1, 50019 Sesto Fiorentino, Italy

All authors
-----------

See GitHub contributors for now. A dedicated Author file
could be set up in the future.
