The fundamental concepts
========================

Object concepts
---------------

Species (Species class)
^^^^^^^^^^^^^^^^^^^^^^^

A chemical species is a type of molecule. A system typically contains multiple instances of any species across its space.
Species can diffuse and be consumed or produced in reactions.

Reaction (Reaction class)
^^^^^^^^^^^^^^^^^^^^^^^^^

A chemical transformation, or reaction, is an event localized in space during which some substrate species molecules are converted to some product species molecules. Reactions are thus represented by their stoechiometric equation :

.. image:: ABC.png
  :align: center

In the example above, this means that when one molecule of species A collide with one molecule of species B, they can form the a complex, which is another molecule of species C (second order reaction), with a reaction rate constant :math:`k_+`, and that the molecule of C can at any moment break down into a molecule of A and a molecule of C, which is the reverse reaction, (first order), whit a rate constant :math:`k_-`.

Network (RDNetwork class)
^^^^^^^^^^^^^^^^^^^^^^^^^

A reaction diffusion network define a set of coupled species and reactions, along with the environments in which they can exist.

Space (RDSpace class)
^^^^^^^^^^^^^^^^^^^^^

A 3D grid of reaction volumes - that we call cells - that constitute the reaction diffusion system space.
The grid is defined by its dimensions (width, height, depth) and the volume of an individual cell.
Each cell can be assigned some environment (reffering to an environment defined in a reaction diffusion network).

System (RDSystem class)
^^^^^^^^^^^^^^^^^^^^^^^

A reaction diffusion system is a system in the physical sense, a system associate a reaction-diffusion network and a space,
and is characterized by a state, corresponding to the distribution of the different network species quantities in each cell of the system space.
A system also defines how its species are chemostated (chemostats are defined species wise and cell wise).

Script (RDScript class)
^^^^^^^^^^^^^^^^^^^^^^^

A reaction diffusion simulation groups together every simulation parameters, that a simulation engine needs to fully carry out a simulation.
It associates a system with the algorithm and sampling parameters.

Engine (RDEngineBase derived class)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Reaction diffusion simulation engines are objects that actually carry out the simulations.
An engine implements some reaction-diffusion algorithm/method.

Trajectory (RDTrajectory class)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A reaction diffusion simulation trajectory is the output returned by a simulation engine once the simulation is complete.
It contains the successive states sampled during the simulations, as well as coontewtual informations such a copy of
the corresponding system or the simulation script.

Other concepts
--------------

Cells
^^^^^^

The cell is the unit component of the system space.
It is a volume of the system space in which species
can transform and diffuse from/to neighbor cells.
In a 3D space, a cell can be seen as a cube with a maximum of 6 neighbors.
For the simulation methods currently implemented, species are considered 
to be distributed homogeneously inside a given cell.

Environments
^^^^^^^^^^^^

An environment is a type of cell.
Some object porperties, such as species initial densities, chemostate and diffusion coefficient,
and reaction possibility, can be set to different values depending on the environment.
This is the feature that enable the creation of inhomogeneous systems.

Chemostats
^^^^^^^^^^

A chemostat refers to a process maintaining a species quantity at a constant level.
Any species can be defined as chemostated, whether globally, only in certain environments, or locally (cell wise).
A chemostated species will have a quantity maintained at constant value during simulations at the location where the chemostat
is defined.
