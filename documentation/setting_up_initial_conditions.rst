Setting up initial conditions
=============================

Setting up the initial state
----------------------------

let us consider the following system
with the reversible reaction A -> B.

*system.json* :

.. code:: json

  {
  "network" : {
    "species" : [
      {"label" : "A", "density" : 1},
      {"label" : "B", "density" : 1}]
    "reactions" : [
      {"eq" : "A -> B", "k+" : 1, "k-" : 1}]
    },
  }

once loaded as a system (a DRSystem object),
it will have an initial state, in which the quantity of each species in each cell
will be determined based on the species density and the cell volume.

However, it is also possible to explicitly set the initial state of the system.
This can be set directly in the JSON/dictionary or done after the system have been loaded.

To set the state explicitly, one must set the "state" attribute of the system JSON/dictionary,
or the state argument of the constructor.

We left the space to default, so w=d=h=1.
the system state will thus be [quantity of A in cell 0, quantity of B in cell 0].
let us say we initially want 5 nmol of A and 3 nmol of B in the unique cell of the system.

using the JSON/dict script, it will be :

.. code:: json

  {
  "network" : {
    "species" : [
      {"label" : "A", "density" : 1},
      {"label" : "B", "density" : 1}]
    "reactions" : [
      {"eq" : "A -> B", "k+" : 1, "k-" : 1}]
    },
  "state" : {"value" : [5, 3], "units" : "µmol"}
  }

.. code:: python

  system = load("system.json")

using RDSystem object constructor :

.. code:: python

  RDSystem(
    network = RDNetwork(
      species = [
        Species("A", density=1),
        Species("B", density=1)
        ]
      species = [
        Reaction("A -> B", kf=1, kr=1)
        ]
      ),
    state = UnitArray([5, 3], "µmol")
    )

doing it after loading the system :

.. code:: python

  system = load_rdsystem("system.json")
  system.state = UnitArray([5, 3], "µmol")

doing it after loading the system, using the set_state method :

.. code:: python

  system = load_rdsystem("system.json")
  system.state.set_state("A", position=(0,0,0), value=5)
  system.state.set_state("B", position=(0,0,0), value=3)

Setting up the chemostates
--------------------------

Chemostats works the same way as the system state.
Default chemostats are generated based on the chstt attribute of the species,
however, it is also possible to specify explicitly the chemostat distribution.
The only difference is that chemostats are booleans, and thus, does not accept UnitValue or UnitArray objects.

.. code:: json

  {
  "network" : {
    "species" : [
      {"label" : "A", "density" : 1},
      {"label" : "B", "density" : 1}]
    "reactions" : [
      {"eq" : "A -> B", "k+" : 1, "k-" : 1}]
    },
  "chemostats" : [1, 0]
  }

.. code:: python

  system = load("system.json")

using RDSystem object constructor :

.. code:: python

  RDSystem(
    network = RDNetwork(
      species = [
        Species("A", density=1),
        Species("B", density=1)
        ]
      species = [
        Reaction("A -> B", kf=1, kr=1)
        ]
      ),
    chemostats = [1, 0]
    )

doing it after loading the system :

.. code:: python

  system = load_rdsystem("system.json")
  system.chemostats = [1, 0]


doing it after loading the system, using the set_chemostat method :

.. code:: python

  system = load_rdsystem("system.json")
  system.state.set_chemostat("A", position=(0,0,0), value=5)
  system.state.set_chemoatst("B", position=(0,0,0), value=3)

changing the species density of a system
----------------------------------------

let us keep the previous example.
We loaded the system as it was, but now want to change the density of A and B.
We load the system and set A and B densities and print the system state :

.. code:: python

  system = load_rdsystem("system.json")
  system.get_species("A").density = 10
  system.get_species("B").density = 20
  print(system.state)

however, it will print : [1., 1.] molecule
For the change to be taken into account, one must regenerate the default system state :

.. code:: python

  system = load_rdsystem("system.json")
  system.get_species("A").density = 10
  system.get_species("B").density = 20
  system.set_default_state()
  print(system.state)

it will now print : [10., 20.] molecule
