Setting up initial conditions
=============================

Setting up the initial state
----------------------------

Let us consider the following system
with the reversible reaction A -> B.

*system.json*:

.. code:: json

  {
  "network" : {
    "species" : [
      {"label" : "A", "density" : 1},
      {"label" : "B", "density" : 1}],
    "reactions" : [
      {"eq" : "A -> B", "k+" : 1, "k-" : 1}]
    }
  }

Once loaded as a system (a DRSystem object),
it has an initial state in which the quantity of each species in each cell
is determined based on the species density and the cell volume.

However, it is also possible to explicitly set the initial state of the system.
This can be done directly in the JSON file/dictionary or with the RDSystem instance.

To set the state explicitly, one must set the *"state"* attribute of the system JSON/dictionary,
the *state* argument of the RDSystem constructor or the
*state* attribute of the RDSystem object.

The *"space*" attribute has been left to default in the example above.
As a consequence, the dimensions of the reaction-diffusion space are w=d=h=1 cell. 
The system state is thus: [quantity of A in cell 0, quantity of B in cell 0].
Let us say we initially want 5 nmol of A and 3 nmol of B in the unique cell of the system.

Using the JSON/dict script, it will be:

.. code:: json

  {
  "network" : {
    "species" : [
      {"label" : "A", "density" : 1},
      {"label" : "B", "density" : 1}],
    "reactions" : [
      {"eq" : "A -> B", "k+" : 1, "k-" : 1}]
    },
  "state" : {"value" : [5, 3], "units" : "µmol"}
  }

.. code:: python

  system = load_rdsystem("system.json")

Using RDSystem object constructor:

.. code:: python

  RDSystem(
    network = RDNetwork(
      species = [
        Species("A", density=1),
        Species("B", density=1)],
      reactions = [
        Reaction("A -> B", kf=1, kr=1)]
      ),
    state = UnitArray([5, 3], "µmol")
    )

Doing it after loading the system:

.. code:: python

  system = load_rdsystem("system.json")
  system.state = UnitArray([5, 3], "µmol")

Doing it after loading the system, using the *set_state* method,
which allows to set the quantity of a given species at a given position:

.. code:: python

  system = load_rdsystem("system.json")
  system.set_state("A", position=(0,0,0), value=5)
  system.set_state("B", position=(0,0,0), value=3)

Setting up the chemostates
--------------------------

Chemostats work the same way as the system state.
Default chemostats are generated based on the chstt attribute of the species,
however, it is also possible to specify explicitly the chemostat distribution.
The only difference is that chemostats are booleans, and thus, do not accept UnitValue or UnitArray objects:

.. code:: json

  {
  "network" : {
    "species" : [
      {"label" : "A", "density" : 1},
      {"label" : "B", "density" : 1}],
    "reactions" : [
      {"eq" : "A -> B", "k+" : 1, "k-" : 1}]
    },
  "chemostats" : [1, 0]
  }

.. code:: python

  system = load_rdsystem("system.json")

Using RDSystem object constructor:

.. code:: python

  RDSystem(
    network = RDNetwork(
      species = [
        Species("A", density=1),
        Species("B", density=1)],
      species = [
        Reaction("A -> B", kf=1, kr=1)]
      ),
    chemostats = [1, 0]
    )

Doing it after loading the system:

.. code:: python

  system = load_rdsystem("system.json")
  system.chemostats = [1, 0]


Doing it after loading the system, using the set_chemostat method:

.. code:: python

  system = load_rdsystem("system.json")
  system.set_chemostat("A", position=(0,0,0), value=5)
  system.set_chemostat("B", position=(0,0,0), value=3)

Changing the Species.density attribute
--------------------------------------

Let us keep the previous example.
The system is loaded, but now we want to change the density of the species A and B.
Just changing A and B's density attributes won't change the system state.
For the change in density to be taken into account, one must also update the default system state:

.. code:: python

  system = load_rdsystem("system.json")
  system.network.get_species("A").density = 10
  system.network.get_species("B").density = 20
  system.set_default_state()
  print(system.state)

It should print : [10., 20.] molecule.

Changing the Species.chstt attribute
----------------------------------------

The Species.chstt attribute states if a chemostat must be
applied for the species.
As for the density attribute, just changing it after the system has been loaded won't
change the actual chemostat map associated with the system.
One must also update the default chemostat maps:

.. code:: python

  system = load_rdsystem("system.json")
  system.network.get_species("A").chstt = 0
  system.network.get_species("B").chstt = 1
  system.set_default_chemostats()
  print(system.chemostats)

It should print : [0 1].
