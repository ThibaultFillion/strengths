Dictionary/JSON documentation
==============================

This page describes how input dictionaries/JSON
expected in functions such as *load_rdsystem*, *rdsystem_from_dict*,
*load_rdnetwork*, *rdnetwork_from_dict*,
*load_rdgridspace*, *rdgridspace_from_dict*, etc. should be written.

Reaction-diffusion simulation script
------------------------------------

Describes the parameters for the simulation of the trajectory of a
reaction diffusion system

"system" :
^^^^^^^^^^

the reaction diffusion system to be simulated.

* dictionary describing the system.

* path to a JSON file containging sucha dictionary.

"t_sample" :
^^^^^^^^^^^^

times at which the algorithm should try to sample the system,
if sampling_policy = "on_t_sample".

* unit array dictionnary

* array

"time_step" :
^^^^^^^^^^^^^

in-simulation time step to be used by the algorithm (if the algorithm accept it).

* number (ie. 1),
  interpreted as a UnitValue (time)

* string (ie. "1 s")
  interpreted as a UnitValue (time)

"sampling_policy" :
^^^^^^^^^^^^^^^^^^^

string indicating the method of smapling to be used.
accepted values are :

* "on_t_sample"
* "on_iteration"
* "on_interval"
* "no_sampling"

"sampling_interval" :
^^^^^^^^^^^^^^^^^^^^^

in-simulation time step to be used for the sampling.

* number (ie. 1),
  interpreted as a UnitValue (time)

* string (ie. "1 s")
  interpreted as a UnitValue (time)

"t_max" :
^^^^^^^^^

in-simulation time at which the simulation should stop.

* "default" (string)
  t_max will be the last value in t_sample

* number (ie. 1),
  interpreted as a UnitValue (time)

* string (ie. "1 s")
  interpreted as a UnitValue (time)

"rng_seed" :
^^^^^^^^^^^^

integer to be used as a seed by the simulation engine's pseudo random number generators.

* number

* null/None

"units" :
^^^^^^^^^

alias "units_system", "units system", "u".

* "default"
  Apply the default unit system = µm, s, molecules

* "inherit"
  In the specific case of the reaction-diffusion system,
  "inherit" is the same as "default".

* units system dictionary

default : "default"

Reaction-diffusion system
-------------------------

Describes a physical reaction-diffusion system (RDSystem object).
Functions matching such dict/json to RDSystem instance are rdsystem_from_dict and load_rdsystem.

"network" :
^^^^^^^^^^^

alias "rdnetwork".

* reaction-diffusion network json path (ie. "network.json")

* reaction-diffusion network dictionary

"space" :
^^^^^^^^^

alias "rdspace".

* space json path (ie. "space.json")

* space dictionary

* None/null
  Space will be a cell grid built with default parameters,
  except for its units system, which is inherited from the parent
  reaction diffusion system.
  (system.space = RDGridSpace(units_system=system.units_system))

default : None/null

"state" :
^^^^^^^^^

* unit array dictionary

* None/null
  A default system state will be generated,
  based on the species densities.

default: None/null

"chstt_map" :
^^^^^^^^^^^^^

* array of chemostats (0=false, 1=true)
  ie. [0, 0, 1, ..., 1]

* path to a file containing the array of chemostats

* None/null
  A default chemostat map will be generated according
  to the species chemostats.

default : None/null

"units" :
^^^^^^^^^

alias "units_system", "units system", "u".

* "default"
  Apply the default unit system = µm, s, molecules

* "inherit"
  In the specific case of the reaction-diffusion system,
  "inherit" is the same as "default".

* units system dictionary

default : "inherit"

Reaction-diffusion network
--------------------------

Describes a physical reaction-diffusion network (RDNetwork object).
Functions matching such dict/json to RDNetwork instance are rdnetwork_from_dict and load_rdnetwork.

"species" :
^^^^^^^^^^^

* array of species dictionaries.

"reactions" :
^^^^^^^^^^^^^

* array of reaction dictionaries.

default : []

"environments" :
^^^^^^^^^^^^^^^^

alias "env".

* array of environment labels (strings).

"units" :
^^^^^^^^^

alias "units_system", "units system", "u".

* "default"
  Apply the default units system = µm, s, molecules

* "inherit"
  The units system is inherited from the reaction diffusion system.
  If the network is not declared inside a system, the default units
  system is applied (see "default", above).

* units system dictionary

default : "inherit"

species
-------

Describes a chemical species (Species object).
The functions matching such a dictionary to a Species instance is species_from_dict.

"label" :
^^^^^^^^^

alias "l".

* species label (string)

"density" :
^^^^^^^^^^^

alias "concentration", "dens", "conc", "C".

* density numerical value in quentity/space^3 in the network units system (number).
  ie. "density" : 1

* density in quentity/space^3 (string).
  ie. "density" : "1 molecule/µm3"

* dictionary associating environment labels (keys) to either
  densities numerical value in quantity/space^3 in the network units system (number)
  and densities in quantity/space^3 (string).
  The "default" key, if used, will design the species density to be applied in
  environment which are not specified in the dictionary. by default, "default" is 0.
  ie. "density" : {"env1" : 1, "env2" : "1 molecule/µm3", "default" : 0}

default : 0

"D" :
^^^^^

alias "diff_coef", "diff coef", "diffusion_coefficient", "diffusion coefficient".

* diffusion coefficient numerical value in space^2/time in the network units system (number).
  ie. "D" : 1

* diffusion coefficient in space^2/time (string).
  ie. "D" : "1 µm2/s"

* dictionary associating environment labels (keys) to either
  diffusion coefficient numerical values in space^2/time in the network units system (number)
  and densities in space^2/time (string).
  The "default" key, if used, will design the diffusion coefficient to be applied in
  environment which are not specified in the dictionary. by default, "default" is 0.
  ie. "D" : {"env1" : 1, "env2" : "1 µm2/time", "default" : 0.1}

default : 0

"chstt" :
^^^^^^^^^

alias "chemostat".

* boolean value indicating if the species must be globally chemostated :
  true/1/True : the species must be chemostated
  false/0/False : the species must not be chemostated
  ie. "chstt" : True (python)
  ie. "chstt" : true (json)

* dictionary associating environment labels (keys) to boolean chemostate values indicating if the species should be chemostated in the given compartment.
  The "default" key, if used, will design the chemostat boolean to be applied in
  environments which are not specified in the dictionary. by default, "default" is false.
  ie. "chstt" : {"env1" : true, "env2" : false, "default" : true}

default : false

"units" :
^^^^^^^^^

alias "units_system", "units system", "u".

* "default"
  Apply the default units system = µm, s, molecules

* "inherit"
  The units system is inherited from the reaction diffusion system.
  If the network is not declared inside a system, the default units
  system is applied (see "default", above).

* units system dictionary

default : "inherit"

reaction
--------

Describes a chemical species (Reaction object).
The functions matching such a dictionary to a Reaction instance is reaction_from_dict.

"label" :
^^^^^^^^^

alias "l".

* species label (string)

* None/null

default : None/null

"stoechiometry" :
^^^^^^^^^^^^^^^^^

alias "sto", "equation", "eq".

* stoechiometric equation string
  ie. "stoechiometry" : "A + 2 B -> C"

"k+" :
^^^^^^

alias "kf".

* forward reaction rate constant numerical value in the network units system (number).
  units dimensions depend on the reaction substrates stoechiometry.
  ie. "k+" : 1

* forward reaction rate constant in the network units system (number).
  units dimensions must be chosen according to the substrate stoechiometry.
  ie. "k+" : "1 s-1"

default : 0

"k-" :
^^^^^^

alias "kr".

same as k+, except it is the backward reaction rate.

default : 0

"environments" :
^^^^^^^^^^^^^^^^

alias "env".

* array of the environment labels (strings) in which the reaction should happen.

* None/null

default : None/null

"units" :
^^^^^^^^^

alias "units_system", "units system", "u".

* "default"
  Apply the default units system = µm, s, molecules

* "inherit"
  The units system is inherited from the reaction diffusion system.
  If the network is not declared inside a system, the default units
  system is applied (see "default", above).

* units system dictionary

default : "inherit"

Grid space
----------

Describes the discrete space in which the reaction and diffusion happens (RDGridSpace object).
The functions matching such a json/dictionary to a RDGridSpace instance are load_rdspace/rdspace_from_dict.

"w" :
^^^^^

alias "width".

* width of the cell grid (integer)

default : 1

"h" :
^^^^^

alias "height".

* height of the cell grid (integer)

default : 1

"d" :
^^^^^

alias "depth".

* depth of the cell grid (integer)

default : 1

"cell_env" :
^^^^^^^^^^^^

alias "cell_environments".

* None/null

* integer

* array

default : 0

"cell_volume" :
^^^^^^^^^^^^^^^

alias "cell_vol".

* numerical value for the volume of an individual cell in space^3 in the space units system (number).
  ie. "cell_vol" : 1

* volume of an individual cell in space^3 (string).
  ie. "cell_vol" : "1 µm3"

default : 1

"units" :
^^^^^^^^^

alias "units_system", "units system", "u".

* "default"
  Apply the default units system = µm, s, molecules

* "inherit"
  The units system is inherited from the reaction diffusion system.
  If the network is not declared inside a system, the default units
  system is applied (see "default", above).

* units system dictionary

default : "inherit"

Units system
------------

Describes a choice of units for space distance, time distance and quantity of matter (UnitsSystem).
The function matching such a dictionary to a Species instance is unitssystem_from_dict.

"space" :
^^^^^^^^^

* space distance unit (string)

default : "µm" (default space units)

"time" :
^^^^^^^^

* time distance unit (string)

default : "s" (default time units)

"quantity" :
^^^^^^^^^^^^

* unit for the quantity of matter (string)

default : "molecule" (default quantity units)

unit array
----------

Describe an array of physical quantities expressed in the same units (UnitArray).
The function matching such a dictionary to a Species instance is unitarray_from_dict.

"value" :
^^^^^^^^^

* array of numerical values (numbers).

* path to a file containing the array of numerical values (string).

"units" :
^^^^^^^^^

* units in which the values are expressed
  ie. "units" : "molecule/µm/s"
