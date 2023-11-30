Array indexing
==============

Multidimensional data is most of the time stored as unidimensional arrays.
There are usually dedicated functions helping to access elements by coordinates,
however it may still be of interest to understand how data are stored in the first place.
This is what this section will do, by explaining how to compute the unidimensional array index
from coordinates.

let N, w, h and d be respectively the number of species in a system, and the system width, height and depth in cell.

Space cell index
----------------

the index i of a cell in a grid space is given by

.. math::

  i = z * w*h + y*w + x

where x, y and z are the grid dimensions.

System state and system chemostats
----------------------------------

The index i is, for a given species at a given position is given by

.. math::

  i = species_index * space_size + space_index

State trajectory data
---------------------

trajectory data are stored in the RDTrajectory.data property.
The index i is, for a given species, in a given sample at a given position is given by

.. math::

  i = sample_index * n_samples*n_species*space_size + species_index*space_size + space_index
