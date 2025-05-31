Array indexing
==============

Multidimensional data is most of the time stored as unidimensional arrays.
There are usually dedicated functions helping to access elements by coordinates,
however, it may still be of interest to understand how data are stored in the first place.
This is what this section will do, by explaining how to compute the unidimensional array index
from coordinates.

Let N, w, h and d be respectively the number of species in a system, and the system width, height and depth in cells.

Space cell index
----------------

The index of a cell at coordinates (x,y,z) of a grid space of dimensions (w, h, d) is given by

.. math::

  \textit{space_index} = z w h + y w + x

System state and system chemostats
----------------------------------

The index for a given species at a given position of a grid space is given by

.. math::

  \textit{state_index} = \textit{species_index} \textit{space_size} + \textit{space_index}

State trajectory data
---------------------

Trajectory data are stored in the RDTrajectory.data property.
The index for a given species, in a given sample at a given position is given by

.. math::

  \textit{trajectory_index} = \textit{sample_index} \textit{state_size} + \textit{state_index}

