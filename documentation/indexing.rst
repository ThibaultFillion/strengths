Array indexing
==============

Multidimensional data is most of the time stored as unidimensional arrays.
There are usually dedicated functions helping to access elements by coordinates,
however, it may still be of interest to understand how data are stored in the first place.
This section documents how data is organized inside system states and other important arrays for STReNGHTS.

Let N, w, h and d be respectively the number of species in a system, and the system width, height and depth in cells.

Space cell index
----------------

The index of a cell at coordinates (x, y, z) of a grid space of dimensions (w, h, d) is given by

.. math::

  \textit{space_index} = z \times w \times h + y \times w + x

System state and system chemostats
----------------------------------

The index for a given species at a given position of a grid space is given by

.. math::

  \textit{state_index} = \textit{species_index} \times \textit{space_size} + \textit{space_index}

State trajectory data
---------------------

Trajectory data are stored in the RDTrajectory.data property.
The index for a given species, in a given sample at a given position is given by

.. math::

  \textit{trajectory_index} = \textit{sample_index} \times \textit{state_size} + \textit{state_index}

