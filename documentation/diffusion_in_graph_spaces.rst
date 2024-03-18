Diffusion in graph spaces
=========================

This section details how diffusion is handled by the different engines
``euler_engine()``, ``gillespie_engine()`` and ``tauleap_engine()`` with systems using RDGraphSpace.
Graph type reaction-diffusion spaces (``RDGraphSpace``) are characterized by cells which can have different volumes and more than 6 neighbors.
They usually results from the coarse-graining of a grid type space (``RDGridSpace``).

First, the diffusion coefficient, for a given species, between two neighbor cells *i* and *j*, :math:`D_{ij}`, is defined 
strictly as defined by Bernstein (Bernstein, 2005)[1], by approximating the cell shape as a cube

.. math::
	
	D_{ij} = \frac{ h_i+h_j }{ \frac{h_i}{D_i} + \frac{h_j}{D_j} };

where :math:`D_i` and :math:`D_j` are the diffusion coefficents of the species in cells *i* and *j*, respectively
and :math:`h_i` and :math:`h_j` are the square roots of the volumes of cells *i* and *j*, respectively.

Next, the reaction rate constants for the diffusive process are defined from :math:`D_{ij}` based on Bernstein's method once more (Bernstein, 2005)[1],
but with a slight modification so that the exchange surface shared by neighbor cells is taken into account

.. math::
	k_{i,j} = \frac{D_{ij} s_{ij}} {V_i d_{ij}}\\
	k_{j,i} = \frac{D_{ij} s_{ij}} {V_j d_{ij}}

with :math:`k_{i,j}` the rate constant for diffusion from cell *i* to cell *j*, and :math:`k_{j,i}` from  *j* to *i*,
where :math:`s_{ij}` is the contact surface betwenn *i* and *j*, :math:`d_{ij}` the distance beween the center of *i* and *j*,
and :math:`V_i` and :math:`V_j` the volumes of cells *i* and *j*. In the case where all nodes are cubes with the same volume,
the expressions of :math:`k_{i,j}` and :math:`k_{j,i}` simplifies to those given by Bernstein (Bernstein, 2005)[1].

References
==========

* [1] : Bernstein, D. (2005). Simulating mesoscopic reaction-diffusion systems using the Gillespie algorithm.
	Physical Review E, 71(4), Article 041103. https://doi.org/10.1103/PhysRevE.71.041103
