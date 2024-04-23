Example 2
=========

The folder contains the simulation script for the second example of the paper.
It contains the following files:

for the animal shaped example:

*  README.rst : This file.
*  animal2_env.npy : Environment map for the animal shaped system.
*  animal2_results.npy : Plotting script for the animal shaped system.
*  animal2_script.npy : Simulation script for the animal shaped system.

for the simple 1D and 2D layouts:

*  patterns.py : Simulation script for the simple 1D and 2D systems.
*  patterns_results.py : Plotting script for the simple 1D and 2D systems.

for the sphere example:

*  env.txt : Environment map for the spheric system.
*  cgmap.txt : Coarse-graining map merging unused cells inside and outside the sphere.
*  sim.py : Simulation script for the spheric system.
*  plotenv.py : Script for plotting the layout of the spheric system.
*  plotsim.py : plotting script for the spheric system.

To run this example, the simulation scripts must be ran fist,
the plotting scripts may then be ran to display the results.

Note: The original results (featured in the paper) were generated with an early version of the package and of those script files.
The scripts was optimized and reorganized a bit compared to those which generated the original results, including splitting some code in different
files and renaming some files, as well as putting the environment map in text file, rather than generating it from another file.

Reaction and diffusion kinetic parameters are guessed.

Similar reaction-diffusion networks were used in the documentation [1]. 

References
----------

* [1] Fillion, T., & Piazza, F. (2024). Building a reaction diffusion system and performing simulations. *Documentation for Strengths*. https://strengths.readthedocs.io/en/latest/index.html Accessed online the 23th of april, 2024.
