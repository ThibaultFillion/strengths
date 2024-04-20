Example 2
=========

The folder contains the simulation script for the first example of the paper.
It contains the following files:

for the animal shaped example:
*  README.rst : This file.
*  animal2_env.npy : Environment map for the animal shaped system.
*  animal2_results.npy : Plotting script for the animal shaped system.
*  animal2_script.npy : Simulation script for the animal shaped system.
and output files:
*  animal2_output(0).json
*  animal2_output(0)_data.npy
*  animal2_output(1).json
*  animal2_output(1)_data.npy
*  animal2_output(2).json
*  animal2_output(2)_data.npy
*  animal2_output(detail).json
*  animal2_output(detail)_data.npy
for the simple 1D and 2D layouts:
*  patterns.py : Simulation script for the simple 1D and 2D systems.
*  patterns_results.py : Plotting script for the simple 1D and 2D systems.
and output files:
*  pattern1D_output_1.json
*  pattern1D_output_1_data.npy
*  pattern2D_output_1.json
*  pattern3D_output_1_data.npy
for the sphere example:
*  env.txt : Environment map for the spheric system.
*  cgmap.txt : Coarse-graining map merging unused cells inside and outside the sphere.
*  sim.py : Simulation script for the spheric system.
*  plotenv.py : Script for plotting the layout of the spheric system.
*  plotsim.py : plotting script for the spheric system.
and output files:
*  sim_output_1.json
*  sim_output_1_data.npy

To run this example, the simulation scripts must be ran fist,
the plotting scripts may be ran to display the results.
It is also possible to directly plot the original results witch are already available in the directory.

Note: the original results (NumPy and JSON files) inside the directory were generated with an early version of the package and of the script files.
The scripts was optimized and reorganized a bit compared to the one which generated the original results, including splitting some code in different
files and renaming some files, as well as putting the environment map in text file, rather than generating it from another file.

Reaction and diffusion kinetic parameters are guessed.
