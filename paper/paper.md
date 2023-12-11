---
title : STReNGTHS, a python package to model and simulate reaction-diffusion systems
authors :
  - name : Thibault Fillion
    affiliation : "1, 2"
  - name : Francesco Piazza
    affiliation : "3, 4"
affiliations :
  - name : Université d'Orléans, UFR Sciences Techniques, Avenue du Parc Floral, BP 6749, 45067 Orléans, France
    index : 1
  - name : Centre de Biophysique Moléculaire (CBM), CNRS-UPR 4301, Rue Charles Sadron, 45071 Orléans, France
    index : 2
  - name : Dipartimento di fisica & Astronomia, Università di Firenze, Via G. Sansone 1, 50019 Sesto Fiorentino, Italy
    index : 3
  - name : INFN sezione di Firenze, Via G. Sansone 1, 50019 Sesto Fiorentino, Italy
    index : 4
bibliography : references.bib
---
# Summary

STReNGTHS is an open-source Python package that provides a simple and intuitive
interface for designing models of discrete 3D heterogeneous reaction-diffusion systems and simulating their trajectories. Different algorithms are available, both stochastic (exact or approximate  solutions of the associated 
master equation) or deterministic (numerical solutions of the corresponding rate equations). 
The acronym stands for "Simulation and modeling Tool for REaction-diffusion Networks in Graphs and Tridimensional Heterogeneous Systems" (STReNGTHS). The simulation algorithms are interfaced through a general abstract interface, which makes it easy to extend STReNGTHS with new algorithms and other features.  It is implemented in python (standard library, 
Numpy [@harris_numpy_2020] and Matplotlib [@hunter_matplotlib_2007]) 
and C++ (standard C++11 or later), and can be easily installed from the Python Package Index (PyPI, https://pypi.org) with (i.e.)

\begin{verbatim}
    pip install strengths
\end{verbatim}

# Statement of need

Biology at the cell scale relies on complex networks of biochemical reactions, 
usually operating across multiple membrane-associated as well as
membrane-less compartments (i.e. plasma membrane, cytoplasm, cellular
organelles, stress granules, etc) and driven far from equilibrium by highly regulated species, 
such as nucleotides (ATP/GTP), amino acids and different 
ions. In order to understand the properties of such reaction-diffusion 
networks, and especially their role in macroscopic *emergent* 
phenomena, convenient, reliable and efficient modeling and simulation tools supporting heterogeneous systems are necessary. Moreover, the choice of the simulation algorithm to be used may depend on the system : deterministic approaches such as ODE integration (rate equations) are effective, but inappropriate for systems that are sensitive to fluctuations, 
such as systems that operate with species present at low copy numbers 
(e.g. certain enzymes, many mRNA species) and/or
in tiny reaction volumes (i.e. within mitochondria), 
for which stochastic methods should be preferred. 
This is why one needs to have multiple, flexible approaches available.
STReNGTHS provides an interface to simulate reaction-diffusion systems 
and manipulate their trajectories, as well as full control
and access to the simulation algorithms themselves.

Reaction-diffusion systems consist primarily 
of a reaction-diffusion network, which is represented by 
the *RDSystem* class. This defines a set of coupled chemical transformations 
that is supplemented with a spatial distribution of 
chemical species, whose features describe
where and how fast each molecule shall diffuse. 
The system space is discrete, and consists of a 3D mesh of individual 
volume elements, which we refer to as *cells*. 
It can be either a regular grid of cubic cells with 
uniform volumes, or an arbitrary network of cells with different volumes, 
which can be obtained by coarse-graining a mesh grid.

In order to account for systems with different compartments, 
SReNGTHS implements the system of reaction-diffusion environments, 
which allow the user to define different types of cells (referred to as environments) 
with specific reactive and diffusive properties. Many 
properties, such as the initial density of species, diffusion coefficients, or reactions 
occurrence, can be defined environment-wise.

Importantly, species can be *chemostatted*, i.e. kept at 
a fixed, prescribed concentration during the simulation, 
globally or only in specific 
environments or cells. Chemostatted species allow one to model 
non-equilibrium conditions existing in living cells that are associated with 
chemical potential baths, such as in the case of the tightly regulated
cytoplasmic levels of ATP or ADP.

![Definition of a simple reaction-diffusion system implementing an 
association reaction over 3 cells using the JSON/dictionary format.
The rates used were: $k_+ = 1$ M$^{-1}$s$^{-1}$, $k_-= 1$ s$^{-1}$. \label{jsonsyntax} ](jsonsyntax.png)

In STReNGTHS reaction-diffusion systems can be defined either 
using python dictionaries or through JSON input files, 
following a specific intuitive syntax, as shown figure \ref{jsonsyntax}.
Simulations are handled by objects called *simulation engines*, 
which offer a general abstract interface for simulation algorithms. 
The *simulate* function warps the engine call to run the whole 
simulation at once. The resulting system trajectory, which is the sequence of 
system states successively sampled during the simulation and the corresponding 
sampling times, is stored in a *RDTrajectory* object. 

So far, STReNGHTS proposes simulation engines implementing the Original 
Gillespie algorithm [@gillespie_exact_1977], the $\tau$-leap approximation to 
the Gillespie algorithm [@gillespie_approximate_2001], and the Euler Method, 
operating both on grid and graph spaces, 
with diffusion handled according to the method described in 
Ref. [@bernstein_simulating_2005]. 

# STReNGTHS and similar tools

There already exist various computational tools to model and simulate
complex reaction-diffusion systems. Some of them are general-purpose 
tools, while others have been designed to handle more specific systems. 
Existing simulation packages include:

- STEPS [@10.3389/neuro.11.015.2009].
  This is a reaction-diffusion program interfaced with python, which 
  uses Gillespie algorithm. It handles simulations in geometries composed of 
  tetrahedral voxels with faces that can represent biological 
  membranes [@10.3389/neuro.11.015.2009].
- Readdy [@hoffmann2019readdy].
  This  is a reaction-diffusion tool with a python interface that uses
  a particle-based approach. An especially interesting feature is that 
  it can deal with complex molecule geometries and reaction patterns, 
  such as polymer dynamics [@hoffmann2019readdy]. 
  A python interface is also available [@hoffmann2019readdy]. 
- MesoRD [@10.1093/bioinformatics/bti431].
  This is a tool that employs a stochastic approach based on the 
  Next Subvolume Method [@Elf2004]. The simulation 
  parameters are defined through XML script files, 
  using the System Biology Markup Language (SBML) 
  format [@10.1093/bioinformatics/bti431]. 
  The software relies on Constructive Solid Geometry (CSG) to define the 
  different reaction-diffusion compartments [@10.1093/bioinformatics/bti431]. 
  It comes with graphical (Windows) 
  and command-line (Unix) user interfaces [@10.1093/bioinformatics/bti431].
- BioNetGen [@10.1093/bioinformatics/btw469].
  This is a modeling and simulation tool 
  that provides a rich scripting language with a 
  rule-based approach [@10.1093/bioinformatics/btw469]. Such an approach enables one 
  to consider systems that may be difficult to apprehend with methods requiring to 
  define explicitly the full reaction network 
  (such as polymerization) [@10.1093/bioinformatics/btw469]. 
  BioNetGen supports both deterministic and stochastic 
  methods [@10.1093/bioinformatics/btw469].

Compared to the aforementioned tools, STReNGTHS is more rudimentary and only 
handles reaction networks with explicitly defined species and reactions, as opposed to 
pattern-based approaches or rule-based approaches used by tools such as 
BioNetGen [@10.1093/bioinformatics/btw469].
Still, it allows one to build in a very intuitive and 
user-friendly way simulations able to describe a vast range of complex systems. 

So far, as opposed with Readdy, STReNGTHS only implements non particle-based methods, 
similar to those proposed by the other software. However, rather than proposing only one 
all-purpose fitting method, STReNGTHS's approach is to display a collection of various 
simulation methods, leaving the choice at the user's discretion. Moreover, simulation 
features can be easily extended using the simulation engine interface. 
In fact, STReNGTHS has been designed to be extended easily.

One of STReNGTHS's key features is the use of *reaction-diffusion environments*, which make 
it easy to design extremely rich system landscapes, i.e. featuring plenty of different 
compartments of arbitrary shape that encode physical and chemical segregation. 

The use of a JSON/dictionary syntax for the definition of reaction-diffusion systems 
brings readability and simplicity to the workflow.

![Example of simulation with STReNGTHS implementing a simple
model  of signal transduction by a single cell. 
(a) Schematic representation of the system. 
(b) Layout of the two different system spaces used,
the $26 \times 26$ mesh grid (left) and its coarse-grained graph version (right). 
(c) Set of coupled stoechiometric equations that compose the reaction network. 
(d) Diffusion coefficients of individual species in the different environments. 
(e) Initial densities of each species in the different environments. 
(f) Time course of the transduced signal: Global trajectory of $Y$ obtained 
from the simulation using the $\tau$-leap 
algorithm [@gillespie_approximate_2001] and the Euler method using both 
system spaces (see b). 
(g) Distribution of the $Y$ species at different times from the $\tau$-leap simulations.
The rates used were: 
$L+R\rightarrow C$: $k_1 = 0.5$ nM$^{-1}$s$^{-1}$. 
$C\rightarrow L + R$: $k_{-1}=0.5 \times 10^{-3}$ s$^{-1}$.
$C\rightarrow R$: $k_2= 10^{-2}$ s$^{-1}$.
$C+X\rightarrow C+Y$: $k_3=1 \ $\ $\mu$M$^{-1}$s$^{-1}$.
$C+Y\rightarrow C+X$: $k_{-3} = 10^{-4}$ $\mu$M$^{-1}$s$^{-1}$.
$Y\rightarrow X, \ k_4 = 10^{-2}$ s$^{-1}$.
The diffusion coefficients for the different species were:
$D_L = 100$ $\mu$m$^2$s$^{-1}$,
$D_R = 0.1$ $\mu$m$^2$s$^{-1}$,
$D_C = 0.1$ $\mu$m$^2$s$^{-1}$,
$D_X = 10$ $\mu$m$^2$s$^{-1}$,
$D_Y = 10$ $\mu$m$^2$s$^{-1}$.
\label{example1}](example1.png)

# Examples

For the first example, let us consider a simple model of signal transduction, 
where some extracellular chemical signal is sensed by a cell, which triggers 
the production of a second messenger, as well as the scavenging of the signal species. 
The network is designed as follows : The extracellular ligand $L$ can bind to a 
plasma membrane receptor $R$ to form a complex $C$. This species catalyzes the 
conversion of the inactive second messenger $X$ to its active form,  $Y$.
The complex $C$ is directly converted back into $R$, which accounts for 
the internalization of the complex, the degradation of the ligand and full 
recycling of the receptor. No degradation of the receptors (either through the proteasomes
or lysosomes) is assumed for simplicity,  although it would be straightforward 
to add additional reactions to implement such reaction channels 
(see Fig. \ref{example1} a,c). 
The model features 3 reaction-diffusion environments, "ext", "cyt" and "mmb", 
accounting, respectively, for the extracellular space, the cytoplasm 
and the interface between these two compartments containing the plasma membrane. 
We use two different system spaces. The first one is a $26 \times 26$ mesh grid consisting of 
1 $\mu$m$^3$ cubic cells, while the second one is a coarse-grained version of the 
first that contains only 85 nodes/cells (for 676 cells in the grid) (Fig. \ref{example1} b). 
The trajectory of the system state is simulated, both for the fully detailed grid and for 
its coarse-grained version, using the $\tau$-leap algorithm [@gillespie_approximate_2001] 
and the Euler method, for a total duration of 1 hour using a time step of 1 ms. 
The global trajectory of $Y$ as well as its distribution at $t=0, 100, 1500$ s 
are plotted in Fig. \ref{example1} (f, g).

![Example of simulations of different pattern-forming reaction-diffusion systems
at increasing level of environmental complexity. 
All simulations are performed with the $\tau$ leap algorithm [@gillespie_approximate_2001].
(a) Description of the chemical reactions and associated rates. 
(c) Evolution of a 1D system in time and space. 
(d) Pseudo-stationary state of a 2D system. 
(e) reaction-diffusion environments for the two inhomogeneous systems. 
(b) Layout of a 3D system where the patterns are forming at the surface of a sphere. 
(f) States of the 3D system described at different time points along a stochastic 
trajectory (concentration of the $A$ species). 
(g) Layout of a 2D system mimicking the shape of an animal.
(h) States of the system (distribution of the $A$ species) at different time points from one simulation.
(i) Patterns formed at t = 3500.0 h resulting from 3 different stochastic simulations.
The rates used were: 
$\emptyset\rightarrow A$: 
$k_1 = 10^{-4}$ molecules$\times\mu$m$^{-3}\times$h$^{-1}$ in $a$ and $c$ and $1.05\times10^{-4}$ molecules$\times\mu$m$^{-3}\times$h$^{-1}$ in $b$.
$A \rightarrow\emptyset$: $k_{-1} = 10^{-3}$ h$^{-1}$.
$\emptyset\rightarrow B$: 
$k_2 = 10^{-4}$ molecules$\times\mu$m$^{-3}\times$h$^{-1}$ in $a$ and $b$ and $1.05\times10^{-4}$ molecules$\times\mu$m$^{-3}\times$h$^{-1}$ in $c$.
$B \rightarrow\emptyset$: $k_{-2} = 0.001$ h$^{-1}$.
$A+2B \rightarrow 3B$: $k_3 = 1$ molecules$^{-2}\times\mu$m$^6\times$h$^{-1}$.
$B+2A \rightarrow 3A$: $k_4 = 1$ molecules$^{-2}\times\mu$m$^6\times$h$^{-1}$.
The diffusion coefficient for the different species were:
$D_A = D_B = 80$ $\mu$m$^2$h$^{-1}$. Reaction and diffusion rates constants are all 0 in compartment ext. 
\label{example2} ](example2.png)

For the second example, we consider a pattern-forming reaction-diffusion 
network that features two competitive auto-catalytic species $A$ and $B$ 
mutually converting into each other (Fig. \ref{example2} (a)). 
We first simulate the evolution of the system in 1 dimension using the 
$\tau$-leap algorithm [@gillespie_approximate_2001]. The corresponding 
spatio-temporal evolution of $A$ concentration over time is reported in 
Fig. \ref{example2} (c) as a 2D heat map. 
It can be observed how the system, starting from a homogeneous state, 
progressively builds up spatial reaction-diffusion patterns. 
We also simulate a square $100 \times 100$ mesh system and plot the 
final distribution of $A$, so that the shape of the pattern can be 
appreciated directly (Fig. \ref{example2} (d)). 

Next, we apply this model to two systems with higher complexity, where the synthesis rates 
of $A$ and $B$ vary depending on the region (Fig. \ref{example2} (b), (e), (g)). 
The first one (Fig. \ref{example2} (e)) represents pattern formation at the surface of a sphere, 
while the second one (Fig. \ref{example2} (g)) illustrates a similar phenomenon of pattern 
formation in a domain that takes the shape of an animal. Panels (f) and (h) in 
Fig. \ref{example2}  demonstrate the evolution of both systems 
(levels of the $A$ species) in time, highlighting the progressive evolution 
of the spatial patterns that form as a result of the subtle combination of 
the underlying auto-catalytic process with the inhomogeneous reaction landscape. 
Figure \ref{example2} (i) additionally illustrates how, due to the stochastic nature 
of the process, different patterns may arise from the same homogeneous initial state.  

# Source code and documentation

STReNGTHS's source code and documentation are distributed under the term of the 
MIT licence and can be found on the dedicated GitHub repository:

\smallskip
\texttt{https://github.com/ThibaultFillion/strengths}
\smallskip

The documentation includes tutorials and an API Reference. The tutorials 
demonstrate how to define reaction-diffusion systems by taking advantage of STReNGTHS's different 
features (environments, chemostats, boundary conditions, etc.) as well as how to carry 
out simulations and post-process the trajectories. The API Reference documents the exposed 
functions and classes. 

# Conclusions and perspectives

STReNGHTS is a new integrated and flexible platform for modeling and simulation of inhomogeneous 
reaction-diffusion systems, which aims to provide an extensible collection of 
stochastic and deterministic simulation engines. It has been designed to be easily 
integrated with new tools and we do hope it will continue to grow.  
Perspectives for future developments include:

- CPU-GPU massively parallel implementations of the existing simulation methods,
- implementing methods with dynamically adaptive time steps,
- implementing methods using smart joint use of stochastic and deterministic methods for faster simulations,
- developing a GUI that would facilitate the design of the reaction-diffusion system layouts. 

# References


