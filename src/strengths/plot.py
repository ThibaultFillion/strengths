import matplotlib.pyplot as plt
import matplotlib.colors as pltcol
from strengths.typechecking import *
from strengths.units import *
from strengths.rdspace import RDGridSpace, RDGraphSpace
import math

"""
plotting utility
"""

def plot_trajectory (output, species, position=None, *, ax=None) :
    """
    Plots the time trajectory of one or more species. 
    It relies on the Matplotlib package.
    
    :param output: simulation output
    :type output: RDOutput
    :param species: labels of the species for which the trajectories should be plotted.
        if given a label, only one trajectory is plotted.
        if given an array, trajectories of all the labels are plotted.
    :param position: position of the cell from which we want the trajectory. 
        if None, the global trajectory for the whole system is plotted instead.
    :type position: None, number, tuple or Coord like
    :param ax: (optional keyword argument) Matplotlib axis to be used (default = None).
    :type species: str or array of str.
    """

    show = False
    if ax is None :
        ax = plt.gca()
        show = True
    
    merge = isnone(position)
    if merge : 
        position=0
    
    if type(species) == str :
        ax.set_title(species + " trajectory")
        ax.set_ylabel("species quantity (" + str(output.data.units) + ")")
        ax.set_xlabel("time (" + str(output.t.units) + ")")
        ax.plot(output.t.value, output.get_trajectory(species, merge=merge, position=position).value)
    else :
        ax.set_title("species trajectories")
        ax.set_ylabel("species quantity (" + str(output.data.units) + ")")
        ax.set_xlabel("time (" + str(output.t.units) + ")")
        for s in species: 
            ax.plot(output.t.value, output.get_trajectory(s, merge=merge, position=position).value, label=s)
        ax.legend(loc="best")
    if show :
        plt.show()
    
def plot_sample_state_2D(output, species, sample, axis="auto", axis_position=0, 
                         environments=None, xmin=None, xmax=None, units_system=None,
                         *, ax=None) : 
    """
    Plots the distribution for the quantity of a species on a plan of the system space.
    It relies on the Matplotlib package.
    
    :param output: simulation output
    :type output: RDOutput   
    :param spacies: label of the species of which the state should be plotted
    :type species: str
    :param sample: index of the sample at chich the state should be taken
    :type sample: int
    :param axis: axis along which the state slice should be taken.
        accepted values are :
            
        * "auto" : the axis that seems the most appropriate is used
        * "x" : the x axis is taken 
        * "y" : the y axis used
        * "z" : the z axis is used
        
    :type axis: str
    :param axis_position: position of the slice along the slice axis
    :type axis_position: int
    :param ax: (optional keyword argument) Matplotlib axis to be used (default = None).
    """

    show = False
    if ax is None :
        ax = plt.gca()
        show = True
        
    if type(output.system.space) != RDGridSpace :
        raise ValueError("plot_sample_state_2D requires the system to have a grid space (RDGridSpace).")
    
    state_ua = 0
    
    if units_system is None :
        state_ua = output.get_state(species, sample)
    else :
        state_ua = output.get_state(species, sample).convert(units_system)
    
    state = state_ua.value.reshape(
        output.system.space.d,
        output.system.space.h,
        output.system.space.w)
    
    system = output.system
    
    env = system.space.cell_env.reshape(
        output.system.space.d,
        output.system.space.h,
        output.system.space.w).copy()
    
    if   axis == "auto" : 
        if   output.system.space.w == min(output.system.space.w, output.system.space.h, output.system.space.d) : axis = "x"
        elif output.system.space.h == min(output.system.space.w, output.system.space.h, output.system.space.d) : axis = "y"
        else :                                                                                                   axis = "z"
    
    sample_time = output.t.get_at(sample)
    if units_system is not None :
        sample_time = sample_time.convert(units_system)
        
    ax.set_title(species + "\n" + axis + " = "+str(axis_position)+" cell"+"\nt = "+str(sample_time))
    
    if   axis == "x" : 
        ax.set_xlabel("y (cell)")
        ax.set_ylabel("z (cell)")
        state = state[:,:,axis_position]
        env   = env  [:,:,axis_position]

    elif axis == "y" : 
        ax.set_xlabel("x (cell)")
        ax.set_ylabel("z (cell)")
        state = state[:,axis_position,:]
        env   = env  [:,axis_position,:]
    elif axis == "z" : 
        ax.set_xlabel("x (cell)")
        ax.set_ylabel("y (cell)")
        state = state[axis_position,:,:]
        env   = env  [axis_position,:,:]
    else : 
        raise ValueError(str(axis) + "is not an accepted axis value. accepted value are \"x\", \"y\", \"z\" and \"auto\"" )
    
    if environments is None :
        pass
    elif isarray(environments) :
        env_indices = [system.network.get_environment_index(e) for e in environments]
        for i in range(len(state)) :
            for j in range(len(state[i])) :
                if env[i, j] not in env_indices : 
                    state[i, j] = math.nan
    
    im = ax.imshow(state, vmin=xmin, vmax=xmax)
    plt.colorbar(im, label="quantity ("+str(state_ua.units)+")")
    
    if show :
        plt.show()

def plot_state_2D(system, species, axis="auto", axis_position=0, *, ax=None) : 
    """
    Plots the distribution for the quantity of a species on a plan of the system space.
    It relies on the Matplotlib package.
    """

    show = False
    if ax is None :
        ax = plt.gca()
        show = True
        
    if type(system.space) != RDGridSpace :
        raise ValueError("plot_sample_state_2D requires the system to have a grid space (RDGridSpace).")

    species_index = system.network.get_species_index(species)

    state = system.state.value.reshape(
        system.network.nspecies(),
        system.space.d,
        system.space.h,
        system.space.w)[species_index]
    
    if   axis == "auto" : 
        if   system.space.w == min(system.space.w, system.space.h, system.space.d) : axis = "x"
        elif system.space.h == min(system.space.w, system.space.h, system.space.d) : axis = "y"
        else :                                                                                                   axis = "z"
    
    if   axis == "x" : 
        ax.set_title(species+"\nx = "+str(axis_position)+" cell")
        ax.set_xlabel("y (cell)")
        ax.set_ylabel("z (cell)")
        im = ax.imshow(state[:,:,axis_position])
        plt.colorbar(im, label="quantity "+str(system.state.units))
    elif axis == "y" : 
        ax.set_title(species+"\ny = "+str(axis_position)+" cell")
        ax.set_xlabel("x (cell)")
        ax.set_ylabel("z (cell)")
        im = ax.imshow(state[:,axis_position,:])
        plt.colorbar(im, label="quantity "+str(system.state.units))
    elif axis == "z" : 
        ax.set_title(species+"\nz = "+str(axis_position)+" cell")
        ax.set_xlabel("x (cell)")
        ax.set_ylabel("y (cell)")
        im = ax.imshow(state[axis_position,:,:])
        plt.colorbar(im, label="quantity ("+str(system.state.units)+")")
    else : 
        raise ValueError(str(axis) + "is not an axxepted axis value. accepted value are \"x\", \"y\", \"z\" and \"auto\"" )
    
    if show : 
        plt.show()        
        
def plot_environments_2D(system, axis="auto", axis_position=0, env_color_dict=None, *, ax=None) : 
    """
    Plots the distribution of the reaction diffusion environments on a plan of the system space.
    env_color_dict allow to define which color should be used to represent each environment.
    If it is ignored, a gray scale will be used.
    It relies on the Matplotlib package.
    """

    show = False
    if ax is None :
        ax = plt.gca()
        show = True
        
    if type(system.space) != RDGridSpace :
        raise ValueError("plot_sample_state_2D requires the system to have a grid space (RDGridSpace).")
        
    envmap = system.space.cell_env.reshape(
        system.space.d,
        system.space.h,
        system.space.w)
    
    colormap_colors = None
    n_env = system.network.nenvironments()
    
    if isdict(env_color_dict) : 
        colormap_colors = [env_color_dict[system.network.environments[i]] for i in range(n_env)]
    elif isarray(env_color_dict) : 
        colormap_colors = env_color_dict
    else:
        colormap_colors = [(i/n_env, i/n_env, i/n_env) for i in range(n_env)]
    
    colormap_ticks =        [i+0.5                          for i in range(n_env)]
    colormap_ticks_labels = [system.network.environments[i] for i in range(n_env)]
    colormap = pltcol.ListedColormap(colormap_colors)
    
    if   axis == "auto" : 
        if   system.space.w == min(system.space.w, system.space.h, system.space.d) : axis = "x"
        elif system.space.h == min(system.space.w, system.space.h, system.space.d) : axis = "y"
        else :                                                                       axis = "z"
    
    if   axis == "x" : 
        ax.set_title("environment map\nx = "+str(axis_position)+" cell")
        ax.set_xlabel("y (cell)")
        ax.set_ylabel("z (cell)")
        im = ax.imshow(envmap[:,:,axis_position], cmap=colormap, vmax=n_env)
        plt.colorbar(im, ticks=colormap_ticks).set_ticklabels(colormap_ticks_labels)
    elif axis == "y" : 
        ax.set_title("environment map\ny = "+str(axis_position)+" cell")
        ax.set_xlabel("x (cell)")
        ax.set_ylabel("z (cell)")
        im = ax.imshow(envmap[:,axis_position,:], cmap=colormap, vmax=n_env)
        plt.colorbar(im, ticks=colormap_ticks).set_ticklabels(colormap_ticks_labels)
    elif axis == "z" : 
        ax.set_title("environment map\nz = "+str(axis_position)+" cell")
        ax.set_xlabel("x (cell)")
        ax.set_ylabel("y (cell)")
        im = ax.imshow(envmap[axis_position,:,:], cmap=colormap, vmax=n_env)
        plt.colorbar(im, ticks=colormap_ticks).set_ticklabels(colormap_ticks_labels)
    else : 
        raise ValueError(str(axis) + "is not an accepted axis value. accepted value are \"x\", \"y\", \"z\" and \"auto\"" )

    if show : 
        plt.show()

def plot_chemostats_2D(system, species, axis="auto", axis_position=0, color_no="white", color_yes="black", *, ax=None) : 
    """
    Plots the distribution of the chemostats for a given species on a plan of the system space.
    color_no and color_yes are colors to be used to indicate the absence or prence of chemostats.
    defaults are "white" and "black".
    It relies on the Matplotlib package.
    """

    show = False
    if ax is None :
        ax = plt.gca()
        show = True
        
    if type(system.space) != RDGridSpace :
        raise ValueError("plot_sample_state_2D requires the system to have a grid space (RDGridSpace).")
        
    chstts = system.chemostats.reshape(
        system.network.nspecies(),
        system.space.d,
        system.space.h,
        system.space.w)

    n = system.network.get_species_index(species)
    
    colormap_colors = [color_no, color_yes]
    colormap_ticks =        [0.5, 1.5]
    colormap_ticks_labels = ["no", "yes"]
    colormap = pltcol.ListedColormap(colormap_colors)
    
    if   axis == "auto" : 
        if   system.space.w == min(system.space.w, system.space.h, system.space.d) : axis = "x"
        elif system.space.h == min(system.space.w, system.space.h, system.space.d) : axis = "y"
        else :                                                                       axis = "z"
    
    if   axis == "x" : 
        ax.set_title("chemostat map\nx = "+str(axis_position)+" cell")
        ax.set_xlabel("y (cell)")
        ax.set_ylabel("z (cell)")
        im = ax.imshow(chstts[n,:,:,axis_position], cmap=colormap, vmax=2)
        plt.colorbar(im, ticks=colormap_ticks, label="chemostated").set_ticklabels(colormap_ticks_labels)
    elif axis == "y" : 
        ax.set_title("chemostat map\ny = "+str(axis_position)+" cell")
        ax.set_xlabel("x (cell)")
        ax.set_ylabel("z (cell)")
        im = ax.imshow(chstts[n,:,axis_position,:], cmap=colormap, vmax=2)
        plt.colorbar(im, ticks=colormap_ticks, label="chemostated").set_ticklabels(colormap_ticks_labels)
    elif axis == "z" : 
        ax.set_title("chemostat map\nz = "+str(axis_position)+" cell")
        ax.set_xlabel("x (cell)")
        ax.set_ylabel("y (cell)")
        im = ax.imshow(chstts[n,axis_position,:,:], cmap=colormap, vmax=2)
        plt.colorbar(im, ticks=colormap_ticks, label="chemostated").set_ticklabels(colormap_ticks_labels)
    else : 
        raise ValueError(str(axis) + "is not an accepted axis value. accepted value are \"x\", \"y\", \"z\" and \"auto\"" )

    if show : 
        plt.show()
