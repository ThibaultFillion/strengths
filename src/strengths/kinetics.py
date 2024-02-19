import numpy as np
import json
import strengths.value_processing as valproc
from strengths.units import *
from strengths.rdnetwork import Reaction
from strengths import value_processing as valproc
from strengths import RDGraphSpace, RDGridSpace

"""
Module that implements computation functions
associated with reaction diffusion systems.
"""

def compute_reaction_rates(system, 
                           reaction, 
                           position=0, 
                           state=None, 
                           units_system=UnitsSystem()) :
    """
    Computes and returns the reaction forward and backward rates for a given reaction in a given position 
    of a given system at a given state.
    
    :param system: system in which the reaction occurs.
    :type system: RDSystem
    :param reaction: reaction of interest. If a Reaction object is given,
        if will be used directly. if an index or a lable is given, the reaction will be
        retrived from the system's network.
    :type reaction: Reaction, number or string
    :param position: position of the cell in which the rates should be computed
    :type position: number, tuple or coord like
    :param state: state of the system. If set to None, system.state is used instead.
    :type state: UnitArray or quantity units dimensions or None
    :param units_system: units_system in which the returned rates should be expressed.
    :type units_system: UnitsSystem
    :return: [forward rate, reverse rate]
    :rtype: list
    """
    
    if isnone(state) :
        state = system.state
    
    if type(reaction) != Reaction :
        reaction_index = system.network.get_reaction_index(reaction)
        if reaction_index is None :
            raise ValueError("No index match the given reaction.")
        if isarray(reaction_index) :
            raise ValueError("Ambiguous reaction argument have multiple index. Try passing a reaction index or label instead?")
        reaction = system.network.reactions[reaction_index]

    position_index = system.get_cell_index(position)
    
    
    ssto = reaction.ssto(system.network.species_labels())
    psto = reaction.psto(system.network.species_labels())
    
    volume = system.space.get_cell_vol_array().get_at(position_index)
    
    if reaction.environments is not None :
        environment_index = system.space.get_cell_env_array()[position_index]
        environment_label = system.network.environments[environment_index]
        if not environment_label in reaction.environments :
            return [ UnitValue(0, Units(sys=units_system, dim=UnitsDimensions(quantity=1, space=0, time=-1))),
                     UnitValue(0, Units(sys=units_system, dim=UnitsDimensions(quantity=1, space=0, time=-1))) ]

    rf = reaction.kf * volume
    rr = reaction.kr * volume
    for i in range(system.network.nspecies()) :
        state_index = system.get_state_index(species=i, position=position_index)
        rf *= (state.get_at(state_index)/volume)**ssto[i]
        rr *= (state.get_at(state_index)/volume)**psto[i]
    return rf.convert(units_system), rr.convert(units_system)


def compute_diffusion_rates(system, 
                            species, 
                            src_position,
                            dst_position,
                            state=None, 
                            units_system=UnitsSystem()) :
    """
    Computes and returns the forward and backward diffusion rates for a 
    given species of a given system at a given state from a given source cell to another destination cell.
    Those are computed from the diffusion coefficient and cell size according to the method presented by David Bernstein (Bernstein, 2015) [#Bernstein2015]_,
    with some modification in the case of graph diffusion (with RDGraphSpace).
         
    :param system: system in which the reaction occurs.
    :type system: RDSystem
    :param species: species of interest.
    :type species: Species, number or string
    :param src_position: position of the diffusion source cell.
    :type src_position: number, tuple or coord like
    :param dst_position: position of the diffusion destination cell.
    :type dst_position: number, tuple or coord like
    :param state: state of the system. If set to None, system.state is used instead.
    :type state: UnitArray or quantity units dimensions or None
    :param units_system: units_system in which the returned rates should be expressed.
    :type units_system: UnitsSystem
    :return: [forward rate,reverse rate]
    :rtype: list
    
    references
    ^^^^^^^^^^
    
    .. [#Bernstein2015] Bernstein, D. (2005). Simulating mesoscopic reaction-diffusion systems using the Gillespie algorithm. Physical Review E, 71(4), Article 041103. https://doi.org/10.1103/PhysRevE.71.041103
    """
    
    if isnone(state) :
        state = system.state
    
    species_index = system.network.get_species_index(species)
        
    src_position_index = system.get_cell_index(src_position)
    dst_position_index = system.get_cell_index(dst_position)
    
    #checking neighbors :
    if type(system.space) == RDGridSpace and not system.space.are_neighbors(src_position_index, dst_position_index) :
           raise ValueError("Cannot compute diffusion rates for non-neighbor cells.")
    if type(system.space) == RDGraphSpace and system.space.get_edge(src_position_index, dst_position_index) is None :
           raise ValueError("Cannot compute diffusion rates for non-neighbor cells.")
           
    species = system.network.species[species_index]
     
    cellenv = system.space.get_cell_env_array()
    envid_i = cellenv[src_position_index]
    envid_j = cellenv[dst_position_index]
    envl_i = system.network.environments[envid_i]
    envl_j = system.network.environments[envid_j]
    
    Di, Dj = 0, 0
    if isdict(species.D):
        Di = valproc.get_value_in_env(value=species.D, environment=envl_i, default=UnitValue(0, "µm2/s"))
        Dj = valproc.get_value_in_env(value=species.D, environment=envl_j, default=UnitValue(0, "µm2/s"))
    else :
        Di = species.D
        Dj = species.D

    if type(system.space) == RDGraphSpace:
        volumes = system.space.get_cell_vol_array()
        Vi = volumes.get_at(src_position_index)
        Vj = volumes.get_at(dst_position_index)
        hi = Vi**(1/3)
        hj = Vj**(1/3)
        surface  = system.space.get_edge(src_position_index, dst_position_index).surface
        distance = system.space.get_edge(src_position_index, dst_position_index).distance 
        
        Dij = UnitValue("0 µm2/s")
        
        if Di.value!=0 and Dj.value!=0:
            Dij=(hi+hj)/(hi/Di + hj/Dj)
        
        kf = Dij * surface / (Vi * distance)
        kr = Dij * surface / (Vj * distance)
                
        src_state_index = system.get_state_index(species=species_index, position=src_position_index)
        dst_state_index = system.get_state_index(species=species_index, position=dst_position_index)
        
        return (kf*state.get_at(src_state_index)).convert(units_system), (kr*state.get_at(dst_state_index)).convert(units_system)

    elif type(system.space) == RDGridSpace:
        h = system.space.cell_vol**(1/3)
        k = UnitValue("0 s-1")
        if Di!=0 and Dj!=0 : 
            k = 2/(h**2 * (1/Di + 1/Dj))
        
        src_state_index = system.get_state_index(species=species_index, position=src_position_index)
        dst_state_index = system.get_state_index(species=species_index, position=dst_position_index)
        
        return (k*state.get_at(src_state_index)).convert(units_system), (k*state.get_at(dst_state_index)).convert(units_system)

    else :
        raise ValueError("Invalid system space type. Must be RDGridSpace or RDGraphSpace.")

def _compute_dspeciesdt_grid(system, 
                       species, 
                       position=0,
                       state=None, 
                       apply_chemostats=True,
                       units_system=UnitsSystem()) :
    
    species_index = system.network.get_species_index(species)
    species_label = system.network.species[species_index].label
    
    d = 0
    for reaction in system.network.reactions :
        rates = compute_reaction_rates(system, reaction, position, state, units_system)
        d += (rates[0] - rates[1]) * (reaction.get_product_stoechiometry(species_label)-reaction.get_substrate_stoechiometry(species_label))
        
    p = system.space.get_cell_coordinates(system.space.get_cell_index(position))
    for c in [[p[0]+1,p[1],p[2]],
              [p[0]-1,p[1],p[2]],              
              [p[0],p[1]+1,p[2]],
              [p[0],p[1]-1,p[2]],
              [p[0],p[1],p[2]+1],
              [p[0],p[1],p[2]-1]] :
        
        if system.space._boundary_conditions["x"] == "periodical" and system.space.w > 1:
            c[0] = (system.space.w+c[0])%system.space.w
        if system.space._boundary_conditions["y"] == "periodical" and system.space.h > 1:
            c[1] = (system.space.h+c[1])%system.space.h
        if system.space._boundary_conditions["z"] == "periodical" and system.space.d > 1:
            c[2] = (system.space.d+c[2])%system.space.d
        
        if system.space.is_within_bounds(c) :
            d_rates = compute_diffusion_rates(system, species, p, c, state, units_system)
            d += (d_rates[1] - d_rates[0])

    if apply_chemostats and system.chemostats[system.space.get_cell_index(position)]:
    	return UnitValue(0, "molecule/s").convert(units_system)
    
    return d.convert(units_system)

def _compute_dspeciesdt_graph(system, 
                       species, 
                       position=0,
                       state=None, 
                       apply_chemostats=True,
                       units_system=UnitsSystem()) :
    
    species_index = system.network.get_species_index(species)
    species_label = system.network.species[species_index].label
    position = system.space.get_cell_index(position)
    
    d = 0
    for reaction in system.network.reactions :
        rates = compute_reaction_rates(system, reaction, position, state, units_system)
        d += (rates[0] - rates[1]) * (reaction.get_product_stoechiometry(species_label)-reaction.get_substrate_stoechiometry(species_label))
            
    neighbor_indices = []
    for j in range(system.space.size()):
        if j != position :
            if system.space.get_edge(position, j) is not None :
                neighbor_indices.append(j)
    
    for j in neighbor_indices :
        d_rates = compute_diffusion_rates(system, species, position, j, state, units_system)
        d += (d_rates[1] - d_rates[0])
    
    if apply_chemostats and system.chemostats[system.space.get_cell_index(position)]:
    	return UnitValue(0, "molecule/s").convert(units_system)
    
    return d.convert(units_system)

def compute_dspeciesdt(system, 
                       species, 
                       position=0,
                       state=None, 
                       apply_chemostats=True,
                       units_system=UnitsSystem()) :
    """
    Computes and returns time derivative of the system state for a given 
    species at a given position.
    
    :param system: system in which the reaction occurs.
    :type system: RDSystem
    :param species: species of interest.
    :type species: Species, number or string
    :param position: position of the cell in which the rates should be computed
    :type position: number, tuple or coord like
    :param state: state of the system. If set to None, system.state is used instead.
    :type state: UnitArray or quantity units dimensions or None
    :param units_system: units_system in which the returned rates should be expressed.
    :type units_system: UnitsSystem
    :return: time derivative of the species quantity at the given position of the system.
    :rtype: UnitValue with quantity per time unit dimensions
    """
    if   type(system.space) == RDGridSpace :
        return _compute_dspeciesdt_grid(
            system = system, 
            species = species, 
            position = position,
            state = state, 
            apply_chemostats = apply_chemostats,
            units_system = units_system)

    elif type(system.space) == RDGraphSpace :
        return _compute_dspeciesdt_graph(
            system = system, 
            species = species, 
            position = position,
            state = state, 
            apply_chemostats = apply_chemostats,
            units_system = units_system)

    else :
        raise ValueError("System space type must be of type RDGridSpace or RDGraphSpace.")

def compute_dstatedt(system, 
                     state=None, 
                     apply_chemostats = True,
                     units_system=UnitsSystem()) :
    """
    Computes and returns time derivative of the system state for all species at all poition.
    The returned UnitArray have the same organisation than state.
    
    :param system: system in which the reaction occurs.
    :type system: RDSystem
    :param state: state of the system. If set to None, system.state is used instead.
    :type state: UnitArray or quantity units dimensions or None
    :param units_system: units_system in which the returned rates should be expressed.
    :type units_system: UnitsSystem
    :return: [dxdt for x in state]
    :rtype: UnitArray with quantity per time unit dimensions
    """
    
    a = []
    for s in range(system.network.nspecies()) :
        for i in range(system.space.size()) :
                a.append(compute_dspeciesdt(system, s, i, state, apply_chemostats, units_system))
    return UnitArray(a, a[0].units)
