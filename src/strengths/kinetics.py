import numpy as np
import json
import strengths.value_processing as valproc
from strengths.units import *

"""
Module that computation functions
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
    :param reaction: reaction of interest.
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
    
    reaction_index = system.network.get_reaction_index(reaction)
    
    position_index = system.get_cell_index(position)
    
    r = system.network.reactions[reaction_index]
    
    ssto = r.ssto(system.network.species_labels())
    psto = r.psto(system.network.species_labels())
    
    volume = system.space.cell_vol.copy()
    
    if not r.environments==None and not system.network.environments[system.cell_env[position_index]] in r.environments :
        rates = [ unitValue(0, Units(sys=units_system, dim=r.kf_units_dimensions)),
                  unitValue(0, Units(sys=units_system, dim=r.kr_units_dimensions)) ]
    else :
        rf, rr = r.kf * volume, r.kr * volume
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
    Those are computed from the diffusion coefficient and cell size according to the method presented by David Bernstein (Bernstein, 2015) [#Bernstein2015]_.
    returns UnitValue(0, "molecule/s").convert(units_system), UnitValue(0, "molecule/s").convert(units_system)
    if positions are out of bounds or are not neighbors
         
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
    
    if (not system.space.is_within_bounds(src_position) or 
        not system.space.is_within_bounds(dst_position) or
        not system.space.are_neighbors(src_position, dst_position)) :
           return UnitValue(0, "molecule/s").convert(units_system), UnitValue(0, "molecule/s").convert(units_system)
    
    src_position_index = system.get_cell_index(src_position)
    dst_position_index = system.get_cell_index(dst_position)
    
    s = system.network.species[species_index]
    
    volume = system.space.cell_vol.copy()

    Di, Dj = 0, 0
    if isdict(s.D):
        Di = s.D.get(system.network.environments[system.space.cell_env[src_position_index]], UnitValue(0, "µm2/s"))
        Dj = s.D.get(system.network.environments[system.space.cell_env[dst_position_index]], UnitValue(0, "µm2/s"))
    else :
        Di = s.D
        Dj = s.D

    h = system.space.cell_vol**(1/3)
    k = UnitValue("0 s-1")
    if Di!=0 and Dj!=0 : 
        k = 2/(h**2 * (1/Di + 1/Dj))
    
    src_state_index = system.get_state_index(species=species_index, position=src_position_index)
    dst_state_index = system.get_state_index(species=species_index, position=dst_position_index)
    
    return (k*state.get_at(src_state_index)).convert(units_system), (k*state.get_at(dst_state_index)).convert(units_system)


def compute_dspeciesdt(system, 
                       species, 
                       position=0,
                       state=None, 
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
    
    s = system.network.species[system.network.get_species_index(species)].label
    d = 0
    for r in system.network.reactions :
        r_rates = compute_reaction_rates(system, r, position, state, units_system)
        d += (r_rates[0] - r_rates[1]) * (r.get_product_stoechiometry(s)-r.get_substrate_stoechiometry(s))
        
    p = system.space.get_cell_coordinates(system.space.get_cell_index(position))
    for c in [[p[0]+1,p[1],p[2]],
              [p[0]-1,p[1],p[2]],              
              [p[0],p[1]+1,p[2]],
              [p[0],p[1]-1,p[2]],
              [p[0],p[1],p[2]+1],
              [p[0],p[1],p[2]-1]] :
        
        if system.space._boundary_conditions["x"] == "periodical" :
            c[0] = (system.space.w+c[0])%system.space.w
        if system.space._boundary_conditions["y"] == "periodical" :
            c[1] = (system.space.h+c[1])%system.space.h
        if system.space._boundary_conditions["z"] == "periodical" :
            c[2] = (system.space.d+c[2])%system.space.d
        
        d_rates = compute_diffusion_rates(system, species, p, c, state, units_system)
        d += (d_rates[1] - d_rates[0])

    return d.convert(units_system)
    

def compute_dstatedt(system, 
                     state=None, 
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
                a.append(compute_dspeciesdt(system, s, i, state, units_system))
    return UnitArray(a, a[0].units)
