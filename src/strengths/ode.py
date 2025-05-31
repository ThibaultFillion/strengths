import copy
import numpy as np
from scipy.integrate import LSODA
from strengths.rdspace import RDGridSpace, RDGraphSpace
from strengths.value_processing import get_value_in_env
from strengths.coarsegrain import grid_to_graph
from strengths.units import *

def _get_diffusion_edges(system):
    graph = None
    if type(system.space)==RDGridSpace:
        graph = grid_to_graph(system.space)
    else:
        graph = system.space
    return copy.deepcopy(graph.edges)

def _compute_kd(system, s, i, j, units_system):
    """
    Compute diffusion rate constants from 
    node i to node j according to reference 1 
    for grid spaces or with an extended method 
    based on this same reference for graph spaces.
    References:
    [1] Bernstein, D. (2005). Simulating mesoscopic 
    reaction-diffusion systems using the Gillespie algorithm. 
    Physical Review E, 71(4), Article 041103. 
    https://doi.org/10.1103/PhysRevE.71.041103
    """
    vi = system.space.get_cell_vol(i).convert(units_system).value
    vj = system.space.get_cell_vol(j).convert(units_system).value
    ei = system.network.environments[system.space.get_cell_env(i)]
    ej = system.network.environments[system.space.get_cell_env(j)]
    Di = get_value_in_env(system.network.species[s].D, ei, UnitValue("0 m2/s")).convert(units_system).value
    Dj = get_value_in_env(system.network.species[s].D, ej, UnitValue("0 m2/s")).convert(units_system).value
    li = vi**(1./3.)
    lj = vj**(1./3.)
    
    if Di == 0 or Dj == 0:
        return 0.
    
    Dij = (li+lj)/((li/Di)+(lj/Dj))
    
    if type(system.space) == RDGridSpace:
        return Dij/(li**2)
    else:
        edge = system.space.get_edge(i, j)
        sij = edge.surface.convert(units_system).value
        dij = edge.distance.convert(units_system).value
        return Dij*sij/(vi*dij)
    
def build_ode_function(system, units_system):
    """
    returns an ODE function with the argument format 
    expected by SciPy's integration functions and classes
    """
    
    # diffusion edges
    edges = _get_diffusion_edges(system)
    
    # dimensions
    L = system.space.size()
    N = len(system.network.species)
    M = len(system.network.reactions)
    O = len(edges)
    
    # reaction rate constants 
    # accounting for node environment and volume
    k_fwd = np.zeros(L*M)
    k_rev = np.zeros(L*M)
    for i in range(L):
        for j in range(M):
            r = system.network.reactions[j]
            v = system.space.get_cell_vol(i).convert(units_system).value
            e = system.network.environments[system.space.get_cell_env(i)]
            k_fwd[i*M+j] = get_value_in_env(r.kf, e, UnitValue(0, Units(units_system, r.kf_units_dimensions()))).convert(units_system).value * v**(1-r.order())
            k_rev[i*M+j] = get_value_in_env(r.kr, e, UnitValue(0, Units(units_system, r.kr_units_dimensions()))).convert(units_system).value * v**(1-r.rorder())
    
    # diffusion edges    
    kd_fwd = np.zeros(O*N)
    kd_rev = np.zeros(O*N)
    for i in range(O):
        for s in range(N):
            kd_fwd[i*N+s] = _compute_kd(system, s, edges[i].i, edges[i].j, units_system)
            kd_rev[i*N+s] = _compute_kd(system, s, edges[i].j, edges[i].i, units_system)
            
    # stoichiometry
    ssto = [r.ssto(system.network.species_labels()) for r in system.network.reactions]
    psto = [r.psto(system.network.species_labels()) for r in system.network.reactions]
    dsto = [r.dsto(system.network.species_labels()) for r in system.network.reactions]
    
    # chemostats
    chstts = 1-system.chemostats
    
    def ode_function(t, x): # x is the state
        dxdt = np.zeros(L*N) # output = state time derivative
        rates = np.zeros(L*M)
        d_rates = np.zeros(O*N)
        
        # compute reaction rates
        for i in range(L):
            for j in range(M):
                rate_fwd = k_fwd[i*M+j]
                rate_rev = k_rev[i*M+j]
                for s in range(N):
                    rate_fwd *= x[s*L+i]**ssto[j][s]
                    rate_rev *= x[s*L+i]**psto[j][s]
                rates[i*M+j] = rate_fwd - rate_rev
        
        for i in range(O):
            for s in range(N):
                d_rates[i*N+s] = -kd_fwd[i*N+s]*x[s*L+edges[i].i] + kd_rev[i*N+s]*x[s*L+edges[i].j]
        
        # compute dxdt
        for i in range(L):
            for j in range(M):
                for s in range(N):
                    dxdt[s*L+i] += rates[i*M+j]*dsto[j][s]*chstts[s*L+i]

        for i in range(O):
            for s in range(N):
                dxdt[s*L+edges[i].i] += d_rates[i*N+s]*chstts[s*L+edges[i].i]
                dxdt[s*L+edges[i].j] -= d_rates[i*N+s]*chstts[s*L+edges[i].j]
        
        return dxdt
                
    return ode_function
