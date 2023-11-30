from strengths.units import *
from strengths.rdsystem import RDSystem
from strengths.rdengine import RDEngineBase
from strengths.rdoutput import  RDTrajectory
from strengths.rdspace import RDGridSpace, RDGraphSpace
from strengths.typechecking import *

import numpy as np
import ctypes

def build_reaction_rate_constant_array(reactions, units_system) :
    """
    return the array of the forward reaction rate constants numerical values.

    :param reactions: array of irreversible forward reactions.
    :type reactions: array of Reaction
    :param units_system: units_system in which the rates should be expressed.
    :type units_system: UnitsSystem
    :rtype: array of int
    """

    return [r.kf.convert(units_system).value for r in reactions]

def build_reaction_environment_boolean_matrix(reactions, environments) :
    """
    return reaction environment boolean matrix in linear form.
    let *m* be such a *n_reactions*x*n_environments* matrix. *r* reffer to a reaction (row), and *e* refer to an environment (column).
    *m_{r,e} = 1*, if reaction *j* can happend in environment *i*, *0* otherwise.
    the linear form is *[m_{0,0},m_{0,1},... m_{0,n_environments},m_{1,0}, ..., m_{n_reactions,n_environments}]*.

    :param reactions: array of irreversible forward reactions.
    :type reactions: array of Reaction
    :param environments: array of environment labels.
    :type environments: array of str
    :rtype: array of int
    """

    nr = len(reactions)
    ne = len(environments)
    r_env = np.zeros(nr*ne, dtype=int)
    for r in range(nr) :
        for e in range(ne) :
            if reactions[r].environments == None :
                r_env[r * ne + e] = 1
            else :
                r_env[r * ne + e] = int(environments[e] in reactions[r].environments)
    return r_env

def build_substrate_stoechiometric_matrix(species, reactions) :
    """
    return substrate stoechiometric matrix in linear form.
    let *m* be such a *n_species*x*n_reactions* matrix. *s* reffer to a species (row), and *r* to a reaction (column).
    *m_{s,r}* is the stoechiometric coefficent associated with *s* in the substrate side of reaction *r*.
    the linear form is *[m_{0,0},m_{0,1},... m_{0,n_reactions},m_{1,0}, ..., m_{n_species,n_reactions}]*.

    :param species: array of species.
    :type species: array of Species
    :param reactions: array of irreversible forward reactions.
    :type reactions: array of Reaction
    :rtype: array of int
    """

    n_species   = len(species)
    n_reactions = len(reactions)
    species_labels = [s.label for s in species]
    sub = np.zeros(n_species*n_reactions, dtype=int)
    for s in range(n_species):
        for r in range(n_reactions):
            sub[s*n_reactions+r] = reactions[r].ssto(species_labels)[s]
    return sub

def build_stoechiometric_difference_matrix(species, reactions) :
    """
    return stoechiometric difference matrix in linear form.
    let *m* be such a *n_species*x*n_reactions* matrix. *s* reffer to a species (row), and *r* to a reaction (column).
    *m_{s,r}* is the difference of stoechiometric coefficent (products-substrates) for species *s* in reaction *r*.
    the linear form is *[m_{0,0},m_{0,1},... m_{0,n_reactions},m_{1,0}, ..., m_{n_species,n_reactions}]*.

    :param species: array of species.
    :type species: array of Species
    :param reactions: array of irreversible forward reactions.
    :type reactions: array of Reaction
    :rtype: array of int
    """

    n_species   = len(species)
    n_reactions = len(reactions)
    species_labels = [s.label for s in species]
    sto = np.zeros(n_species*n_reactions, dtype=int)
    for s in range(n_species):
        for r in range(n_reactions):
            sto[s*n_reactions+r] = reactions[r].dsto(species_labels)[s]
    return sto

def build_diff_coef_environment_matrix(species, environments, units_system) :
    """
    return the species diffusion coefficient matrix in linear form.
    let *m* be such a *n_species*x*n_environments* matrix. *s* reffer to a species (row), and *e* refer to an environment (column).
    *m_{s,e}* is the diffusion coefficent of species *s* in the enviroment *e*.
    the linear form is *[m_{0,0},m_{0,1},... m_{0,n_environments},m_{1,0}, ..., m_{n_species,n_environments}]*.

    :param species: array of species.
    :type species: array of Species
    :param environments: array of environment labels.
    :type environments: array of str
    :rtype: array of int
    """

    n_species=len(species)
    n_env=len(environments)
    D = np.zeros(n_species*n_env, dtype=float)
    for s in range(n_species):
        for e in range(n_env):
            D[s*n_env+e] = valproc.get_value_in_env(
                value = species[s].D, 
                environment = environments[e], 
                default = UnitValue(0, "µm2/s")).convert(units_system).value
    return D

def make_ctypes_array(a, t) :
    """
    make a ctype array of type t from the array a.
    """

    l = len(a)
    b = (t*l)()
    for i in range(l):
        b[i] = a[i]
    return b

class LibRDEngine(RDEngineBase) :
    """
    Engine internally relying on compiled dynamic libraries / DLLs.
    """

    def __init__(self, lib, option="", description = "", requires_molecules=False) :
        self._lib = lib
        self._requires_molecules = requires_molecules
        self._simulation_unfinished = 1
        
        lib.GetProgress.restype = ctypes.c_double
        lib.GetT.restype = ctypes.c_double
        super(LibRDEngine, self).__init__(option, description)

    def setup(self, script) :
        
        self._script = script.copy()
        
        units_system = script.units_system.copy()

        if self._requires_molecules : 
            units_system.quantity = "molecule"
            
        self._units_system = units_system
            
            
        species = script.system.network.species

        reactions = []

        for r in script.system.network.reactions :
            rf, rr =r.split()
            reactions.append(rf)
            reactions.append(rr)

        environments = script.system.network.environments
        
        if type(script.system.space) == RDGridSpace :
            self._setup_grid(script, units_system, species, reactions, environments)
        elif type(script.system.space) == RDGraphSpace :
            self._setup_graph(script, units_system, species, reactions, environments)
        else :
            raise TypeError("unsupported space type.")
        
    def _setup_graph(self, script, units_system, species, reactions, environments) :

        res = self._lib.InitializeGraph(
            #n_nodes
                ctypes.c_int(len(script.system.space.nodes)),

            #n_species
                ctypes.c_int(len(species)),
                
            #n_reactions
                ctypes.c_int(len(reactions)),
                
            #n_env
                ctypes.c_int(len(environments)),

            #n_edges
                ctypes.c_int(len(script.system.space.edges)),

            #edge_i
                make_ctypes_array([edge.i for edge in script.system.space.edges], ctypes.c_int),
            
            #edge_j
                make_ctypes_array([edge.j for edge in script.system.space.edges], ctypes.c_int), 
                
            #edge_sfc
                make_ctypes_array(
                    UnitArray(
                        [edge.surface for edge in script.system.space.edges], 
                        Units(sys=units_system, dim=surface_units_dimensions())
                        ).value,
                    ctypes.c_double),
            
            #edge_dst
                make_ctypes_array(
                    UnitArray(
                        [edge.distance for edge in script.system.space.edges], 
                        Units(sys=units_system, dim=space_units_dimensions())
                        ).value,
                    ctypes.c_double),

            #cell_state
                make_ctypes_array(script.system.state.convert(units_system).value, ctypes.c_double),
                
            #cell_chstt
                make_ctypes_array(script.system.chemostats,                         ctypes.c_int),
                
            #cell_env
                make_ctypes_array(script.system.space.get_cell_env_array(),                    ctypes.c_int),
                
            #cell_vol
                make_ctypes_array(script.system.space.get_cell_vol_array().convert(units_system).value, ctypes.c_double),
                
            #k
                make_ctypes_array(build_reaction_rate_constant_array(reactions, units_system),             ctypes.c_double),
                
            #sub
                make_ctypes_array(build_substrate_stoechiometric_matrix(species, reactions),               ctypes.c_int),
                
            #sto
                make_ctypes_array(build_stoechiometric_difference_matrix(species, reactions),              ctypes.c_int),
                
            #r_env
                make_ctypes_array(build_reaction_environment_boolean_matrix(reactions, environments),      ctypes.c_int),
                
            #D
                make_ctypes_array(build_diff_coef_environment_matrix(species, environments, units_system), ctypes.c_double),
                                
            #n_sample
                ctypes.c_int(len(script.t_sample)),
                
            #t_sample
                make_ctypes_array(script.t_sample.convert(units_system).value,   ctypes.c_double),
                
            #sampling_policy
                ctypes.c_char_p(script.sampling_policy.encode()),                

            #sampling_interval
                ctypes.c_double(script.sampling_interval.convert(units_system).value),

            #t_max
                ctypes.c_double(script.t_max.convert(units_system).value),                

            #time_step
                ctypes.c_double(script.time_step.convert(units_system).value),
                
            #seed
                ctypes.c_int(script.rng_seed),
                
            #option
                ctypes.c_char_p(self.option.encode())
                )

        if   res == 1 :
            raise Exception("Invalid option argument : \""+engine.get_option()+"\".")
        elif res == 2 :
            raise Exception("Invalid boundary conditions.")
            
    def _setup_grid(self, script, units_system, species, reactions, environments) :
        
        res = self._lib.Initialize3D(
            #w
                ctypes.c_int(script.system.space.w),
                
            #h
                ctypes.c_int(script.system.space.h),
                
            #d
                ctypes.c_int(script.system.space.d),
                
            #n_species
                ctypes.c_int(len(species)),
                
            #n_reactions
                ctypes.c_int(len(reactions)),
                
            #n_env
                ctypes.c_int(len(environments)),
                
            #cell_state
                make_ctypes_array(script.system.state.convert(units_system).value, ctypes.c_double),
                
            #cell_chstt
                make_ctypes_array(script.system.chemostats,                         ctypes.c_int),
                
            #cell_env
                make_ctypes_array(script.system.space.cell_env,                    ctypes.c_int),
                
            #cell_vol
                ctypes.c_double(script.system.space.cell_vol.convert(units_system).value),
                
            #k
                make_ctypes_array(build_reaction_rate_constant_array(reactions, units_system),             ctypes.c_double),
                
            #sub
                make_ctypes_array(build_substrate_stoechiometric_matrix(species, reactions),               ctypes.c_int),
                
            #sto
                make_ctypes_array(build_stoechiometric_difference_matrix(species, reactions),              ctypes.c_int),
                
            #r_env
                make_ctypes_array(build_reaction_environment_boolean_matrix(reactions, environments),      ctypes.c_int),
                
            #D
                make_ctypes_array(build_diff_coef_environment_matrix(species, environments, units_system), ctypes.c_double),
                
            #boundary_conditions_x
                ctypes.c_char_p((script.system.space.get_boundary_conditions()["x"]).encode()),
                
            #boundary_conditions_y
                ctypes.c_char_p((script.system.space.get_boundary_conditions()["y"]).encode()),
                
            #boundary_conditions_z
                ctypes.c_char_p((script.system.space.get_boundary_conditions()["z"]).encode()),                
                
            #n_sample
                ctypes.c_int(len(script.t_sample)),
                
            #t_sample
                make_ctypes_array(script.t_sample.convert(units_system).value,   ctypes.c_double),
                
            #sampling_policy
                ctypes.c_char_p(script.sampling_policy.encode()),                

            #sampling_interval
                ctypes.c_double(script.sampling_interval.convert(units_system).value),

            #t_max
                ctypes.c_double(script.t_max.convert(units_system).value),                
                
            #time_step
                ctypes.c_double(script.time_step.convert(units_system).value),
                
            #seed
                ctypes.c_int(script.rng_seed),
                
            #option
                ctypes.c_char_p(self.option.encode())
                )

        if   res == 1 :
            raise Exception("Invalid option argument : \""+engine.get_option()+"\".")
        elif res == 2 :
            raise Exception("Invalid boundary conditions.")
            
    def run(self, breathe_dt) :
        
        self._simulation_unfinished = self._lib.Run(breathe_dt)
        return bool(self._simulation_unfinished)

    def iterate(self) :
        
        self._simulation_unfinished = self._lib.Iterate()
        return bool(self._simulation_unfinished)


    def iterate_n(self, n_iterations) :
        
        self._simulation_unfinished = self._lib.IterateN(n_iterations)
        return bool(self._simulation_unfinished)


    def get_progress(self) :
        
        progress = self._lib.GetProgress()
        return float(progress)

    def sample(self) :
        
        self._lib.Sample()
        
    def is_complete(self) : 
        
        return not bool(self._simulation_unfinished)
    
    def _count_samples(self) :
        return self._lib.GetNSamples()

    def _get_t_sample(self) :
        n_sample = self._count_samples()
        t_sample_ = (n_sample*ctypes.c_double)()
        self._lib.GetTSample(t_sample_)
        t_sample = np.zeros(n_sample)
        for i in range(n_sample):
            t_sample[i] = t_sample_[i]
            
        return UnitArray(value=t_sample, 
                         units=Units(
                             sys=self._units_system ,
                             dim=time_units_dimensions()),
                         check_value=False
                         ).convert(self._script.units_system)

    def _get_data(self) :
        n_sample = self._count_samples()
        data_len = n_sample*self._script.system.state_size()
        
        data_ = (data_len*ctypes.c_double)()
        self._lib.GetOutput(data_)
        data = np.zeros(data_len)
        for i in range(data_len):
            data[i] = data_[i]
            
        return UnitArray(value=data, 
                         units=Units(
                             sys=self._units_system ,
                             dim=quantity_units_dimensions()),
                         check_value=False
                         ).convert(self._script.units_system)
    
    def get_output(self) :
        return RDTrajectory(
            data = self._get_data(), 
            t_sample = self._get_t_sample(), 
            system = self._script.system,
            script = self._script,
            engine_description = self.description, 
            engine_option = self.option
            )
        
    def finalize(self) :
        
        self._lib.Finalize()
