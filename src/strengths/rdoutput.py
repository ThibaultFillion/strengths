from strengths.units import *
from strengths.rdsystem import RDSystem, rdsystem_from_dict, rdsystem_to_dict
from strengths.rdscript import RDScript, rdscript_from_dict, rdscript_to_dict
from strengths.rdspace import RDGridSpace, RDGraphSpace
from strengths.typechecking import *

import json
import copy
import random

class RDTrajectory :
    """
    Trajectory of a reaction-diffusion system.

    :param data: trajectory data [sample index, species index, cell index]
    :type data: UnitArray with quantity units dimensions
    :param t_sample: system time associated with each sample of the trajectory
    :type t_sample: UnitArray with time units dimensions
    :param system: reaction diffusion system associated with the trajectory
    :type system: RDSystem    

    :param script: simulation script associated with the trajectory.
    :type script: RDScript ot None
    :param engine_description: description of the engine used for the simulation.
    :type engine_description: str or None
    :param engine_option: option used with the engine used for the simulation.
    :type engine_option: str or None
    :param cg_map: coarse graining index map or None. if not None, it means
        that the data results from an uncoarsegraining operation, with cgmap.
    :type cg_map: array of int or None
    """

    def __init__ (self, 
                  data, 
                  t_sample, 
                  system,
                  
                  script = None, 
                  engine_description = None, 
                  engine_option = None,
                  cgmap = None
                  ):
        """
        constructor
        """
        self._data = data.copy()
        self._t = t_sample.copy()
        self._system = system.copy()
        
        if isnone(script) :
            self._script = None
        else :            
            self._script = script.copy()
            
        self._engine_description = engine_description
        self._engine_option = engine_option
        self._cgmap = cgmap

    @property
    def t(self) :
        """
        sampling time points.
        """

        return self._t

    @property
    def data(self) :
        """
        array of the system states sucessively sampled at the time points in self.t.
        """

        return self._data

    @property
    def system(self) :
        """
        reaction diffusion system.
        """

        return self._system
    
    @property
    def script(self) :
        """
        reaction diffusion simulation script.
        """

        return self._script

    @property
    def engine_description(self) :
        """
        Description of the engine used for the simulation.
        """
        
        return self._engine_description

    @property
    def engine_option(self) :
        """
        Option used for the engine used for the simulation.
        """
        
        return self._engine_option

    @property
    def cgmap(self) :
        """
        Coarse graining index map which has been used to uncoarsegrain the trajectory.
        """
        
        return self._cgmap
    
    def ncells(self):
        """
        Returns the number of species.
        """

        return self.system.space.size()

    def nspecies(self) :
        """
        Returns the number of species.
        """

        return self.system.network.nspecies()

    def nsamples(self) :
        """
        Returns the number of samples.
        """

        return len(self.t)

    def get_trajectory(self, species, position=0, merge=False) :
        """
        Returns the trajectory of a given species. if merge=False, it is the trajectory at the given position,
        otherwise, it is the global trajectory in the whole system.

        :param species: species for which trajectory should be returned.
        :type species: Species, number or str
        :param position: position of the cell in which the trajectory should be taken.
            it is ignored if merge=True.
        :type position: number, tuple or Coord like.
        :param merge: if True, position is ignored, and the global trajectory of the species is trturned.
        :type merge: bool
        :returns: local ot global trajectory of the species.
        :rtype: UnitArray of quantity units dimensions
        """

        species_index = self.system.network.get_species_index(species)
        if isnone(species_index) :
            raise Exception("Undefined species \""+species+"\".")

        cell_index = self.system.space.get_cell_index(position)

        if not merge :
            return UnitArray(self.data.value.reshape((self.nsamples(), self.nspecies(), self.ncells()))[:,species_index, cell_index], self.data.units, check_value=False)
        else : #merge
            return UnitArray([sum(state) for state in self.data.value.reshape((self.nsamples(), self.nspecies(), self.ncells()))[:,species_index, :]], self.data.units, check_value=False)

    def get_state(self, species, sample) :
        """
        Returns the state of a given species at a given sample index.
        if species is None, the whole state for the given sample index is returned.

        :param species: species, species label or species index
        :type species: Species, int or str or None
        :param sample: sample index
        :type sample: int
        :returns: sampled species state at the given sample index.
        """
        if isnone(species) :
    
            sample_index = sample
    
            return UnitArray(
                self.data.value.reshape((self.nsamples(), self.nspecies()*self.ncells()))[sample_index, : ], self.data.units, check_value=False)
        else:
            species_index = self.system.network.get_species_index(species)
            if isnone(species_index) :
                raise Exception("Undefined species \""+species+"\".")
    
            sample_index = sample
    
            return UnitArray(
                self.data.value.reshape((self.nsamples(), self.nspecies(), self.ncells()))[sample_index, species_index, :], self.data.units, check_value=False)

    def get_trajectory_point(self, species, sample, position) :
        """
        Returns the value of a given species trajectory, at a given sample, at a given position.

        :param species: species, species label or index
        :type species: Species, int or str
        :param sample: sample index
        :type sample: int
        :param position: position in the cell space.
            can be a linear index (int) or a class with x,y,z properties/members representing spatial cooridnates in a MeghGrid.
        :returns: sample value at the given position for the given species (number).
        """

        species_index = self.system.network.get_species_index(species)
        if isnone(species_index) :
            raise Exception("Undefined species \""+species+"\".")

        cell_index = self.system.space.get_cell_index(position)

        sample_index = sample

        return self.data.get_at(sample_index*self.nspecies()*self.ncells() + species_index*self.ncells() + cell_index)
    
    def _get_sample_index_closest(self, t) :

        if len(self.t)==0:
            return None

        if t<=self.t.get_at(0) :
            return 0
        
        if t>=self.t.get_at(self.nsamples()-1) :
            return self.nsamples()-1
        
        for i in range(self.nsamples()-1) :
            if t>=self.t.get_at(i) and t<self.t.get_at(i+1):
                dt0 = t-self.t.get_at(i)
                dt1 = self.t.get_at(i+1)-t
                
                if dt0<=dt1 : 
                    return i
                else : 
                    return i+1

    def _get_sample_index_infeq(self, t) :

        if len(self.t)==0:
            return None

        if t<self.t.get_at(0) :
            return None
        
        if t>=self.t.get_at(self.nsamples()-1) :
            return self.nsamples()-1
        
        for i in range(self.nsamples()-1) :
            if t>=self.t.get_at(i) and t<self.t.get_at(i+1):
                return i

    def _get_sample_index_supeq(self, t) :

        if len(self.t)==0:
            return None

        if t<=self.t.get_at(0) :
            return 0
        
        if t==self.t.get_at(self.nsamples()-1) :
            return self.nsamples()-1

        if t>self.t.get_at(self.nsamples()-1) :
            return None
        
        for i in range(self.nsamples()-1) :
            if t>self.t.get_at(i) and t<=self.t.get_at(i+1):
                return i+1

    def get_sample_index(self, t, policy="closest") :
        """
        return the sample index i which self.t[i] is the closest from t.
        t must be a UnitValue (or its string representation) in time units (ie. "s", "min", etc.) and  policy must be a string with one of the following values :
            
        * if policy = "closest", the closest index is returned, or None if there is no sample.
        * if policy = "supeq", the closest index for which t[i]>=t is returned. None is returned if there is not such value.
        * if policy = "infeq", the closest index for which t[i]<=t is returned. None is returned if there is not such value.
        """
        t = UnitValue(t, self.t.units, convert=True)
        if policy not in ["closest", "sup", "inf", "supeq", "infeq"] :
            raise ValueError("invalid policy value. accepted values are : \"closest\", \"supeq\" and \"infeq\".")
        found=False
        
        if policy=="closest" : return self._get_sample_index_closest(t)
        if policy=="infeq"   : return self._get_sample_index_infeq(t)
        if policy=="supeq"   : return self._get_sample_index_supeq(t)
        
def save_rdtrajectory(so, path, separate_data=True) :
    """
    Saves a simulation output as a file.

    :param so: simulation output to be saved.
    :type so: RDOutput
    :param path: path of the output file to be created or replaced.
        the .json extension suffix is added if absent.
    :type path: str
    :param separate_data: specifies if the data should be saved in a separate file.
        if true, the trajectory data are saved in a different file using the numpy.save function [#numpy_save]_.
        This makes the saving and loading faster, especially for large simulation outputs.
        if filename.json is the name of the json file,
        the data are saved as filename_data.npy (NPY format [#numpy_npy]_).
    :type separate_data: bool
    """
    # references :
    # .. [#numpy_save] Numpy Developers. numpy 1.25 API reference : numpy.save. (consulted on september 05, 2023). https://numpy.org/doc/stable/reference/generated/numpy.save.html#numpy.save
    # .. [#numpy_npy] Numpy Developers. numpy 1.25 API reference : numpy.lib.format # NPY format. (consulted on september 05, 2023). https://numpy.org/doc/stable/reference/generated/numpy.lib.format.html#npy-format
    if type(so) != RDTrajectory :
        raise TypeError("so must be a SimulationOutput.")
    if not isstr(path) :
        raise TypeError("path must be a str.")
        
    json_path = filepath.append_extension_if_missing(path, ".json")
    data_path = filepath.remove_extension_if_existing(path, ".json") + "_data.npy"


    d = {
        "script" : rdscript_to_dict(so.script),
        "system" : rdsystem_to_dict(so.system),
        "data" : {"value" : filepath.get_last_element(data_path), "units" : str(so.data.units)},
        "t_sample" : unitarray_to_dict(so.t),
        "engine_description" : so.engine_description,
        "engine_option" : so.engine_option
        }
    
    if so.cgmap is not None :
        d["cgmap"] = list(so.cgmap)
    
    if separate_data :
        d["data"] = {"value" : filepath.get_last_element(data_path), "units" : str(so.data.units)}
        np.save(data_path, so.data.value)

    else :
        d["data"] = unitarray_to_dict(so.data)
        
    f = open(json_path, "w", encoding="utf-8")
    json.dump(d, f, indent = 4)

def load_rdtrajectory(path) :
    """
    Load a simulation output from a file.

    :param path: path of the output file to be created or replaced.
    :type path: str
    :return: loaded simulation output.
    :rtype: RDOutput
    """

    f = open(path, "r", encoding="utf-8")
    d = json.load(f)

    script = rdscript_from_dict(d["script"], base_path=filepath.get_base_path(path))
    system = rdsystem_from_dict(d["system"], base_path=filepath.get_base_path(path))
    data = unitarray_from_dict(d["data"], base_path=filepath.get_base_path(path))
    t_sample = unitarray_from_dict(d["t_sample"])
    engine_description = d["engine_description"]
    engine_option = d["engine_option"]
    cgmap = d.get("cgmap", None)
        
    return RDTrajectory(
        data = data, 
        t_sample = t_sample, 
        system = system,
        script = script,
        engine_description = engine_description,
        engine_option = engine_option,
        cgmap = cgmap
        )
