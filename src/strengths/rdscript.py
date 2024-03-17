from strengths.units import *
from strengths.rdsystem import RDSystem, rdsystem_from_dict, rdsystem_to_dict, load_rdsystem, save_rdsystem
from strengths.typechecking import *

import json
import copy
import random

class RDScript :
    """
    Class that gather simulation parameters.
    Each argument sets the property with the same name.
    Default values are : 
        
    * time_step : 1e-3,
    * t_max : "default",
    * sampling_policy : "on_t_sample",
    * sampling_interval : 1,
    * rng_seed : None,
    * units_system : UnitsSystem()
    """
    
    def __init__(self,
                system,
                t_sample,
                time_step = 1e-3,
                t_max = "default",
                sampling_policy = "on_t_sample",
                sampling_interval = 1,
                rng_seed = None,
                units_system = UnitsSystem()
                ) :
        
            self.units_system = units_system
            self.system = system
            self.t_sample = t_sample
            self.time_step = time_step
            self.t_max = t_max
            self.sampling_policy = sampling_policy
            self.sampling_interval = sampling_interval
            self.rng_seed = rng_seed

    @property
    def system(self) :
        """
        reaction diffusion system of which the time trajectory should be simulated.
        """
        
        return self._system
    
    @system.setter
    def system(self, system) :        
        if type(system) != RDSystem :
            raise TypeError("system must be a DRSystem object.")
        self._system = system.copy()
    
    @property
    def t_sample(self) :
        """
        Expected sampling times (used if sampling_policy=="on_t_sample").
        """
        
        return self._t_sample


    @t_sample.setter
    def t_sample(self, t_sample) :
        self._t_sample = UnitArray(t_sample, Units(sys=self.units_system, dim=time_units_dimensions()), convert=False)
    
    @property
    def time_step(self) :
        """
        Time step to be used, if necessary.
        """
        
        return self._time_step

    @time_step.setter
    def time_step(self, time_step) :        
        self._time_step = UnitValue(time_step, Units(sys=self.units_system, dim=time_units_dimensions()), convert=False)
        
    @property
    def t_max(self) :
        """
        Time at which the simulation should stop.
        if set to "default", t_max will be the last element of t_sample.
        """
        
        if self._t_max=="default" :
            return self.t_sample.get_at(len(self.t_sample)-1)
        else :
            return self._t_max.copy()

    @t_max.setter
    def t_max(self, t_max) :
        if t_max == "default" :
            self._t_max = t_max
        else :            
            self._t_max = UnitValue(t_max, Units(sys=self.units_system, dim=time_units_dimensions()), convert=False)
        
    @property
    def sampling_policy(self) :
        """
        String specifying how sampling should be handled by the simulation engine.
        accepted sampling policies are :
            
        * "on_t_sample" : (default) sampling is done accorging to t_sample (as close as possible).
        * "on_iteration" : sampling is done at t=0, then at every iteration.
        * "on_interval" : sampling is done at t=0, then at the given time interval.
        * "no_sampling" : no sampling is managed by the engine. the sample() method should be used for manual sampling.
        """
        
        return self._sampling_policy

    @sampling_policy.setter
    def sampling_policy(self, sampling_policy) :
        
        if not isstr(sampling_policy) : 
            raise TypeError("sampling_policy must be a string")
        
        if sampling_policy not in ["on_t_sample" , 
                                   "on_iteration", 
                                   "on_interval" , 
                                   "no_sampling" ] :
            raise ValueError("\""+sampling_policy+"\" is not a valid sampling policy. accepted values are : \"on_t_sample\", \"on_iteration\", \"on_interval\" and \"no_sampling\".")
        
        self._sampling_policy = sampling_policy

    @property
    def sampling_interval(self) :
        """
        In-simulation time interval at which the system state should be sampled,
        if the sampling policy is "on_interval".
        """
        
        return self._sampling_interval

    @sampling_interval.setter
    def sampling_interval(self, sampling_interval) :
        self._sampling_interval = UnitValue(sampling_interval, Units(sys=self.units_system, dim=time_units_dimensions()), convert=False)
    
    
    @property
    def rng_seed(self) : 
        """
        Seed to be used for the engine's pseudo random number generator, if any.
        """

        return self._rng_seed

    @rng_seed.setter
    def rng_seed(self, rng_seed) :         
        if rng_seed == None :
            rng_seed = random.randint(0, (2**32)-1)
        else :
            rng_seed = int(rng_seed)
            
        self._rng_seed = rng_seed

    @property
    def units_system(self) :
        """
        Default units system to be used for numerical inputs. This is also the units system in which the system
        trajectory (sampled states and times) will be expressed.
        """

        return self._units_system

    @units_system.setter
    def units_system(self, units_system) :
        
        if isdict(units_system) :
            self._units_system = unitssystem_from_dict(units_system)
        elif type(units_system) == UnitsSystem :
            self._units_system = units_system.copy()
        else :
            raise TypeError("units_system must be a dict or an instance of UnitsSystem.")
    
    def copy(self) : 
        """
        Returns a deep copy of the object.
        """
        
        return copy.deepcopy(self)

def rdscript_from_dict(d, base_path=None) :
    """
    Creates a RDScript object from a dictionary.
    """
    
    d = valproc.process_input_dict_keys(d, [
        ["system"],
        ["t_sample"],
        ["time_step", "time step", "dt"],
        ["t_max", "tmax"],
        ["sampling_policy", "sampling policy"],
        ["sampling_interval", "sampling interval"],
        ["rng_seed", "rng seed", "seed"],
        ["units", "units_system", "units system", "u"]
        ])
    
    da = {}

    da["units_system"] = valproc.retrive_units_system_from_dict(
        d=d,
        default="default",
        parent_units_system=UnitsSystem()
        )
    
    if "system" in d :
        rds = d["system"]
        
        if isdict(rds) : 
            rds = rdsystem_from_dict(rds, da["units_system"], base_path)
        elif isstr(rds) : 
            rds = filepath.get_path_with_base(rds, base_path)
            rds = load_rdsystem(rds, da["units_system"])
        else :
            raise TypeError("system must be a dictionary.")
        da["system"] = rds
    else :
        raise ValueError("script system is missing")

    if "t_sample" in d :
        t_sample = d["t_sample"]
        
        if isdict(t_sample) : 
            t_sample = unitarray_from_dict(t_sample, base_path=base_path)
        da["t_sample"] = t_sample
    else :
        raise ValueError("script t_sample is missing")

    if "time_step"         in d : da["time_step"]         = d["time_step"]
    if "t_max"             in d : da["t_max"]             = d["t_max"]
    if "sampling_policy"   in d : da["sampling_policy"]   = d["sampling_policy"]
    if "sampling_interval" in d : da["sampling_interval"] = d["sampling_interval"]
    if "rng_seed"          in d : da["rng_seed"]          = d["rng_seed"]
    
    return RDScript(**da)
        
def rdscript_to_dict(script) :
    """
    Creates a dictionary from a RDScript object.
    """
    
    d = {
        "system"            : rdsystem_to_dict(script.system),
        "t_sample"          : unitarray_to_dict(script.t_sample),
        "time_step"         : str(script.time_step),
        "t_max"             : str(script.t_max),
        "sampling_policy"   : script.sampling_policy,
        "sampling_interval" : str(script.sampling_interval),
        "rng_seed"          : script.rng_seed,
        "units"             : unitssystem_to_dict(script.units_system)
        }
    
    return d

def load_rdscript(path) :
    """
    Loads a RDScript object from a JSON file.
    """
    
    f = open(path, "r", encoding="utf-8")
    d = json.load(f)
    f.close()
    return rdscript_from_dict(d, base_path=filepath.get_base_path(path))

def save_rdscript(script) :
    """
    Saves a RDScript object as a JSON file.
    """
    
    d = rdscript_to_dict(script)
    f = open(path, "w", encoding="utf-8")
    json.dumps(d, f, indent = 4)
    f.close()
