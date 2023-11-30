from strengths.units import *
from strengths.rdsystem import RDSystem
from strengths.rdengine import RDEngineBase
from strengths.rdoutput import RDOutput
import strengths.kinetics
from strengths.typechecking import *

import numpy as np
import ctypes
import random
import time
import math

class KineticsRDEngine(RDEngineBase) :
    """
    Engine internally relying on the strengths.kinetics module.
    """

    def __init__(self) :

        self._incomplete = True
        
        super(KineticsRDEngine, self).__init__(option="", description="")

    def setup(self,
              script
                   # system,
                   # t_sample,
                   # *,
                   # time_step,
                   # t_max,
                   # sampling_policy,
                   # sampling_interval,
                   # rng_seed,
                   # units_system
                   ) :

        system = script.system
        t_sample = script.t_sample
        time_step = script.time_step
        t_max = script.t_max
        sampling_policy = script.sampling_policy
        sampling_interval = script.sampling_interval
        rng_seed = script.rng_seed
        units_system = script.units_system
        
        # #####################################################################################
        
        self._script = script
        self._units_system      = units_system 
        self._system            = system.copy()
        self._time_step         = UnitValue(time_step, Units(sys=units_system, dim=time_units_dimensions()))
        self._t_sample          = t_sample
        self._t_max             = UnitValue(t_max, Units(sys=units_system, dim=time_units_dimensions()))
        self._sampling_interval = UnitValue(sampling_interval, Units(sys=units_system, dim=time_units_dimensions()))
        self._sampling_policy   = sampling_policy
        self._state = self._system.state.convert(units_system)
        self._time = UnitValue(0, Units(sys=units_system, dim=time_units_dimensions()))
        self._sampled_states = []
        self._sampled_times = []
        self._last_sti_ratio = -1
        self._t_sample_pos = 0
        self._rng_seed = rng_seed
        
        self._perform_sampling_step()
        
    def run(self, breathe_dt) :
        
        t0 = time.time()
        
        while(1000*(time.time() - t0) < breathe_dt) :
            if not self.iterate() :
                break
            
        return self._incomplete

    def _perform_sampling_step(self) :
        
        if   self._sampling_policy == "on_t_sample" : 
            while self._t_sample_pos<len(self._t_sample) and self._time >= self._t_sample.get_at(self._t_sample_pos) :
                self.sample()
                self._t_sample_pos += 1
        
        elif self._sampling_policy == "on_iteration" : 
            self.sample()
            
        elif self._sampling_policy == "on_interval" : 
            tsi_ratio = math.floor((self._time/self._sampling_interval).value)
            if tsi_ratio > self._last_tsi_ratio :
                sample()
                self._last_sti_ratio = tsi_ratio
                
        elif self._sampling_policy == "on_sampling" : 
            pass
        
    
    def iterate(self) :
        if not self._incomplete :
            return self._incomplete
        
        self._state += strengths.kinetics.compute_dstatedt(self._system, self._state, self._units_system) * self._time_step

        self._time += self._time_step
        
        self._perform_sampling_step()
        
        if self._t_max >= 0 and self._time > self._t_max :  
            self._incomplete = False;
        
        return self._incomplete

    def iterate_n(self, n_iterations) :
        
        for i in range(n_iterations) : 
            self.iterate()

        return self._incomplete

    def get_progress(self) :
        
        progress = 0.0
        
        if self._t_max>0 :
            progress = 100*(self._time/self._t_max).value
            
        return progress
    
    def is_complete(self) : 
        
        return not self._incomplete
    
    def _count_samples(self) :
        
        return len(self._sampled_times)
    
    def _get_t_sample(self) :
        
        return UnitArray(self._sampled_times, Units(sys=self._units_system, dim=time_units_dimensions()))

    def _get_data(self) :
        
        return UnitArray(np.concatenate(tuple([s.convert(self._units_system).value for s in self._sampled_states])),
                         Units(sys=self._units_system, dim=quantity_units_dimensions()))
    
    def sample(self) : 
        
        self._sampled_states.append(self._state.copy())
        self._sampled_times.append(self._time.copy())
        
    def get_output(self) :
        
        return RDOutput(
            data = self._get_data(), 
            t_sample = self._get_t_sample(), 
            script = self._script,
            engine_description = self.description, 
            engine_option = self.option
            )
        
    def finalize(self) :
        
        pass
