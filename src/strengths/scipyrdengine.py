import scipy
import time as python_time
from strengths.sampler import create_sampler
from strengths.ode import build_ode_function
from strengths.units import *
from strengths.rdengine import *
from strengths.rdoutput import RDTrajectory

class ScipyRDEngine(RDEngineBase):
    """
    Engine relying on SciPy's ODE solvers
    """
    def __init__(self, integrator_class=scipy.integrate.LSODA):
        RDEngineBase.__init__(
            self,
            "",
            "IntegrateEngine"
            )
        self.integrator_class = integrator_class
        
    def setup(self, script):
        self._terminated = False
        self._script = script.copy()
        self._units_system = self._script.units_system.copy()
        self._t_max = self._script.t_max.convert(self._units_system).value
        self.sampler = create_sampler(
            self._script.sampling_policy,
            self._script.t_sample.convert(self._units_system).value,
            self._script.sampling_interval.convert(self._units_system).value,
            self._t_max
            )
        system = self._script.system.copy()
        ode_function = build_ode_function(
            system, 
            self._units_system
            )
        self.integrator = self.integrator_class(
            ode_function,
            t0 = 0,
            y0 = system.state.convert(self._units_system).value,
            t_bound = self._t_max
            )
        if self.sampler.requires_sample(0):
            self.sampler.sample(self.integrator.t, self.integrator.y)

    def run(self, breathe_dt):
        t_start = python_time.time_ns()
        while 1:
            self.iterate()
            if (python_time.time_ns()-t_start)/1e6 >= breathe_dt or self._terminated:
                break
        return self._ongoing()

    def iterate(self):
        if self._terminated:
            return False
        s = self.integrator.step()
        if self.integrator.status == "failed":
            print(s)
            self._terminated = True
            return False
        if self.sampler.requires_sample(self.integrator.t):
            self.sampler.sample(self.integrator.t, self.integrator.y)
        if self.integrator.status == "finished":
            self._terminated = True
        return self._ongoing()
        
    def iterate_n(self, n_iterations):
        for i in range(n_iterations):
            self.iteration()
        return self._ongoing()

    def get_progress(self):
        return 100*self.integrator.t/self._t_max

    def sample(self):
        self.sampler.sample(self.integrator.t, self.integrator.y)

    def get_output(self):
        data = UnitArray(
            np.array(self.sampler.x).flatten(), 
            Units(self._units_system, quantity_units_dimensions())
            )
        output = RDTrajectory(
            data = data,
            t_sample = UnitArray(
                self.sampler.t, 
                Units(self._units_system, time_units_dimensions())
                ),
            system = self._script.system.copy(),
            script = self._script.copy(),
            engine_description = self.description,
            engine_option = self.option
            )
        return output
        
    def finalize(self):
        pass
    
    def _ongoing(self):
        return (self.integrator.t<=self._t_max) and (not self._terminated)
