from strengths.units import *
from strengths.rdscript import RDScript
from strengths.typechecking import *

class RDEngineBase :
    """
    base class for a reaction-diffusion system trajectory simulation engine.

    :param option: engine option (default = "").
    :type option: str
    :param description: description of the engine (ie. its name, and which algorithms it implements) that will be featured in the simulation outputs for context (default = "").
    :type description: str        
    """

    def __init__(self, option = "", description = "") :
        """
        constructor.
        """
        
        self._option = option
        self._description = description

    @property
    def option(self):
        """
        Option currently associated with the engine (str).
        Should be set before calling initialize_3D.
        """
        
        return self._option

    @property
    def description(self):
        """
        Description of the engine (str).
        """
        
        return self._description
    
    def setup(self, script) :
        """
        setup the engine/ initialize the simulation.
        
        :param script: Simulation parameters
        :type script:  RDScript
        """
        
        raise NotImplementedError("")

    def run(self, breathe_dt) :
        """
        Makes the engine keep iterating for a real-time duration of *breathe_dt* or until the simulation is finished.
        Should be called after initialize_3D and before get_output and finalize.
        
        :param breathe_dt: maximal real-time duration for which the engine should keep running in ms.
        :type breathe_dt: int
        :returns: True if the simulation is incomplete, False otherwise.
        :rtype: bool
        """
        
        raise NotImplementedError("")

    def iterate(self) :
        """
        Makes the engine iterate once.
        Should be called after initialize_3D and before get_output and finalize.
        
        :returns: True if the simulation is incomplete, False otherwise.
        :rtype: bool
        """
        
        raise NotImplementedError("")


    def iterate_n(self, n_iterations) :
        """
        Makes the engine iterate *n_itrations* times or until the simulation is finished.
        Should be called after initialize_3D and before get_output and finalize.
        
        :param n_iterations: number of iterations requested.
        :type n_iterations: int
        :returns: True if the simulation is incomplete, False otherwise.
        :rtype: bool
        """
        
        raise NotImplementedError("")

    def get_progress(self) :
        """
        Returns the percentage of simulation progress.
        Should be called after initialize_3D and before finalize.
        
        :rtype: number
        """
        
        raise NotImplementedError("")

    def sample(self) :
        """
        sample the current system state.
        """
        
        raise NotImplementedError("")
        
    def is_complete(self) : 
        """
        Tells if the simulation is complete.

        :returns: True if the simulation is complete, False otherwise.
        :rtype: bool
        """
        
        raise NotImplementedError("")
    
    def get_output(self) :
        """
        Returns the system trajectory data array.
        Should be called after initialize_3D and before finalize.

        :rtype: array
        """
        
        raise NotImplementedError("")

    def finalize(self) :
        """
        Finalize the engine use.
        Should be called last.
        """
        
        raise NotImplementedError("")

