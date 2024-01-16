import platform
import strengths
import ctypes
import pathlib

from strengths.librdengine import LibRDEngine

def _get_strengths_path():
        path = ""
        for i in strengths.__path__ :
            path += i
        return path

def _get_engine_path() :
    # Assumes the engine library binary file is the 
    # only non directory item in the strengths_engine 
    # directory.

    strnpath = _get_strengths_path()
    sedirpath = strnpath + "/engines/strengths_engine"
    path = pathlib.Path(sedirpath)
    content = path.iterdir()
    files = []
    for i in content :
        if not i.is_dir() : 
            files.append(i)
    if len(files) == 0 :
        raise RuntimeError("It seems that the engine shared library is missing. Maybe it needs to be compiled form source. For more information on how to build the package from source, please refer to the documentation.")
    return str(files[0])

def gillespie_engine():
    """
    Engine using the original Gillespie algorithm (Gillespie, 1977) [#Gillespie1977]_.
    Diffusion is treated as a first order reaction according to Bernstein's method (Bernstein, 2005) [#Bernstein2005]_.
    """
    # references :
    # .. [#Bernstein2005] Bernstein, D. (2005). Simulating mesoscopic reaction-diffusion systems using the Gillespie algorithm. Physical Review E, 71(4), Article 041103. https://doi.org/10.1103/PhysRevE.71.041103
    # .. [#Gillespie1977] Gillespie, D. T. (1977). Exact stochastic simulation of coupled chemical reactions. The Journal of Physical Chemistry, 81(25), 2340-2361. https://doi.org/10.1021/j100540a008
    
    path = _get_engine_path()
    return LibRDEngine(
        ctypes.CDLL(path),
        option="gillespie",
        description="description",
        requires_molecules=True
        )

def tauleap_engine():  
    """
    Engine using the Gillespie tau leap method (Gillespie, 2001) [#Gillespie2001]_ with a static time step.
    Diffusion is treated as a first order reaction according to Bernstein's method (Bernstein, 2005) [#Bernstein2005]_.
    """
    # references :
    # .. [#Bernstein2005] Bernstein, D. (2005). Simulating mesoscopic reaction-diffusion systems using the Gillespie algorithm. Physical Review E, 71(4), Article 041103. https://doi.org/10.1103/PhysRevE.71.041103
    # .. [#Gillespie2001] Gillespie, D. T. (2001). Approximate accelerated stochastic simulation of chemically reacting systems. The Journal of Chemical Physics, 115(4), 1716-1733. https://doi.org/10.1063/1.1378322
    
    path = _get_engine_path()
    return LibRDEngine(
        ctypes.CDLL(path),
        option = "tauleap",
        description="description",
        requires_molecules=True
        )

def euler_engine():
    """
    Engine using a simple Euler method with a static time step.
    Diffusion is treated as a first order reaction according to Bernstein's method (Bernstein, 2005) [#Bernstein2005]_.
    """
    # references :
    # .. [#Bernstein2005] Bernstein, D. (2005). Simulating mesoscopic reaction-diffusion systems using the Gillespie algorithm. Physical Review E, 71(4), Article 041103. https://doi.org/10.1103/PhysRevE.71.041103

    path = _get_engine_path()
    return LibRDEngine(
        ctypes.CDLL(path),
        option = "euler",
        description = "description",
        requires_molecules=False
        )

def default_engine():
    """
    Default engine (among those above) used by functions such as simulate.
    """
        
    return euler_engine()
