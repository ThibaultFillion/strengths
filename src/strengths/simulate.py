from strengths.rdsystem import RDSystem, rdsystem_from_dict, rdsystem_to_dict
from strengths.units import *
from strengths.rdscript import RDScript
from strengths import engine_collection
from strengths.coarsegrain import coarsegrain_system, uncoarsegrain_trajectory

def simulate_script(
        script,
        engine,
        print_progress = False,
        cgmap = None
        ) :
    """
    Simulates the trajectory of a reaction diffusion system, using a given engine.

    :param script: the script of the simulation to be run.
    :type script: RDScript

    :param engine: the simulation engine that should handle the simulation.
        see the documentation or the engine_collection submodule for more information
        on the engines preinstalled with strengths.
    :type engine: RDSimulationEngineBase derived class

    :param print_progress: if true, the progression of the simulation (percentage) is printed
        at frequently.
    :type print_progress: bool
    
    :param cgmap: optionnal coarse graining index map.
    :type cgmap: array of int or None
    
    :return: system trajectory.
    :rtype: RDTrajectory
    """
    if cgmap is None :
        #initialization phase
        res = engine.setup(script)
    
        if res == 1 :
            raise Exception("Invalid option argument : \""+engine.get_option()+"\".")
        elif res == 2 :
            raise Exception("Invalid boundary conditions.")
            
        if print_progress :
            print("0 %", end="")
    
        # loop phase
        continue_simulation = True
        while continue_simulation :
            continue_simulation = engine.run(1000)
            if print_progress :
                print("\r" + str(engine.get_progress()) + " %", end="")
    
        # output phase
        output = engine.get_output()
    
        # finalize
        engine.finalize()
            
        return output
    else : 
        cgscript = script.copy()
        cgscript.system = coarsegrain_system(cgscript.system, cgmap)
        cgoutput = simulate_script(cgscript, engine, print_progress, None)
        output = uncoarsegrain_trajectory(cgoutput, script.system, cgmap)
        
        return output


def simulate(
        system,
        t_sample,
        engine = None,
        print_progress = False,
        cgmap = None,
        **script_keyword_arguments
        ) :
    """
    Simulates the trajectory of a reaction diffusion system.

    :param engine: the simulation engine that should handle the simulation.
        see the documentation or the engine_collection submodule for more information
        on the engines preinstalled with strengths. if None (default), engine_collection.default_engine()
        is used.
    :type engine: RDSimulationEngineBase derived class

    :param print_progress: if true, the progression of the simulation (percentage) is printed
        at frequently.
    :type print_progress: bool
    
    :param cgmap: optionnal coarse graining index map.
    :type cgmap: array of int or None
    
    Other parameters corresond to the remaining RDScript properties, and have the same default values :
        
    * time_step = 1e-3
    * sampling_policy = "on_t_sample"
    * sampling_interval = 1
    * t_max = "default"
    * rng_seed = None
    * units_system = UnitsSystem()

    :return: system trajectory
    :rtype: RDTrajectory
    """

    if engine is None :
        engine = engine_collection.default_engine() 

    d = dict(script_keyword_arguments)
    d["system"] = system
    d["t_sample"] = t_sample
    
    script = RDScript(**d)

    output = simulate_script(script = script, 
                             engine = engine, 
                             print_progress = print_progress,
                             cgmap = cgmap)
        
    return output
