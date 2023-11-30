import sys
import numpy
sys.path.append("../src/")
from strengths import *

# tests ###############################################

def generate_rds() :
    rds = {
        "network" : {
            "species" : [
                {"label" : "A", "density" : 1000}, 
                {"label" : "B"}
                ], 
            "reactions":[
                {"eq" : "A -> B", "k+" : 1, "k-" : 1}
                ]
            },
        "space" : {"w" : 1, "h" : 1, "d" : 1}
        }
    return rdsystem_from_dict(rds)   
    
def test_save_output_load_output() :
    rds = generate_rds()
    out = simulate(rds, t_sample=np.linspace(0, 100, 100), time_step=0.1)
    
    save_rdtrajectory(out, "test_output_files/out.json")
    out2 = load_rdtrajectory("test_output_files/out.json")

    save_rdtrajectory(out, "test_output_files/out2.json", separate_data=False)
    out3 = load_rdtrajectory("test_output_files/out2.json")

def test_output_traj_t0_matching_system_state() :
    rds = generate_rds()
    out = simulate(rds, t_sample=np.linspace(0, 100, 100), time_step=0.1)
    assert numpy.allclose(list(out.data.value[0:out.system.state_size()]), 
                          list(out.system.state.value))
    
    out = simulate(rds, t_sample=np.linspace(0, 100, 100), time_step=0.1, units_system=UnitsSystem(quantity="mol"))
    assert numpy.allclose(list(out.data.value[0:out.system.state_size()]), 
                          list(out.system.state.convert("mol").value))
    
    out = simulate(rds, t_sample=np.linspace(0, 100, 100), time_step=0.1, engine=engine_collection.gillespie_engine(),units_system=UnitsSystem(quantity="mol"))
    assert numpy.allclose(list(out.data.value[0:out.system.state_size()]), 
                          list(out.system.state.convert("mol").value))
    
    out = simulate(rds, t_sample=np.linspace(0, 100, 100), time_step=0.1, engine=engine_collection.tauleap_engine(),units_system=UnitsSystem(quantity="mol"))
    assert numpy.allclose(list(out.data.value[0:out.system.state_size()]), 
                          list(out.system.state.convert("mol").value))
    
    out = simulate(rds, t_sample=np.linspace(0, 100, 100), time_step=0.1, units_system=UnitsSystem(quantity="mmol", time="ms", space="m"))
    assert numpy.allclose(list(out.data.value[0:out.system.state_size()]), 
                          list(out.system.state.convert("mmol").value))

    out = simulate(rds, t_sample=UnitArray(np.linspace(0, 100, 100), "s"), time_step=0.1, units_system=UnitsSystem(quantity="µmol", time="ms", space="m"))
    assert numpy.allclose(list(out.data.value[0:out.system.state_size()]), 
                          list(out.system.state.convert("µmol").value))
    
# run all ###############################################

def run_all_tests() :
    test_save_output_load_output()
    test_output_traj_t0_matching_system_state()