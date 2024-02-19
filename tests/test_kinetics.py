import sys
sys.path.append("../src")
from strengths import *
from strengths.constants import avogadro_number
import strengths.kinetics as strnkin
import numpy as np

# test kintetics module

##################################
# test reaction rate
# should be done both for rdnetworks and graphs

def test_compute_reaction_rates_dimensions_for_different_reaction_orders_with_default_units_default_units(spacetype) :
    #checks that rates are in molecule/s, no mater the reaction order
 
    system_dict = {
        "network" : {
            "species" : [
                    {"label" : "A"},
                    {"label" : "B"},
                    {"label" : "C"},
                    {"label" : "D"},
                    {"label" : "E"},
                    {"label" : "F"}
                ],
            "reactions" : [
                    {"eq" : " -> ", "kf" : 2,  "kr" : 3},
                    {"eq" : "0 A + 0 B  -> 0 C", "kf" : 2,  "kr" : 3},
                    {"eq" : "A -> B", "kf" : 2,  "kr" : 3},
                    {"eq" : "A + B -> C", "kf" : 2,  "kr" : 3},
                    {"eq" : "A + B + C-> C + D + E + F", "kf" : 2,  "kr" : 3},
                    {"eq" : "A + B -> C", "kf" : 2,  "kr" : 3},
                    {"eq" : "2 C + 4 B -> 6 C + 0 E + 2 F", "kf" : 2,  "kr" : 3}
                ]
            }
        }
    if spacetype=="graph" :
        system_dict["space"] = {
            "type" : "graph",
            "nodes" : [{}],
            "edges" : []}

    system = rdsystem_from_dict(system_dict)

    state = UnitArray([1,1,1,1,1,1], "molecule")

    for i in range(len(system.network.reactions)) :
        rates = strnkin.compute_reaction_rates(
            system = system,
            reaction = system.network.reactions[i],
            position = 0,
            state = state,
            units_system = UnitsSystem())
        for rate in rates :
            assert rate.units.dim == UnitsDimensions(quantity=1, space=0, time=-1)

def test_compute_reaction_rates_value(spacetype) :
    #checks that rates computed values are correct
    kf = 12
    kr = 0.5
    V = 1.23e-3

    A = 5
    B = 90
    C = 1
    D = 548
    E = 0.33
    F = 1000
    
    NA = avogadro_number()
    state = UnitArray([A, B, C, D, E, F], "mol")    

    system_dict = {
        "network" : {
            "species" : [
                    {"label" : "A"},
                    {"label" : "B"},
                    {"label" : "C"},
                    {"label" : "D"},
                    {"label" : "E"},
                    {"label" : "F"}
                ],
            "reactions" : [
                    {"eq" : "2 C + 4 B + 1 A -> 6 D + 0 E + 3 F", "kf" : str(kf)+" µm18/molecule6/min",  "kr" : kr}
                ]
            }
        }

    if spacetype=="graph" :
        system_dict["space"] = {
            "type" : "graph",
            "nodes" : [{"vol" : str(V)+" mL"}],
            "edges" : []}
    else :
        system_dict["space"] = {"cell_vol" : str(V)+" mL"}
        
    system = rdsystem_from_dict(system_dict)

    V *= 1e12
    rf = kf/60 * V * NA*A/V * (NA*C/V)**2 * (NA*B/V)**4 # molecule/s
    rr = kr * V * (NA*D/V)**6 * (NA*F/V)**3             # molecule/s
    
    rates = strnkin.compute_reaction_rates(
        system = system,
        reaction = system.network.reactions[0],
        position = 0,
        state = state,
        units_system = UnitsSystem())

    assert np.isclose(rates[0].convert("molecule/s").value, rf)
    assert np.isclose(rates[1].convert("molecule/s").value, rr)

def test_compute_reaction_rates_with_environments(spacetype) :
    #checks that reaction env are taken into account

    state = UnitArray([1, 1, 1, 1], "molecule")        
    system_dict = {
        "network" : {
            "environments" : ["a", "b"],
            "species" : [
                    {"label" : "A"},
                    {"label" : "B"}
                ],
            "reactions" : [
                    {"eq" : "A -> B", "kf":1, "kr":1, "env":["b"]}
                ]
            }
        }

    if spacetype=="graph" :
        system_dict["space"] = {
            "type" : "graph",
            "nodes" : [{"env" : 0}, {"env" : 1}],
            "edges" : []}
    else :
        system_dict["space"] = {"w":2, "environments":[0, 1]}

    system = rdsystem_from_dict(system_dict)
    
    rates = strnkin.compute_reaction_rates(
        system = system,
        reaction = system.network.reactions[0],
        position = 0,
        state = state,
        units_system = UnitsSystem())

    assert rates[0].convert("molecule/s").value == 0
    assert rates[1].convert("molecule/s").value == 0

    rates = strnkin.compute_reaction_rates(
        system = system,
        reaction = system.network.reactions[0],
        position = 1,
        state = state,
        units_system = UnitsSystem())

    assert rates[0].convert("molecule/s").value == 1
    assert rates[1].convert("molecule/s").value == 1

def test_compute_reaction_rates_in_inhomogeneous_state(spacetype) :
    #checks that reaction env are taken into account

    state = UnitArray([1, 2, 3, 4], "molecule")        
    system_dict = {
        "network" : {
            "species" : [
                    {"label" : "A"},
                    {"label" : "B"}
                ],
            "reactions" : [
                    {"eq" : "A -> B", "kf":1, "kr":1}
                ]
            }
        }

    if spacetype=="graph" :
        system_dict["space"] = {
            "type" : "graph",
            "nodes" : [{}, {}],
            "edges" : []}
    else :
        system_dict["space"] = {"w":2}
        
    system = rdsystem_from_dict(system_dict)
    
    rates = strnkin.compute_reaction_rates(
        system = system,
        reaction = system.network.reactions[0],
        position = 0,
        state = state,
        units_system = UnitsSystem())

    assert rates[0].convert("molecule/s").value == 1
    assert rates[1].convert("molecule/s").value == 3

    rates = strnkin.compute_reaction_rates(
        system = system,
        reaction = system.network.reactions[0],
        position = 1,
        state = state,
        units_system = UnitsSystem())

    assert rates[0].convert("molecule/s").value == 2
    assert rates[1].convert("molecule/s").value == 4

#####################################

def test_compute_diffusion_rates(spacetype) :

    system_dict = {
        "network" : {
            "species" : [{"label" : "B"},
                         {"label" : "A", "D" : 1},
                         {"label" : "C"}],
            "reactions" : []
            }
        }
    
    if spacetype == "graph" :
        system_dict["space"] = {
            "type" : "graph",
            "nodes" : [{},{}],
            "edges" : [{"nodes" : [0, 1], "surface" : 1, "distance" : 1}]}
    else :
        system_dict["space"] = {"w" : 2}
    
    system = rdsystem_from_dict(system_dict)
    #                  B  B  A  A  C  C
    #                  0  1  0  1  0  1   
    state = UnitArray([0, 0, 30, 40, 0, 0], "molecule")

    rates = strnkin.compute_diffusion_rates(
        system = system,
        species = "A",
        src_position = 0,
        dst_position = 1,
        state = state,
        units_system = UnitsSystem()
        )
    
    assert rates[0].value == 30
    assert rates[1].value == 40

def test_compute_diffusion_rates_with_different_env(spacetype) :

    system_dict = {
        "network" : {
            "environments" : ["a", "b"],
            "species" : [{"label" : "A", "D" : {"a" : 1, "b" : 0}}],
            "reactions" : []
            }
        }
    
    if spacetype == "graph" :
        system_dict["space"] = {
            "type" : "graph",
            "nodes" : [{"env":0},{"env":1}],
            "edges" : [{"nodes" : [0, 1], "surface" : 1, "distance" : 1}]}
    else :
        system_dict["space"] = {"w" : 2, "cell_env" : [0, 1]}
    
    system = rdsystem_from_dict(system_dict)
    state = UnitArray([30, 40], "molecule")

    rates = strnkin.compute_diffusion_rates(
        system = system,
        species = "A",
        src_position = 0,
        dst_position = 1,
        state = state,
        units_system = UnitsSystem()
        )

    assert rates[0].value == 0
    assert rates[1].value == 0

def test_compute_diffusion_rates_with_different_units(spacetype) :

    system_dict = {
        "network" : {
            "species" : [{"label" : "A", "D" : "1 m2/h"}],
            "reactions" : []
            }
        }
    
    if spacetype == "graph" :
        system_dict["space"] = {
            "type" : "graph",
            "nodes" : [{"vol" : "1 nm3"},{"vol" : "1 nm3"}],
            "edges" : [{"nodes" : [0, 1], "surface" : "1 nm2", "distance" : "1 nm"}]}
    else :
        system_dict["space"] = {"w" : 2, "cell_vol" : "1 nm3"}
    
    system = rdsystem_from_dict(system_dict)
    #                  B  B  A  A  C  C
    #                  0  1  0  1  0  1   
    state = UnitArray([30, 40], "mol")

    rates = strnkin.compute_diffusion_rates(
        system = system,
        species = "A",
        src_position = 0,
        dst_position = 1,
        state = state,
        units_system = UnitsSystem()
        )

    assert np.isclose(rates[0].value, 1e18/(60*60) * 30 * avogadro_number())
    assert np.isclose(rates[1].value, 1e18/(60*60) * 40 * avogadro_number())

#####################################

def test_compute_dspeciesdt(spacetype) :
    
    system_dict = {
        "network" : {
            "species" : [
                {"label" : "A", "D" : 1},
                {"label" : "B", "D" : 0.5}
                ],
            "reactions" : [
                {"eq" : "A -> B", "kf" : 1, "kr" : 1}
                ]
            }
        }

    if spacetype == "graph" : 
        system_dict["space"] = {
            "type" : "graph",
            "nodes" : [{}, {}],
            "edges" : [{"nodes" : [0,1]}]
            }
    else :
        system_dict["space"] = {
            "w" : 2
            }

    system = rdsystem_from_dict(system_dict)
    #                  A   A   B   B
    #                  0   1   0   1
    state = UnitArray([10, 20, 1,  2], "molecule")
    
    dA = strnkin.compute_dspeciesdt(
        system = system,
        species = "A",
        position = 0,
        state = state,
        units_system = UnitsSystem()
        )

    assert dA.value == (-10*1 + 1*1   -   10*1 + 20*1)
    assert dA.units.dim == UnitsDimensions(space=0, time=-1, quantity=1)

def test_compute_dspeciesdt_with_boundary_conditions() :
    
    system_dict = {
        "network" : {
            "species" : [
                {"label" : "A", "D" : 1},
                {"label" : "B", "D" : 0.5}
                ],
            "reactions" : [
                {"eq" : "A -> B", "kf" : 1, "kr" : 1}
                ]
            },
        "space" : {
            "w" : 2,
            "boundary_conditions" : {"x" : "periodical"}
            }
        }

    system = rdsystem_from_dict(system_dict)
    #                  A   A   B   B
    #                  0   1   0   1
    state = UnitArray([10, 20, 1,  2], "molecule")
    
    dA = strnkin.compute_dspeciesdt(
        system = system,
        species = "A",
        position = 0,
        state = state,
        units_system = UnitsSystem()
        )

    assert dA.value == (-10*1 + 1*1   +   2*(-10*1 + 20*1))
    
    #############"
    
    system_dict = {
        "network" : {
            "species" : [
                {"label" : "A", "D" : 1},
                {"label" : "B", "D" : 0.5}
                ],
            "reactions" : [
                {"eq" : "A -> B", "kf" : 1, "kr" : 1}
                ]
            },
        "space" : {
            "w" : 2,
            "boundary_conditions" : {"y" : "periodical"}
            }
        }

    system = rdsystem_from_dict(system_dict)
    #                  A   A   B   B
    #                  0   1   0   1
    state = UnitArray([10, 20, 1,  2], "molecule")
    
    dA = strnkin.compute_dspeciesdt(
        system = system,
        species = "A",
        position = 0,
        state = state,
        units_system = UnitsSystem()
        )

    assert dA.value == (-10*1 + 1*1   +   (-10*1 + 20*1))

def run_all_tests() :
  test_compute_reaction_rates_dimensions_for_different_reaction_orders_with_default_units_default_units("grid")
  test_compute_reaction_rates_value("grid")
  test_compute_reaction_rates_with_environments("grid")
  test_compute_reaction_rates_in_inhomogeneous_state("grid")
  
  test_compute_reaction_rates_dimensions_for_different_reaction_orders_with_default_units_default_units("graph")
  test_compute_reaction_rates_value("graph")
  test_compute_reaction_rates_with_environments("graph")
  test_compute_reaction_rates_in_inhomogeneous_state("graph")
  
  ############
  
  test_compute_diffusion_rates("grid")
  test_compute_diffusion_rates_with_different_env("grid")
  test_compute_diffusion_rates_with_different_units("grid")
  
  test_compute_diffusion_rates("graph")
  test_compute_diffusion_rates_with_different_env("graph")
  test_compute_diffusion_rates_with_different_units("graph")
  
  ############
  
  test_compute_dspeciesdt("grid")
  test_compute_dspeciesdt_with_boundary_conditions()
  
  test_compute_dspeciesdt("graph")
