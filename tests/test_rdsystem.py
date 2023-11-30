import sys
sys.path.append("../src/")
from strengths import *

class SomeCoordClass :
    def __init__(self,x,y,z) :
        self.x = x
        self.y = y
        self.z = z

def make_test_rds_1() :
    d_rds = {"network" : {
        "species" : 
            [
            {"label" : "A", "density" : "10.5 molecule/µm3", "chstt" : {
                 "cytoplasm" : True,
                 "membrane" : 1,
                 "outside" : 0}
                 },
            {"label" : "B", "density" : {
                    "cytoplasm" : "10 molecule/µm3",
                    "membrane" : "5 molecule/µm3",
                    "outside" : "2 molecule/µm3"},
                "chstt" : {
                    "cytoplasm" : False,
                    "membrane" : 1.0,
                    "outside" : 0.0}
                }
            ],
        "reactions" : 
            [
            {"label" : "r", "stoechiometry" : "A -> ", "k+" : "1 s-1", "k-" : 1}
            ],
        "environments" : 
            [
            "membrane",
            "cytoplasm",
            "outside"
            ]
        },
        "space" : {
            "w" : 8,
            "h" : 8,
            "d" : 8,
            "cell_volume" : "10 µm3",
            "cell_env"    : 0 }
        }
    return rdsystem_from_dict(d_rds)

def make_test_rds_2() :
    d_rds = {"network" : {
        "species" : 
            [
            {"label" : "A", "density" : "10.5 molecule/µm3", "chstt" : {
                 "cytoplasm" : True,
                 "membrane" : 1,
                 "outside" : 0}
                 },
            {"label" : "B", "density" : {
                    "cytoplasm" : "10 molecule/µm3",
                    "membrane" : "5 molecule/µm3",
                    "outside" : "2 molecule/µm3"},
                "chstt" : {
                    "cytoplasm" : False,
                    "membrane" : 1.0,
                    "outside" : 0.0}
                }
            ],
        "reactions" : 
            [
            {"label" : "r", "stoechiometry" : "A -> ", "k+" : "1 s-1", "k-" : 1}
            ],
        "environments" : 
            [
            "membrane",
            "cytoplasm",
            "outside"
            ]
        },
        "space" : {
                "w" : 6,
                "h" : 1,
                "d" : 1,
                "cell_volume" : "10 µm3",
                "cell_env"    : [0,0,1,1,2,2] }
        }
    return rdsystem_from_dict(d_rds)

# tests ############################################

def test_rds_space_size() :
    rds = make_test_rds_1()
    assert rds.space.size() == 8*8*8

def test_generate_space_cell_env() :
    rds = make_test_rds_1()
    for i in range(rds.space.size()) :
        assert rds.space.cell_env[i] == 0

def test_rds_space_cell_env_set_from_dict() :
    rds = make_test_rds_2()
    assert list(rds.space.cell_env) == [0,0,1,1,2,2]
    
def test_generate_species_state_default() :
    rds = make_test_rds_1()
    cell_state = generate_species_state(rds.network.get_species("A"), rds.network, rds.space, rds.units_system)
    for i in range(rds.space.size()) :
        assert cell_state.value[i] == 105
    cell_state = generate_species_state(rds.network.get_species("B"), rds.network, rds.space, rds.units_system)
    for i in range(rds.space.size()) :
        assert cell_state.value[i] == 50

    rds = make_test_rds_2()
    cell_state = list(generate_species_state(rds.network.get_species("A"), rds.network, rds.space, rds.units_system).value)
    assert cell_state == [105,105,105,105,105,105]

    cell_state = list(generate_species_state(rds.network.get_species("B"), rds.network, rds.space, rds.units_system).value)
    assert cell_state == [50,50,100,100,20,20]

def test_generate_species_chstt_map_default() :
    rds = make_test_rds_1()
    cell_chstt = generate_species_chemostats(rds.network.get_species("A"), rds.network, rds.space)
    for i in range(rds.space.size()) :
        assert cell_chstt[i] == True
    cell_chstt = generate_species_chemostats(rds.network.get_species("B"), rds.network, rds.space)
    for i in range(rds.space.size()) :
        assert cell_chstt[i] == True

def test_default_state() :
    rds = {
        "network" : {"species" : [{"label" : "A"}], "reactions":[]},
        "space" : {"w" : 2, "h" : 3, "d" : 4}
        }
    rds = rdsystem_from_dict(rds)
    assert list(rds.state.value) == [0.0, 0.0,
                                     0.0, 0.0,
                                     0.0, 0.0,
                                     
                                     0.0, 0.0,
                                     0.0, 0.0,
                                     0.0, 0.0,
                                                                      
                                     0.0, 0.0,
                                     0.0, 0.0,
                                     0.0, 0.0,

                                     0.0, 0.0,
                                     0.0, 0.0,
                                     0.0, 0.0
                                     ]

    rds = {
        "network" : {"species" : [{"label" : "A", "density" : 3.0}], "reactions":[]},
        "space" : {"w" : 2, "h" : 3, "d" : 4}
        }
    rds = rdsystem_from_dict(rds)
    assert list(rds.state.value) == [3.0, 3.0,
                                     3.0, 3.0,
                                     3.0, 3.0,
                                     
                                     3.0, 3.0,
                                     3.0, 3.0,
                                     3.0, 3.0,
                                                                      
                                     3.0, 3.0,
                                     3.0, 3.0,
                                     3.0, 3.0,

                                     3.0, 3.0,
                                     3.0, 3.0,
                                     3.0, 3.0
                                     ]

    rds = {
        "network" : {"species" : [{"label" : "A", "density" : 3.0}, {"label" : "B", "density" : 0.5}], "reactions":[]},
        "space" : {"w" : 2, "h" : 3, "d" : 4}
        }
    rds = rdsystem_from_dict(rds)
    assert list(rds.state.value) == [3.0, 3.0,
                                     3.0, 3.0,
                                     3.0, 3.0,
                                     
                                     3.0, 3.0,
                                     3.0, 3.0,
                                     3.0, 3.0,
                                                                      
                                     3.0, 3.0,
                                     3.0, 3.0,
                                     3.0, 3.0,

                                     3.0, 3.0,
                                     3.0, 3.0,
                                     3.0, 3.0,
                                     
                                     0.5, 0.5,
                                     0.5, 0.5,
                                     0.5, 0.5,
                                    
                                     0.5, 0.5,
                                     0.5, 0.5,
                                     0.5, 0.5,
                                                                     
                                     0.5, 0.5,
                                     0.5, 0.5,
                                     0.5, 0.5,

                                     0.5, 0.5,
                                     0.5, 0.5,
                                     0.5, 0.5
                                     ]

def test_set_state() :
    rds = {
        "network" : {"species" : [{"label" : "A"}], "reactions":[]},
        "space" : {"w" : 2, "h" : 3, "d" : 4}
        }
    rds = rdsystem_from_dict(rds)
    rds.set_state("A", 5, 3.0)
    assert list(rds.state.value) == [0.0, 0.0,
                                     0.0, 0.0,
                                     0.0, 3.0,
                                     
                                     0.0, 0.0,
                                     0.0, 0.0,
                                     0.0, 0.0,
                                                                      
                                     0.0, 0.0,
                                     0.0, 0.0,
                                     0.0, 0.0,

                                     0.0, 0.0,
                                     0.0, 0.0,
                                     0.0, 0.0
                                     ]

    rds = {
        "network" : {"species" : [{"label" : "A"}], "reactions":[]},
        "space" : {"w" : 2, "h" : 3, "d" : 4}
        }
    rds = rdsystem_from_dict(rds)
    rds.set_state("A", (0,1,2), 3.0)
    assert list(rds.state.value) == [0.0, 0.0,
                                     0.0, 0.0,
                                     0.0, 0.0,
                                     
                                     0.0, 0.0,
                                     0.0, 0.0,
                                     0.0, 0.0,
                                                                      
                                     0.0, 0.0,
                                     3.0, 0.0,
                                     0.0, 0.0,

                                     0.0, 0.0,
                                     0.0, 0.0,
                                     0.0, 0.0
                                     ]

    rds = {
        "network" : {"species" : [{"label" : "A"}], "reactions":[]},
        "space" : {"w" : 2, "h" : 3, "d" : 4}
        }
    rds = rdsystem_from_dict(rds)
    rds.set_state("A", SomeCoordClass(0,1,2), 3.0)
    assert list(rds.state.value) == [0.0, 0.0,
                                     0.0, 0.0,
                                     0.0, 0.0,
                                     
                                     0.0, 0.0,
                                     0.0, 0.0,
                                     0.0, 0.0,
                                                                      
                                     0.0, 0.0,
                                     3.0, 0.0,
                                     0.0, 0.0,

                                     0.0, 0.0,
                                     0.0, 0.0,
                                     0.0, 0.0
                                     ]
    
# run all tests ############################################

def run_all_tests() : 
    test_rds_space_size()
    test_generate_space_cell_env()
    test_rds_space_cell_env_set_from_dict()
    test_generate_species_state_default()
    test_generate_species_chstt_map_default()
    test_default_state() 
    test_set_state() 
