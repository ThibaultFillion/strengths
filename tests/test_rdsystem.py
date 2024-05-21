import sys
sys.path.append("../src/")
from strengths import *
from strengths.constants import avogadro_number

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

def _test_rdnetwork_apply_reaction(spacetype):

    spacedict = {
        "w" : 3
        }

    if spacetype == "graph" :
        spacedict = {
            "type" : "graph",
            "nodes" : [{},{},{}],
            "edges" : []
            }

    system = rdsystem_from_dict({
        "network" : {
            "species" : [
                {"label" : "D", "density" : 0},
                {"label" : "A", "density" : 5},
                {"label" : "B", "density" : 5},
                {"label" : "C", "density" : 0}
                ],
            "reactions" : [
                {"eq" : "A + B -> C", "label" : "1"},
                {"eq" : "C -> D",     "label" : "2"}
                ]
            },
        "space" : spacedict
        })

    # state = [0,0,0,  5,5,5,  5,5,5,  0,0,0]

    # basic case with default arguments
    state = system.apply_reaction("1")
    assert list(state.value) == [0,0,0,  4,5,5,  4,5,5,  1,0,0]
    assert list(system.state.value) == [0,0,0,  5,5,5,  5,5,5,  0,0,0]

    # negative number
    state = system.apply_reaction("1", n=-2)
    assert list(state.value) == [0,0,0,  7,5,5,  7,5,5,  -2,0,0]

    # float number
    state = system.apply_reaction("1", n=1.5)
    assert list(state.value) == [0,0,0,  3.5,5,5,  3.5,5,5,  1.5,0,0]

    # chaning reaction and position
    state = system.apply_reaction("2", position=2)
    assert list(state.value) == [0,0,1,  5,5,5,  5,5,5,  0,0,-1]

    # reaction by index
    state = system.apply_reaction(1, position=2)
    assert list(state.value) == [0,0,1,  5,5,5,  5,5,5,  0,0,-1]

    # chaning number
    state = system.apply_reaction("2", position=2, n=3)
    assert list(state.value) == [0,0,3,  5,5,5,  5,5,5,  0,0,-3]

    # update
    system.apply_reaction("2", position=2, update=True)
    assert list(system.state.value) == [0,0,1,  5,5,5,  5,5,5,  0,0,-1]

    # custom state
    state = system.apply_reaction("2", position=2,
        state=UnitArray([0,0,6,  0,0,0,  0,0,0,  0,0,2], "molecule"))
    assert list(state.value) == [0,0,7,  0,0,0,  0,0,0,  0,0,1]

    # custom state and chemostats
    state = system.apply_reaction("2", position=2,
        state=UnitArray([0,0,6,  0,0,0,  0,0,0,  0,0,2], "molecule"),
        chemostats=[0,0,0,  0,0,0,  0,0,0,  0,0,1])
    assert list(state.value) == [0,0,7,  0,0,0,  0,0,0,  0,0,2]

    #different units
    n=1e23
    state = system.apply_reaction("2", position=2, n=n,
        state=UnitArray([0,0,6,  0,0,0,  0,0,0,  0,0,2], "mol"))
    assert list(state.value) == [0,0,6+n/avogadro_number(),  0,0,0,  0,0,0,  0,0,2-n/avogadro_number()]

def test_rdnetwork_apply_reaction__graph():
    _test_rdnetwork_apply_reaction("graph")

def test_rdnetwork_apply_reaction__grid():
    _test_rdnetwork_apply_reaction("grid")