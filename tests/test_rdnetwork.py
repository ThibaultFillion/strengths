import sys
sys.path.append("../src/")
from strengths import *

def make_test_rdn() :
    d = {
        "units" : {"space":"nm", "time":"min", "quantity":"kmol"},
        "species" : 
            [
            {"label" : "A"},
            {"label" : "B"},
            {"label" : "C"}
            ],
        "reactions" :
            [
            {"stoechiometry": "A -> B ", "label" : "A"},
            {"stoechiometry": "A -> B ", "label" : "B"},
            {"stoechiometry": "A -> B ", "label" : "C"}
            ]
        }
    return rdnetwork_from_dict(d)

# tests ######################################################

def test_rdn_species_labels() :
    rn = make_test_rdn()
    assert rn.get_species("A").label == "A"
    assert rn.get_species("B").label == "B"
    assert rn.get_species("C").label == "C"

def test_rdn_reactions_labels() : 
    rn = make_test_rdn()
    assert rn.get_reaction("A").label == "A"
    assert rn.get_reaction("B").label == "B"
    assert rn.get_reaction("C").label == "C"

def test_rdn_get_reaction_index_from_number() :
    rn = make_test_rdn()
    assert rn.get_reaction_index(0) == 0
    assert rn.get_reaction_index(1) != 0
    assert rn.get_reaction_index(1) == 1
    assert rn.get_reaction_index(1.5) == 1
    assert rn.get_reaction_index(-1) == None
    assert rn.get_reaction_index(3) == None

def test_rdn_get_reaction_index_from_str() :
    rn = make_test_rdn()
    assert rn.get_reaction_index("A") == 0
    assert rn.get_reaction_index("B") != 0
    assert rn.get_reaction_index("B") == 1
    assert rn.get_reaction_index("D") == None

def test_rdn_get_reaction_index_from_reaction() :
    rn = make_test_rdn()
    assert rn.get_reaction_index(Reaction("->", 1, 1, label="A")) == 0
    assert rn.get_reaction_index(Reaction("->", 1, 1, label="A")) != 1
    assert rn.get_reaction_index(Reaction("->", 1, 1, label="B")) == 1
    assert rn.get_reaction_index(Reaction("->", 1, 1, label="D")) == None

def test_rdn_get_spêcies_index_from_number() :
    rn = make_test_rdn()
    assert rn.get_species_index(0) == 0
    assert rn.get_species_index(1) != 0
    assert rn.get_species_index(1) == 1
    assert rn.get_species_index(1.5) == 1
    assert rn.get_species_index(-1) == None
    assert rn.get_species_index(3) == None

def test_rdn_get_reaction_index_from_str() :
    rn = make_test_rdn()
    assert rn.get_species_index("A") == 0
    assert rn.get_species_index("B") != 0
    assert rn.get_species_index("B") == 1
    assert rn.get_species_index("D") == None

def test_rdn_get_reaction_index_from_str() :
    rn = make_test_rdn()
    assert rn.get_species_index(Species("A")) == 0
    assert rn.get_species_index(Species("A")) != 1
    assert rn.get_species_index(Species("B")) == 1
    assert rn.get_species_index(Species("D")) == None

def test_reaction_kf_default_units_dimensions() :
    r = Reaction("A + 2 B -> C", kf=1, kr=1)
    assert r.kf_units_dimensions() == UnitsDimensions(space=6, time=-1, quantity=-2)
    assert r.kr_units_dimensions() == UnitsDimensions(space=0, time=-1, quantity=0)

    r = Reaction("3 A + 2 B -> 2 C", kf=1, kr=1)
    assert r.kf_units_dimensions() == UnitsDimensions(space=12, time=-1, quantity=-4)
    assert r.kr_units_dimensions() == UnitsDimensions(space=3,  time=-1, quantity=-1)

def test_reaction_split() :    
    r = Reaction("A -> 2 B", 
                 kf=1.5, 
                 kr=0.3, 
                 label="test", 
                 environments=["a"], 
                 units_system=UnitsSystem(space="dm", time="ns", quantity="nmol"))
    rf, rr = r.split()

    assert rf.to_string().strip() == "A -> 2 B"
    assert rf.kf == r.kf
    assert rf.kr == 0
    assert rf.environments == r.environments
    assert isnone(rf.label)
    assert rf.units_system == r.units_system

    assert rr.to_string().strip() == "2 B -> A"
    assert rr.kf == r.kr
    assert rr.kr == 0
    assert rr.environments == r.environments
    assert isnone(rr.label)
    assert rr.units_system == r.units_system

def test_reaction_dict() :
    
    rd = {"eq" : "A+B->C", "k+" : 1, "k-" : 1, "env" : ["a", "b"], "label":"label"}  
    
    r = reaction_from_dict(rd)
    
    assert r.kf == UnitValue("1 µm3/molecule/s")
    assert r.kr == UnitValue("1 s-1")
    assert r.to_string().strip() == "A + B -> C"
    assert r.environments == ("a", "b")
    assert r.label == "label"
    
    rd = reaction_to_dict(r)

    assert rd["stoechiometry"].strip() == "A + B -> C"
    assert UnitValue(rd["k+"]) == UnitValue("1 µm3.molecule-1.s-1")
    assert UnitValue(rd["k-"]) == UnitValue("1 s-1")
    assert rd["environments"] == ("a", "b")
    assert rd["units"] == {"space" : "µm", "time" : "s", "quantity" : "molecule"}
    assert rd["label"] == "label"

def test_reaction_dict_default() :
    
    rd = reaction_to_dict(Reaction("A+B->C", kf=10))

    assert rd["stoechiometry"].strip() == "A + B -> C"
    assert UnitValue(rd["k+"]) == UnitValue("10 µm3.molecule-1.s-1")
    assert UnitValue(rd["k-"]) == UnitValue("0 s-1")
    assert rd["environments"] == None
    assert rd["units"] == {"space" : "µm", "time" : "s", "quantity" : "molecule"}
    assert rd["label"] == None

def test_species_dict() :
    
    sd = {"label" : "A", 
          "density" : {"a" : 1, "b" : "2 µM"}, 
          "D" : {"a" : "5 m2/min", "b" : 0}, 
          "chstt" : {"a" : True, "b" : False}
          }  
    
    s = species_from_dict(sd)
    
    assert s.label == "A"
    assert s.density == {"a" : UnitValue("1 molecule/µm3"), "b" : UnitValue("2 µM")}
    assert s.D       == {"a" : UnitValue("5 m2/min"), "b" : UnitValue("0 µm2/s")}
    assert s.chstt == {"a" : True, "b" : False}
    
    sd = species_to_dict(s)
    
    assert sd["label"] == "A"
    assert sd["density"] == {"a" : str(UnitValue("1 molecule/µm3")), "b" : str(UnitValue("2 µM"))}
    assert sd["D"]       == {"a" : str(UnitValue("5 m2/min")), "b" : str(UnitValue("0 µm2/s"))}
    assert sd["chstt"]   == {"a" : True, "b" : False}
    assert sd["units"] == {"space" : "µm", "time" : "s", "quantity" : "molecule"}

def test_species_dict_default() :
    
    sd = species_to_dict(Species("A"))
    
    assert sd["label"] == "A"
    assert sd["density"] == str(UnitValue("0 molecule/µm3"))
    assert sd["D"]       == str(UnitValue("0 µm2/s"))
    assert sd["chstt"]   == False
    assert sd["units"] == {"space" : "µm", "time" : "s", "quantity" : "molecule"}

def test_get_env_index() :
    d = {
        "units" : {"space":"nm", "time":"min", "quantity":"kmol"},
        "species" : 
            [
            {"label" : "A"}
            ],
        "reactions" : [],
        "environments" : ["a", "b", "e", "g"]
        }
        
    rdn = rdnetwork_from_dict(d)
    
    assert rdn.get_environment_index(-1) is None
    assert rdn.get_environment_index(0) == 0
    assert rdn.get_environment_index(0.5) == 0
    assert rdn.get_environment_index(2) == 2
    assert rdn.get_environment_index(3) == 3
    assert rdn.get_environment_index(3.9) == 3
    assert rdn.get_environment_index(4) is None
    
    assert rdn.get_environment_index("f") is None
    assert rdn.get_environment_index("a") == 0
    assert rdn.get_environment_index("g") == 3
    assert rdn.get_environment_index("b") == 1

# run all #############################################################

def run_all_tests() : 
    test_rdn_species_labels() 
    test_rdn_reactions_labels() 
    test_rdn_get_reaction_index_from_number() 
    test_rdn_get_reaction_index_from_str() 
    test_rdn_get_reaction_index_from_reaction() 
    test_rdn_get_spêcies_index_from_number() 
    test_rdn_get_reaction_index_from_str() 
    test_rdn_get_reaction_index_from_str() 
    test_reaction_kf_default_units_dimensions() 
    test_reaction_split()
    test_reaction_dict()
    test_reaction_dict_default()
    test_species_dict()
    test_species_dict_default()
    test_get_env_index()