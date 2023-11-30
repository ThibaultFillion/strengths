import sys
import numpy
import ctypes
sys.path.append("../src/")
from strengths import *
from strengths.librdengine import *

# tests ###############################################

def test_build_reaction_rate_constant_array() :
    
    reactions = [Reaction("A + B -> C", 1, 2), 
                 Reaction("A + B -> C", 3, 4), 
                 Reaction("A + B -> C", 5, 6)]
    k = build_reaction_rate_constant_array(reactions, units_system=UnitsSystem())
    assert list(k) == [1,3,5]

    reactions = [Reaction("A -> C", "1 ms-1", 2), 
                 Reaction("A + B -> C", 3, 4), 
                 Reaction("A + B -> C", "5 m3/mol/min", 6)]
    k = build_reaction_rate_constant_array(reactions, units_system=UnitsSystem())
    assert list(k) == [1000, 3, UnitValue("5 m3/mol/min").convert("µm3/molecule/s").value]

    reactions = [Reaction("A + B -> C", 1, 2), 
                 Reaction("A + B -> C", 3, 4), 
                 Reaction("A + B -> C", 5, 6)]
    k = build_reaction_rate_constant_array(reactions, units_system=UnitsSystem(space="km", time="h", quantity="µmol"))
    assert list(k) == list(UnitArray([1,3,5], "µm3/molecule/s").convert("km3/µmol/h").value)
    
def test_build_reaction_environment_boolean_matrix() :
    environments = ["a", "b", "c"]
    reactions = [Reaction("A->"), 
                 Reaction("A->"), 
                 Reaction("A->")]
    m = build_reaction_environment_boolean_matrix(reactions, environments)
    assert list(m) == [1,1,1, 1,1,1, 1,1,1]

    environments = ["a", "b", "c"]
    reactions = [Reaction("A->", environments=[]), 
                 Reaction("A->"), 
                 Reaction("A->")]
    m = build_reaction_environment_boolean_matrix(reactions, environments)
    assert list(m) == [0,0,0, 1,1,1, 1,1,1]


    environments = ["a", "b", "c"]
    reactions = [Reaction("A->", environments=["a"]), 
                 Reaction("A->", environments=["b"]), 
                 Reaction("A->", environments=["c", "a"])]
    m = build_reaction_environment_boolean_matrix(reactions, environments)
    assert list(m) == [1,0,0, 0,1,0, 1,0,1]
    
    environments = ["c", "a", "b"]
    reactions = [Reaction("A->", environments=["a"]), 
                 Reaction("A->", environments=["b"]), 
                 Reaction("A->", environments=["c", "a"])]
    m = build_reaction_environment_boolean_matrix(reactions, environments)
    assert list(m) == [0,1,0, 0,0,1, 1,1,0]

def test_build_substrate_stoechiometric_matrix() :
    species = [
        Species("A"),
        Species("B"),
        Species("C"),
        Species("D")]
    reactions = [
        Reaction("A + 4 B -> C"),
        Reaction("2 A -> D"),
        Reaction("7 A + 8 B + 1 D -> C")]
    
    m = build_substrate_stoechiometric_matrix(species, reactions)
    assert list(m) == [
        1,2,7,
        4,0,8,
        0,0,0,
        0,0,1]

def test_build_stoechiometric_difference_matrix() :
    species = [
        Species("A"),
        Species("B"),
        Species("C"),
        Species("D")]
    reactions = [
        Reaction("A + 4 B -> C"),
        Reaction("2 A + 9 C-> 5 D + A + 6 C"),
        Reaction("7 A + 8 B + 1 D -> C")]
    
    m = build_stoechiometric_difference_matrix(species, reactions)
    assert list(m) == [
        -1 ,-1, -7,
        -4 , 0, -8,
         1 ,-3,  1,
         0 , 5, -1]

def test_build_diff_coef_environment_matrix() :

    environments = ["a", "b"]
    species = [
        Species("A"),
        Species("B"),
        Species("C"),
        Species("D")]

    m = build_diff_coef_environment_matrix(species, environments, units_system=UnitsSystem())
    assert list(m) == [
        0,0,
        0,0,
        0,0,
        0,0]
    
    environments = ["a", "b", "c"]
    species = [
        Species("A", D=1),
        Species("B", D=2),
        Species("C", D=3),
        Species("D", D=4)]

    m = build_diff_coef_environment_matrix(species, environments, units_system=UnitsSystem())
    assert list(m) == [
        1,1,1,
        2,2,2,
        3,3,3,
        4,4,4]


    environments = ["a", "b", "c"]
    species = [
        Species("A", D="1 m2/s"),
        Species("B", D=2),
        Species("C", D=3),
        Species("D", D=4)]

    m = build_diff_coef_environment_matrix(species, environments, units_system=UnitsSystem())
    assert list(m) == [
        UnitValue("1 m2/s").convert(UnitsSystem()).value, UnitValue("1 m2/s").convert(UnitsSystem()).value, UnitValue("1 m2/s").convert(UnitsSystem()).value,
        2,2,2,
        3,3,3,
        4,4,4]

    environments = ["a", "b", "c"]
    species = [
        Species("A", D=1),
        Species("B", D=2),
        Species("C", D=3),
        Species("D", D=4)]

    units_system=UnitsSystem(space="cm", time="min",  quantity="mmol")
    m = build_diff_coef_environment_matrix(species, environments, units_system)
    assert list(m) == [
        UnitValue("1 µm2/s").convert(units_system),
        UnitValue("1 µm2/s").convert(units_system),
        UnitValue("1 µm2/s").convert(units_system),
        
        UnitValue("2 µm2/s").convert(units_system),
        UnitValue("2 µm2/s").convert(units_system),
        UnitValue("2 µm2/s").convert(units_system),
        
        UnitValue("3 µm2/s").convert(units_system),
        UnitValue("3 µm2/s").convert(units_system),
        UnitValue("3 µm2/s").convert(units_system),
        
        UnitValue("4 µm2/s").convert(units_system),
        UnitValue("4 µm2/s").convert(units_system),
        UnitValue("4 µm2/s").convert(units_system)
        ]

    environments = ["a", "b", "c"]
    species = [
        Species("A", D={"a" : 1, "b" : 2, "c" : 3}),
        Species("B", D={"a" : 4, "b" : 5, "c" : 6}),
        Species("C", D={"a" : 7, "b" : 8, "c" : 9}),
        Species("D", D={"a" : 10, "b" : 11, "c" : 12})]

    units_system=UnitsSystem(space="cm", time="min",  quantity="mmol")
    m = build_diff_coef_environment_matrix(species, environments, units_system)
    assert list(m) == [
        UnitValue("1 µm2/s").convert(units_system),
        UnitValue("2 µm2/s").convert(units_system),
        UnitValue("3 µm2/s").convert(units_system),
        
        UnitValue("4 µm2/s").convert(units_system),
        UnitValue("5 µm2/s").convert(units_system),
        UnitValue("6 µm2/s").convert(units_system),
        
        UnitValue("7 µm2/s").convert(units_system),
        UnitValue("8 µm2/s").convert(units_system),
        UnitValue("9 µm2/s").convert(units_system),
        
        UnitValue("10 µm2/s").convert(units_system),
        UnitValue("11 µm2/s").convert(units_system),
        UnitValue("12 µm2/s").convert(units_system)
        ]

    environments = ["a", "b", "c"]
    species = [
        Species("A", D={"a" : 1, "b" : 2, "c" : 3}),
        Species("B", D={"a" : 4, "b" : "5 km2/h", "c" : 6}),
        Species("C", D={"a" : 7, "b" : 8, "c" : 9}),
        Species("D", D={"a" : 10, "b" : "11 mm2/min", "c" : 12})]

    units_system=UnitsSystem(space="cm", time="min",  quantity="mmol")
    m = build_diff_coef_environment_matrix(species, environments, units_system)
    assert list(m) == [
        UnitValue("1 µm2/s").convert(units_system),
        UnitValue("2 µm2/s").convert(units_system),
        UnitValue("3 µm2/s").convert(units_system),
        
        UnitValue("4 µm2/s").convert(units_system),
        UnitValue("5 km2/h").convert(units_system),
        UnitValue("6 µm2/s").convert(units_system),
        
        UnitValue("7 µm2/s").convert(units_system),
        UnitValue("8 µm2/s").convert(units_system),
        UnitValue("9 µm2/s").convert(units_system),
        
        UnitValue("10 µm2/s").convert(units_system),
        UnitValue("11 mm2/min").convert(units_system),
        UnitValue("12 µm2/s").convert(units_system)
        ]    

    environments = ["a", "b", "c"]
    species = [
        Species("A", D={"a" : 1, "b" : 2, "c" : 3}),
        Species("B", D={"a" : 4, "b" : 5, "c" : 6}),
        Species("C", D={"a" : 7, "b" : 8, "c" : 9}),
        Species("D", D={"a" : 10, "b" : 11, "c" : 12})]

    m = build_diff_coef_environment_matrix(species, environments, units_system=UnitsSystem())
    assert list(m) == [
        1,2,3,
        4,5,6,
        7,8,9,
        10,11,12
        ]

def test_make_ctypes_array() :
    
    a = [1,2,3]
    ca = make_ctypes_array(a, ctypes.c_int)
    for i in range(3) :
        assert a[i] == ca[i]

    a = [1,2,3]
    ca = make_ctypes_array(a, ctypes.c_double)
    for i in range(3) :
        assert a[i] == ca[i]

    a = [0,-2e-7,numpy.exp(3)]
    ca = make_ctypes_array(a, ctypes.c_double)
    for i in range(3) :
        assert a[i] == ca[i]
        
# run all ###############################################

def run_all_tests() :
    test_build_reaction_rate_constant_array()
    test_build_reaction_environment_boolean_matrix()
    test_build_substrate_stoechiometric_matrix()
    test_build_stoechiometric_difference_matrix()
    test_build_diff_coef_environment_matrix()
    test_make_ctypes_array()