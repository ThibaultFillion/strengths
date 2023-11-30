import sys
sys.path.append("../src/")
from strengths import *

# tests ###############################################

def test_load_rds_multifile_unspecified_units_systems() :
    sys = load_rdsystem("test_json_files/1/system.json")
    assert list(sys.space.cell_env) == [0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8]
    assert sys.space.cell_vol == UnitValue("1 µm3")
    
    assert sys.network.species[0].density == UnitValue("2 molecule/µm3")
    assert sys.network.species[1].density == UnitValue("4 molecule/µm3")
    assert sys.network.species[2].density == UnitValue("6 M")
    
    assert sys.network.species[0].D == UnitValue("3 µm2/s")
    assert sys.network.species[1].D == UnitValue("5 µm2/s")
    assert sys.network.species[2].D == UnitValue("1 cm2/min")
    
    assert sys.network.reactions[0].kf == UnitValue("1 µm3/molecule/s")
    assert sys.network.reactions[1].kf == UnitValue("3 µm3/molecule/s")
    assert sys.network.reactions[2].kf == UnitValue("5 M-1.s-1")

    assert sys.network.reactions[0].kr == UnitValue("2 s-1")
    assert sys.network.reactions[1].kr == UnitValue("4 s-1")
    assert sys.network.reactions[2].kr == UnitValue("6 min-1")

def test_load_rds_multifile_units_system_inherited_from_rds_only() :
    sys = load_rdsystem("test_json_files/2/system.json")
    assert list(sys.space.cell_env) == [0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8]
    assert sys.space.cell_vol == UnitValue("1 m3")
    
    assert sys.network.species[0].density == UnitValue("2 mol/m3")
    assert sys.network.species[1].density == UnitValue("4 mol/m3")
    assert sys.network.species[2].density == UnitValue("6 M")
    
    assert sys.network.species[0].D == UnitValue("3 m2/h")
    assert sys.network.species[1].D == UnitValue("5 m2/h")
    assert sys.network.species[2].D == UnitValue("1 cm2/min")
    
    assert sys.network.reactions[0].kf == UnitValue("1 m3/mol/h")
    assert sys.network.reactions[1].kf == UnitValue("3 m3/mol/h")
    assert sys.network.reactions[2].kf == UnitValue("5 M-1.s-1")

    assert sys.network.reactions[0].kr == UnitValue("2 h-1")
    assert sys.network.reactions[1].kr == UnitValue("4 h-1")
    assert sys.network.reactions[2].kr == UnitValue("6 min-1")

def test_load_rds_multifile_fully_heterogenous_units_systems() :
    sys = load_rdsystem("test_json_files/3/system.json")
    assert list(sys.space.cell_env) == [0,1,2,3,4,5,6,7,8,0,1,2,3,4,5,6,7,8]
    assert sys.space.cell_vol == UnitValue("1 cm3")
    
    assert sys.network.species[0].density == UnitValue("2 pmol/pm3")
    assert sys.network.species[1].density == UnitValue("4 µmol/cm3")
    assert sys.network.species[2].density == UnitValue("6 M")
    
    assert sys.network.species[0].D == UnitValue("3 pm2/ps")
    assert sys.network.species[1].D == UnitValue("5 cm2/µs")
    assert sys.network.species[2].D == UnitValue("1 cm2/min")
    
    assert sys.network.reactions[0].kf == UnitValue("1 nm3/kmol/s")
    assert sys.network.reactions[1].kf == UnitValue("3 fm3/cmol/min")
    assert sys.network.reactions[2].kf == UnitValue("5 M-1.s-1")

    assert sys.network.reactions[0].kr == UnitValue("2 s-1")
    assert sys.network.reactions[1].kr == UnitValue("4 min-1")
    assert sys.network.reactions[2].kr == UnitValue("6 min-1")
                
# run all ###############################################

def run_all_tests() :
    test_load_rds_multifile_unspecified_units_systems()
    test_load_rds_multifile_units_system_inherited_from_rds_only()
    test_load_rds_multifile_fully_heterogenous_units_systems()