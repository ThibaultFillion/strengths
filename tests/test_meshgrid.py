import sys
sys.path.append("../src/")
from strengths import *

def make_test_rdspace() : 
    cell_vol = 1.3
    d = {
        "w" : 2,
        "h" : 3,
        "d" : 4,
        "cell_volume" : cell_vol
        }
    return rdspace_from_dict(d)

class CoordLike:
    def __init__(self, x=0, y=0, z=0):
        self.x = x
        self.y = y
        self.z = z
    def __eq__(self, v) : 
        return (v.x == self.x and 
                v.y == self.y and 
                v.z == self.z)

# tests ##########################################################

def test_mg_dimensions_type() : 
    mg = make_test_rdspace()
    assert type(mg.w) == int and type(mg.h) == int and type(mg.d) == int

def test_mg_dimensions_value() : 
    mg = make_test_rdspace()
    assert mg.w == 2 and mg.h == 3 and mg.d == 4

def test_mg_size() : 
    mg = make_test_rdspace()
    assert mg.size() == 4*2*3

def test_mg_len() : 
    mg = make_test_rdspace()
    assert len(mg.cell_env) == mg.size()

def test_mg_cell_env_type() :
    mg = make_test_rdspace()
    assert isarray(mg.cell_env)

def test_mg_cell_env_array_dtype() : 
    mg = make_test_rdspace()
    assert mg.cell_env.dtype == int

def test_mg_cell_vol() :
    mg = make_test_rdspace()
    assert mg.cell_vol == UnitValue(1.3, "µm3")

def test_mg_get_cell_coordinates_as_tuple():
    mg = make_test_rdspace()
    assert mg.get_cell_coordinates(1) == (1,0,0)
    assert mg.get_cell_coordinates(1.5) == (1,0,0)
    assert mg.get_cell_coordinates(0.5) == (0,0,0)
    assert mg.get_cell_coordinates(2.5) == (0,1,0)

def test_mg_get_cell_coordinates_as_coordlike():
    mg = make_test_rdspace()
    assert mg.get_cell_coordinates(1, return_type=CoordLike) == CoordLike(1,0,0)
    assert mg.get_cell_coordinates(1.5, return_type=CoordLike) == CoordLike(1,0,0)
    assert mg.get_cell_coordinates(0.5, return_type=CoordLike) == CoordLike(0,0,0)
    assert mg.get_cell_coordinates(2.5, return_type=CoordLike) == CoordLike(0,1,0)

def test_mg_get_cell_index_from_int():
    mg = make_test_rdspace()
    assert mg.get_cell_index(2) == 2
    
def test_mg_get_cell_index_from_float():
    mg = make_test_rdspace()
    assert mg.get_cell_index(2.5) == 2

def test_mg_get_cell_index_from_tuple():
    mg = make_test_rdspace()
    assert mg.get_cell_index((1,0,0)) == 1
    assert mg.get_cell_index((0,1,0)) == 2
    assert mg.get_cell_index((0,0,1)) == 2*3

def test_mg_get_cell_index_from_coord_like():
    mg = make_test_rdspace()
    assert mg.get_cell_index(CoordLike(1,0,0)) == 1
    assert mg.get_cell_index(CoordLike(0,1,0)) == 2
    assert mg.get_cell_index(CoordLike(0,0,1)) == 2*3

def test_periodicity():
    mg = rdspace_from_dict({"w":3,"h":3,"d":3})
    assert(mg.are_neighbors((0,0,0), (2,0,0))==False)

    mg = rdspace_from_dict({"w":3,"h":3,"d":3})
    assert(mg.are_neighbors((1,0,0), (1,2,0))==False)    

    mg = rdspace_from_dict({"w":3,"h":3,"d":3,"boundary_conditions" : {"x" : "periodical"}})
    assert(mg.are_neighbors((0,0,0), (2,0,0))==True)

    mg = rdspace_from_dict({"w":3,"h":3,"d":3,"boundary_conditions" : {"y" : "periodical"}})
    assert(mg.are_neighbors((1,0,0), (1,2,0))==True)    
    
# run all #############################################################

def run_all_tests() : 
    test_mg_dimensions_type() 
    test_mg_dimensions_value() 
    test_mg_size() 
    test_mg_len()
    test_mg_cell_env_type()
    test_mg_cell_env_array_dtype()
    test_mg_cell_vol()
    test_mg_get_cell_coordinates_as_tuple()
    test_mg_get_cell_coordinates_as_coordlike()
    test_mg_get_cell_index_from_int()
    test_mg_get_cell_index_from_float()
    test_mg_get_cell_index_from_tuple()
    test_mg_get_cell_index_from_coord_like()
    test_periodicity()