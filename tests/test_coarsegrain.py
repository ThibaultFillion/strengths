import sys
import numpy
import ctypes
sys.path.append("../src/")
from strengths import *
from strengths.librdengine import *
from strengths.coarsegrain import *
# tests ###############################################

def test_grid_to_graph() :
    grid = RDGridSpace(w=2, h=3, d=4)
    graph = grid_to_graph(grid)
    
    assert grid.size() == graph.size()
    assert len(graph.edges) == 7*4 + 6*3
    for node in graph.nodes : 
        assert node.volume == grid.cell_vol
    
    grid = RDGridSpace(w=2, h=2, d=2, cell_env=[1,2,3,4,5,6,7,8])
    graph = grid_to_graph(grid)
    
    assert grid.size() == graph.size()
    assert len(graph.edges) == 2*4 + 4
    assert graph.get_edge(0, 1) is not None
    assert graph.get_edge(0, 2) is not None
    assert graph.get_edge(1, 3) is not None
    assert graph.get_edge(3, 2) is not None
    
    assert graph.get_edge(4, 5) is not None
    assert graph.get_edge(4, 6) is not None
    assert graph.get_edge(5, 7) is not None
    assert graph.get_edge(7, 6) is not None
    
    assert graph.get_edge(0, 4) is not None
    assert graph.get_edge(1, 5) is not None
    assert graph.get_edge(2, 6) is not None
    assert graph.get_edge(3, 7) is not None

    assert graph.get_cell_env_array() == list([1,2,3,4,5,6,7,8])

def test_grid_coarsegrain() :
    
    index_map = [
        0,0,1,1,1,
        0,0,0,1,1,
        0,0,0,2,2,
        3,3,2,2,2,
        3,3,2,2,2,
        ]
    env = [
        0,0,1,1,1,
        0,0,0,1,1,
        0,0,0,2,2,
        3,3,2,2,2,
        3,3,2,2,2,
        ]
    grid = RDGridSpace(w=5, h=5, d=1, cell_env=env, cell_vol="1 m3")
    graph = coarsegrain_grid(grid, index_map)
    
    assert graph.size() == 4
    
    assert graph.nodes[0].environment == 0
    assert graph.nodes[1].environment == 1
    assert graph.nodes[2].environment == 2
    assert graph.nodes[3].environment == 3

    assert np.isclose(graph.nodes[0].volume.convert("m3").value, 8)
    assert np.isclose(graph.nodes[1].volume.convert("m3").value, 5)
    assert np.isclose(graph.nodes[2].volume.convert("m3").value, 8)
    assert np.isclose(graph.nodes[3].volume.convert("m3").value, 4)
    
    assert graph.get_edge(0, 1) is not None
    assert graph.get_edge(0, 3) is not None
    assert graph.get_edge(0, 2) is not None
    assert graph.get_edge(3, 2) is not None
    assert graph.get_edge(1, 2) is not None
    assert graph.get_edge(1, 3) is None

    assert np.isclose(graph.get_edge(0, 1).surface.convert("m2").value, 3)
    assert np.isclose(graph.get_edge(0, 3).surface.convert("m2").value, 2)
    assert np.isclose(graph.get_edge(0, 2).surface.convert("m2").value, 2)
    assert np.isclose(graph.get_edge(3, 2).surface.convert("m2").value, 2)
    assert np.isclose(graph.get_edge(1, 2).surface.convert("m2").value, 2)


    x0 = (0+1+2+0+1+2+0+1)/8
    y0 = (0+1+2+0+1+2+1+2)/8
    
    x1 = (2+3+4+3+4)/5
    y1 = (0+0+1+0+1)/5
    
    d01 = ((x0-x1)**2 + (y0-y1)**2)**(1/2)
    
    assert np.isclose(graph.get_edge(0, 1).distance.convert("m").value, d01)    

def test_grid_coarsegrain_rm() :
    space = RDGridSpace(2,2,1, [2,1,4,4], 1)
    cgspace = coarsegrain_grid(space, [0,1,-1,-1])
    assert cgspace.size() == 2
    assert list(cgspace.get_cell_env_array()) == [2,1]
    assert list(cgspace.get_cell_vol_array().value) == [1,1]

def test_system_coarsegrain() :
    system = RDSystem(
        network = RDNetwork(
            species = [
                Species("A"),
                Species("B")
                ],
            reactions = [
                ],
            environments = [
                "a","b","c","d","e"
                ]
            ),
        space = RDGridSpace(
            w = 3,
            h = 3,
            cell_env = [0,1,1,
                        0,2,4,
                        3,2,4]
            ),
        state = [
            1,2,3,
            4,5,6,
            7,8,9,
                 
            2,3,4,
            5,6,7,
            8,9,10
            ],
        chemostats = [
            1,0,0, 
            0,1,0, 
            0,1,0,
            
            0,1,0, 
            0,0,0, 
            0,0,0
            ]
        )
    
    cgmap = [0,1,1,
             0,2,4,
             3,2,4]
    
    cgsystem = coarsegrain_system(system, cgmap)  
    assert cgsystem.space.size() == 5
    assert list(cgsystem.space.get_cell_env_array()) == [0,1,2,3,4]
    assert list(cgsystem.space.get_cell_vol_array().value) == [2,2,2,1,2]
    
    assert list(cgsystem.state.value) == [5,5,13,7,15,
                                          7,7,15,8,17]
    
    assert list(cgsystem.chemostats) == [1,0,1,0,0,
                                         0,1,0,0,0]    

def test_system_coarsegrain_rm() :
    system = RDSystem(
        network = RDNetwork(
            species = [
                Species("A"),
                Species("B")
                ],
            reactions = [
                ],
            environments = [
                "a","b","c","d","e"
                ]
            ),
        space = RDGridSpace(
            w = 3,
            h = 3,
            cell_env = [0,1,1,
                        0,2,4,
                        3,2,4]
            ),
        state = [
            1,2,3,
            4,5,6,
            7,8,9,
                 
            2,3,4,
            5,6,7,
            8,9,10
            ],
        chemostats = [
            1,0,0, 
            0,1,0, 
            0,1,0,
            
            0,1,0, 
            0,0,0, 
            0,0,0
            ]
        )
    
    cgmap = [0,1,1,
             0,2,-1,
             3,2,-1]
    
    cgsystem = coarsegrain_system(system, cgmap)  
    assert cgsystem.space.size() == 4
    assert list(cgsystem.space.get_cell_env_array()) == [0,1,2,3]
    assert list(cgsystem.space.get_cell_vol_array().value) == [2,2,2,1]
    
    assert list(cgsystem.state.value) == [5,5,13,7,
                                          7,7,15,8]
    
    assert list(cgsystem.chemostats) == [1,0,1,0,
                                         0,1,0,0]    

def test_ucg_traj_rm() :
    system = RDSystem(
        network = RDNetwork(
            species = [
                Species("A"),
                Species("B")
                ],
            reactions = [
                ],
            environments = [
                "a","b","c","d","e"
                ]
            ),
        space = RDGridSpace(
            w = 3,
            h = 3,
            cell_env = [0,1,1,
                        0,2,4,
                        3,2,4]
            ),
        state = [
            1,2,2,
            1,3,5,
            4,3,5,
                 
            2,3,3,
            2,4,6,
            5,4,6
            ],
        chemostats = [
            1,0,0, 
            0,1,0, 
            0,1,0,
            
            0,1,0, 
            0,0,0, 
            0,0,0
            ]
        )
    
    cgmap = [0,1,1,
             0,2,-1,
             3,2,-1]
    
    cgsystem = coarsegrain_system(system, cgmap)  
    data = UnitArray([cgsystem.state.value[i%(4*2)] for i in range(4*2*4)], cgsystem.state.units)
    traj = RDTrajectory(data, UnitArray([0,1,2,3], "s"), cgsystem)  
    ucgtraj = uncoarsegrain_trajectory(traj, system, cgmap)
    
    assert list(ucgtraj.data.value) == [
            1,2,2,
            1,3,0,
            4,3,0,                 
            2,3,3,
            2,4,0,
            5,4,0,
            
            1,2,2,
            1,3,0,
            4,3,0,                 
            2,3,3,
            2,4,0,
            5,4,0,
            
            1,2,2,
            1,3,0,
            4,3,0,                 
            2,3,3,
            2,4,0,
            5,4,0,
            
            1,2,2,
            1,3,0,
            4,3,0,                 
            2,3,3,
            2,4,0,
            5,4,0
            ]
        
# run all ###############################################

def run_all_tests() :
    test_grid_to_graph()
    test_grid_coarsegrain()
    test_grid_coarsegrain_rm()
    test_system_coarsegrain()
    test_system_coarsegrain_rm()
    test_ucg_traj_rm()