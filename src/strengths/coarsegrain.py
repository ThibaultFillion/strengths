from strengths.units import *
from strengths.rdsystem import RDSystem
from strengths.rdengine import RDEngineBase
from strengths.rdoutput import  RDTrajectory
from strengths.rdspace import RDGridSpace, RDGraphSpace, RDGraphSpaceNode, RDGraphSpaceEdge
from strengths.typechecking import *

def grid_to_graph(grid) : 
    """
    Returns a RDGraphSpace equivalent to the input RDGridSpace grid.
    """
        
    nodes = []
    edges = []
    
    edge_dst = (grid.cell_vol)**(1/3)
    edge_sfc = edge_dst**2
    
    for i in range(grid.size()) :
        nodes.append(RDGraphSpaceNode(volume=grid.cell_vol.copy(), environment=grid.cell_env[i], units_system=grid.units_system))
    
    for z in range(grid.d) :
        for y in range(grid.h) :
            for x in range(grid.w) :
                if x<grid.w-1 : edges.append(RDGraphSpaceEdge(i=grid.get_cell_index((x,y,z)), j=grid.get_cell_index((x+1,y,z)), surface=edge_sfc, distance=edge_dst, units_system=grid.units_system))
                if y<grid.h-1 : edges.append(RDGraphSpaceEdge(i=grid.get_cell_index((x,y,z)), j=grid.get_cell_index((x,y+1,z)), surface=edge_sfc, distance=edge_dst, units_system=grid.units_system))
                if z<grid.d-1 : edges.append(RDGraphSpaceEdge(i=grid.get_cell_index((x,y,z)), j=grid.get_cell_index((x,y,z+1)), surface=edge_sfc, distance=edge_dst, units_system=grid.units_system))

    if grid.get_boundary_conditions()["x"] == "periodical" :
        for z in range(grid.d) :
            for y in range(grid.h) :
                edges.append(RDGraphSpaceEdge(i=grid.get_cell_index((grid.w-1,y,z)), j=grid.get_cell_index((0,y,z)), surface=edge_sfc, distance=edge_dst, units_system=grid.units_system))

    if grid.get_boundary_conditions()["y"] == "periodical" :
        for z in range(grid.d) :
            for x in range(grid.w) :
                edges.append(RDGraphSpaceEdge(i=grid.get_cell_index((x,grid.h-1,z)), j=grid.get_cell_index((x,0,z)), surface=edge_sfc, distance=edge_dst, units_system=grid.units_system))

    if grid.get_boundary_conditions()["z"] == "periodical" :
        for y in range(grid.h) :
            for x in range(grid.w) :
                edges.append(RDGraphSpaceEdge(i=grid.get_cell_index((x,y,grid.d-1)), j=grid.get_cell_index((x,y,0)), surface=edge_sfc, distance=edge_dst, units_system=grid.units_system))
                
    graph = RDGraphSpace(nodes=nodes, edges=edges, units_system=grid.units_system)
    return graph

def check_index_map_validity(im, space) :
    """
    check if the index map im associated with space is valid. Raises an exception if it is not.
    A valid index map contains only positive integers, and must contain every integer
    between 0 and max(im). The index map size must mat that of space, and all input cells with a same
    coarse grained index should also share the same environment index.
    """
    
    #check length
    if len(im) != space.size() :
        raise ValueError("index map size must match the corresponding space size.")
    
    #check type
    for i in im :
        if type(i) != int :
            raise TypeError("Coarse graining index map should only contain integers.")
            
    im_max = max(im)
    im_min = min(im)
    
    #check min
    if im_min<-1 :
        raise ValueError("Coarse graining index map should not contain negative values except -1.")
    
    #check max
    if im_max<0 : 
        raise ValueError("Coarse graining to an empty graph in invalid.")
    
    #check the presence of every positive integer between 0 and im_max
    for i in range(0, im_max) :
        if i not in im :
            raise ValueError("Missing index "+str(i)+"in coarse graining index map.")

    #check that coarse graining doesnt mix different environments
    env = space.get_cell_env_array()
    env_out = [-2 for i in range(min(im), max(im)+1)]
    for i in range(space.size()) :
        if env_out[im[i]] == -2 :
            env_out[im[i]] = env[i]
        elif env_out[im[i]] == env[i] : 
            pass
        else :
            raise ValueError("output node "+str(im[i])+" contains different environments.")
    
    
def coarsegrain_grid(grid, index_map) :
    """
    grid is the RDGirdSpace to be coarse_grained as a smaller graph space. grid must not have periodical boundary conditions.
    index_map is an array with the same size as the input space, which associated to every cell of the input graph an index of the coarse grained output graph space.
    It returns a RDGraph. index_map additionnaly support the negative index -1. Cells indexed with -1 are excluded from the output graph.
    """

    #validity check
    if grid.get_boundary_conditions() != {"x" : "reflecting", "y" : "reflecting", "z" : "reflecting"} :
        raise ValueError("conversion of graphs with non reflecting boundary conditions are not supported yet.")
        
    grid_cell_edge = grid.cell_vol**(1/3)

    #store input gird cell positions
    in_node_pos = []
    for z in range(grid.d) :
        for y in range(grid.h) :
            for x in range(grid.w) :
                in_node_pos.append([x*grid_cell_edge,y*grid_cell_edge,z*grid_cell_edge])
            
    space = grid_to_graph(grid)
    
    #validity check
    check_index_map_validity(index_map, space)

    n_cell_out = max(index_map)+1 #number of nodes in the output graph
    
    nodes = [RDGraphSpaceNode(volume=0, environment=0, units_system=space.units_system) for i in range(n_cell_out)]
    node_pos = [[0,0,0] for i in range(n_cell_out)]
    node_ncg = [0 for i in range(n_cell_out)]
    
    # loop that sets the nodes informations
    for i in range(space.size()):
        if index_map[i] != -1 : 
            nodes[index_map[i]].volume += space.nodes[i].volume
            nodes[index_map[i]].environment = space.nodes[i].environment
            node_pos[index_map[i]][0] += in_node_pos[i][0]
            node_pos[index_map[i]][1] += in_node_pos[i][1]
            node_pos[index_map[i]][2] += in_node_pos[i][2]
            node_ncg[index_map[i]] += 1
    
    #completes the averaging of nodes center positions
    for i in range(n_cell_out):
        node_pos[i][0] /= node_ncg[i]
        node_pos[i][1] /= node_ncg[i]
        node_pos[i][2] /= node_ncg[i]
        
    edges = []
    out_edge_coords = []

    #adding the edges and computing their surfaces
    for edge in space.edges :
        i = index_map[edge.i]
        j = index_map[edge.j]
        c = (min(i,j), max(i,j))

        if i == j : 
            continue
        
        if i == -1 or j == -1 : 
            continue
        
        #if (i, j) have already been added to the output graph
        if c in out_edge_coords :
            for out_edge in edges :
                if out_edge.i==c[0] and out_edge.j == c[1] :
                    out_edge.surface += edge.surface
                    break
        else :
            edges.append(RDGraphSpaceEdge(i=c[0], j=c[1], surface=edge.surface, distance=0, units_system=space.units_system))
            out_edge_coords.append(c)

    #sets the edges distances
    for edge in edges :
        distance = (
        (node_pos[edge.i][0] - node_pos[edge.j][0])**2 +
        (node_pos[edge.i][1] - node_pos[edge.j][1])**2 +
        (node_pos[edge.i][2] - node_pos[edge.j][2])**2 )**(1/2)
        edge.distance = distance
    
    return RDGraphSpace(nodes=nodes, edges=edges, units_system=space.units_system)

def coarsegrain_system(system, index_map) :
    """
    Returns a coarse grained version of system. system must be a RDSystem with a RDGridSpace space,
    and index_map must be a valid index map.
    """
    
    cgspace = coarsegrain_grid(system.space, index_map)
    nspecies = system.network.nspecies()
    
    cgstate = [0 for i in range(cgspace.size()*nspecies)]
    cgchstt = [0 for i in range(cgspace.size()*nspecies)]
    
    #iterating over the input space cells, to aggregate chemostats and quantities
    for i in range(system.space.size()) :
        if index_map[i] != -1 :
            for s in range(nspecies) :
                cgstate[s * cgspace.size() + index_map[i]] += system.state.value[s * system.space.size() + i]
                cgchstt[s * cgspace.size() + index_map[i]] += system.chemostats [s * system.space.size() + i]
    
    for i in range(len(cgchstt)) :
        cgchstt[i] = int(min(cgchstt[i], 1)) 
    
    cgstate = UnitArray(np.array(cgstate), system.state.units, check_value=False)
    
    cgsystem = system.copy()
    cgsystem.space = cgspace
    cgsystem.state = cgstate
    cgsystem.chemostats = cgchstt
    
    return cgsystem

def uncoarsegrain_trajectory_data(trajectory, ncg_space, index_map) :
    """
    Return an uncoarsegrained version of a trajectory data.
    """
    
    cg_space = trajectory.system.space
    cg_nodes = [[] for i in range(cg_space.size())]
    
    for i in range(ncg_space.size()) : 
        if index_map[i] != -1 :
            cg_nodes[index_map[i]].append(i)
    
    state_size = trajectory.system.network.nspecies() * ncg_space.size()
    data = np.zeros(state_size*trajectory.nsamples())
    
    in_state = trajectory.data.value.reshape((trajectory.nsamples(), trajectory.system.network.nspecies(), cg_space.size()))
    for n in range(trajectory.nsamples()) :
        for s in range(trajectory.system.network.nspecies()) :
            for node_index in range(len(cg_nodes)) :
                for j in cg_nodes[node_index] :
                    data[n*state_size + s*ncg_space.size() + j] = in_state[n, s, node_index]/len(cg_nodes[node_index])

    return UnitArray(data, trajectory.data.units, check_value=False)

def uncoarsegrain_trajectory(trajectory, ncg_system, index_map) :
    """
    Return an uncoarsegrained version of a trajectory.
    """
    
    return RDTrajectory(
        data = uncoarsegrain_trajectory_data(trajectory, ncg_system.space, index_map),
        t_sample = trajectory.t,
        system = ncg_system,
        script = trajectory.script,
        engine_description = trajectory.engine_description,
        engine_option = trajectory.engine_option,
        cgmap = index_map
        )
    
    pass