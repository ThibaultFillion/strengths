import numpy as np
from strengths.typechecking import *
import strengths.value_processing as valproc
from strengths.units import *
from strengths import filepath
from strengths import text_array_rw
import json

class RDGraphSpaceNode : 
    """
    Represents a node from a RDGraphSpace.
    """
    def __init__(self, volume=1, environment=0, units_system=UnitsSystem()) : 
        self.units_system = units_system
        self.volume = volume
        self.environment = environment
 
    @property
    def volume(self) : 
        """
        volume of the node.
        """
        
        return self._volume

    @volume.setter
    def volume(self, v) : 
        self._volume = UnitValue(v, Units(sys=self.units_system, dim=volume_units_dimensions()), convert=False)
    
    @property
    def environment(self) : 
        """
        index of the reaction_diffusion network environment associated with the node/cell.
        """

        return self._environment

    @environment.setter
    def environment(self, v) : 
        self._environment = int(v)

    @property 
    def units_system(self): 
        """
        default units system used when value that require units are given without (:py:class:`UnitsSystem`).
        can be defined from a :py:class:`UnitsSystem` or a :py:class:`dict`.
        ie. 
        
        .. code:: python
        
            rdspace.units_system = UnitsSystem(space="µm", time="s", quantity="molecule")
            rdspace.units_system = {"space"="µm", "time"="s", "quantity"="molecule"}
        """
        
        return self._units_system
     
    @units_system.setter
    def units_system(self, units_system): 
        if isdict(units_system) :
            self._units_system = unitssystem_from_dict(units_system)
        elif type(units_system) == UnitsSystem :
            self._units_system = units_system.copy()
        else :
            raise TypeError("units_system must be a dict or an instance of UnitsSystem.")
            
class RDGraphSpaceEdge : 
    """
    Represents an edge from a RDGraphSpace.
    """
    
    def __init__(self, i, j, surface=1, distance=1, units_system=UnitsSystem()) : 
        self.units_system = units_system
        self._i = int(i)
        self._j = int(j)
        self.surface = surface
        self.distance = distance
 
    @property
    def i(self) : 
        """
        index of the edge's first node (int).
        """
        
        return self._i

    @property
    def j(self) : 
        """
        index of the edge's second node (int).
        """
        
        return self._j      

    @property
    def surface(self) : 
        """
        surface of contact between the linked nodes.
        """
        
        return self._surface

    @surface.setter
    def surface(self, v) : 
        self._surface = UnitValue(v, Units(sys=self.units_system, dim=surface_units_dimensions()), convert=False)
    
    @property
    def distance(self) : 
        """
        distance between the center of the linked nodes.
        """
        
        return self._distance

    @distance.setter
    def distance(self, v) : 
        self._distance = UnitValue(v, Units(sys=self.units_system, dim=space_units_dimensions()), convert=False)

    @property 
    def units_system(self): 
        """
        default units system used when value that require units are given without (:py:class:`UnitsSystem`).
        can be defined from a :py:class:`UnitsSystem` or a :py:class:`dict`.
        ie. 
        
        .. code:: python
        
            rdspace.units_system = UnitsSystem(space="µm", time="s", quantity="molecule")
            rdspace.units_system = {"space"="µm", "time"="s", "quantity"="molecule"}
        """
        
        return self._units_system
     
    @units_system.setter
    def units_system(self, units_system): 
        if isdict(units_system) :
            self._units_system = unitssystem_from_dict(units_system)
        elif type(units_system) == UnitsSystem :
            self._units_system = units_system.copy()
        else :
            raise TypeError("units_system must be a dict or an instance of UnitsSystem.")
            
class RDGraphSpace : 
    """
    represents a graph serving as a reaction-diffusion space. It is a generalization of the cell grid (RDGridSpace)
    as, each cell (node) have its own volume (cell_vol may be an array) and an arbitrary number of neighbors (defined by the edges).
    
    :param nodes: nodes of the graph
    :type nodes: array of RDGraphSpaceNode
        
    :param edges: edges of the graph
    :type edges: array of RDGraphSpaceEdge
    """

    def __init__(
            self, 
            nodes=[],
            edges=[],
            units_system=UnitsSystem()) : 
        
        """
        constructor
        """
        
        self.units_system = units_system
        
        for node in nodes : 
            if type(node) != RDGraphSpaceNode :
                raise TypeError("RDGraphSpace nodes must be RDGraphSpaceNode instances.")
        
        for edge in edges : 
            if type(edge) != RDGraphSpaceEdge :
                raise TypeError("RDGraphSpace edges must be RDGraphSpaceEdge instances.")
                
        self._nodes = tuple(nodes)
        self._edges = tuple(edges)
        
    def check(self) :
        """
        checks for errors in the graph sructures (duplicated edges or invalid node indices).
        """
        
        # chech that there is no duplicated edge (i,j)
        couples = []
        for edge in self.edges : 
            couple = [min(edge.i, edge.j), max(edge.i, edge.j)]
            if couple in couples :
                raise ValueError("duplicated edge (" + str(edge.i) + "," + str(edge.j) + ").")
            else :
                couples.append(couple)
                
        # check that edges node indices are not out of bound
        for edge in self.edges :
            if edge.i<0 or edge.i>=self.size() :
                raise ValueError(str(edge.i) + "in edge (" + str(edge.i) + "," + str(edge.j) + ") is not a valid node index.")
            if edge.j<0 or edge.j>=self.size() :
                raise ValueError(str(edge.j) + "in edge (" + str(edge.i) + "," + str(edge.j) + ") is not a valid node index.")
        
    def size(self) :
        """
        returns the number of cells in the grid (int) (number of nodes).
        """

        return len(self.nodes) 
    
    @property
    def nodes(self) :        
        return self._nodes

    @property
    def edges(self) :
        return self._edges
    
    def get_cell_vol_array(self) :
        """
        Aims to be part of the general reaction-diffusion space interface.
        Returns an array containting the volume of each cell.
        """
        
        return UnitArray([node.volume for node in self.nodes], Units(sys=self.units_system, dim=volume_units_dimensions()))
    
    def get_cell_env_array(self) :
        """
        Aims to be part of the general reaction-diffusion space interface.
        Returns an array containting the environment index associated with each cell.
        """
        
        return [node.environment for node in self.nodes]

    def get_cell_index(self, position) :
        """
        Returns position if it is a valid cell index. raise an Exception otherwise.
        The returned value is int(position) it it is in the [0, self.size()) range.
        """
        
        cell_index = int(position)

        if cell_index<0 or cell_index>=self.size() :
            raise ValueError("position "+str(cell_index)+" is not a valid graph node index.")
        
        return cell_index
    
    def get_edge(self, i, j) : 
        """
        return the edge (i,j) or (j,i) if it exists, None otherwise.
        """
        
        for edge in self.edges : 
            if (edge.i==i and edge.j==j) or (edge.i==j and edge.j==i) :
                return edge
        return None
    
    @property 
    def units_system(self): 
        """
        default units system used when value that require units are given without (:py:class:`UnitsSystem`).
        can be defined from a :py:class:`UnitsSystem` or a :py:class:`dict`.
        ie. 
        
        .. code:: python
        
            rdspace.units_system = UnitsSystem(space="µm", time="s", quantity="molecule")
            rdspace.units_system = {"space"="µm", "time"="s", "quantity"="molecule"}
        """
        
        return self._units_system
     
    @units_system.setter
    def units_system(self, units_system): 
        if isdict(units_system) :
            self._units_system = unitssystem_from_dict(units_system)
        elif type(units_system) == UnitsSystem :
            self._units_system = units_system.copy()
        else :
            raise TypeError("units_system must be a dict or an instance of UnitsSystem.")

    def copy(self) :
        """
        returns a deepcopy of the instance.
        
        .. code:: python
        
            instance.copy()
            
            # is equivalent to
            # import copy
        
            copy.deepcopy(instance)
        
        """
        
        return cpy.deepcopy(self)

def rdgraphspacenode_from_dict(d, parent_units_system = UnitsSystem()) : 
    """
    Builds a RDGridSpaceNode from a dict.
    """    

    d = valproc.process_input_dict_keys(d, [
                ["volume", "vol"],
                ["environment", "env"],
                ["units", "units_system", "units system", "u"]
                
            ]
        )
    
    da = {}

    da["units_system"] = valproc.retrive_units_system_from_dict(
        d = d, 
        default = "inherit",
        parent_units_system = parent_units_system)    
    
    if "volume" in d : 
        da["volume"] = d["volume"]
    
    if "environment" in d : 
        da["environment"] = d["environment"]
    
    return RDGraphSpaceNode(**da)

def rdgraphspacenode_to_dict(node, parent_units_system) : 
    """
    Builds a dict from a RDGridSpaceNode.
    """    
    
    d = { "volume"      : str(node.volume),
          "environment" : node.environment
         }
    
    if node.units_system != parent_units_system :
        d["units"] = unitssystem_to_dict(node.units_system)

    return d

def rdgraphspaceedge_from_dict(d, parent_units_system = UnitsSystem()) : 
    """
    Builds a RDGridSpaceEdge from a dict.
    """    

    d = valproc.process_input_dict_keys(d, [
                ["nodes"],
                ["surface"],
                ["distance"],
                ["units", "units_system", "units system", "u"]
                
            ]
        )
    
    da = {}

    da["units_system"] = valproc.retrive_units_system_from_dict(
        d = d, 
        default = "inherit",
        parent_units_system = parent_units_system)    
    
    if "surface" in d : 
        da["surface"] = d["surface"]
    
    if "distance" in d : 
        da["distance"] = d["distance"]
    
    da["i"] = int(d["nodes"][0])
    da["j"] = int(d["nodes"][1])
    
    return RDGraphSpaceEdge(**da)

def rdgraphspaceedge_to_dict(edge, parent_units_system) : 
    """
    Builds a dict from a RDGridSpaceEdge.
    """    
    
    d = { "nodes"    : [edge.i, edge.j], 
          "surface"  : str(edge.surface),
          "distance" : str(edge.distance)
         }
    
    if edge.units_system != parent_units_system :
        d["units"] = unitssystem_to_dict(node.units_system)

    return d

def rdgraphspace_from_dict(d, parent_units_system = UnitsSystem(), base_path=None) : 
    """
    Builds a RDGridSpace from a dict.
    """    

    d = valproc.process_input_dict_keys(d, [
                ["type"],
                ["nodes"],
                ["edges"],
                ["units", "units_system", "units system", "u"]
                
            ]
        )
    
    da = {}

    da["units_system"] = valproc.retrive_units_system_from_dict(
        d = d, 
        default = "inherit",
        parent_units_system = parent_units_system)    
    
    da["nodes"] = [rdgraphspacenode_from_dict(i, da["units_system"]) for i in d["nodes"]]
    da["edges"] = [rdgraphspaceedge_from_dict(i, da["units_system"]) for i in d["edges"]] 
    
    return RDGraphSpace(**da)

def rdgraphspace_to_dict(space) : 
    """
    Builds a dict from a RDGridSpace.
    """    
    
    d = { "type" : "graph",
          "nodes" : [rdgraphspacenode_to_dict(i, space.units_system) for i in space.nodes],
          "edges" : [rdgraphspaceedge_to_dict(i, space.units_system) for i in space.edges],
          "units" : unitssystem_to_dict(space.units_system)
         }

    return d