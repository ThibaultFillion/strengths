import numpy as np
import strengths.value_processing as valproc
from strengths.units import UnitsSystem
from strengths import filepath
from strengths.rdgridspace import RDGridSpace, rdgridspace_from_dict, rdgridspace_to_dict
from strengths.rdgraphspace import RDGraphSpace, RDGraphSpaceNode, RDGraphSpaceEdge, rdgraphspace_from_dict, rdgraphspace_to_dict
import json

"""
Module that contains the implementation of the Reaction Diffusion Space classes, and related functions.
"""

def rdspace_from_dict(d, parent_units_system = UnitsSystem(), base_path=None) : 
    """
    Builds a Reaction Diffudion Space (RDGridSpace or RDGraphSpace) from a dict.
    """        
    
    if "type" not in d : 
        d["type"] = "grid"
    
    if d["type"] == "grid" : 
        return rdgridspace_from_dict(d, 
                                     parent_units_system=parent_units_system, 
                                     base_path=base_path)
    elif d["type"] == "graph" :
        return rdgraphspace_from_dict(d, 
                                      parent_units_system=parent_units_system, 
                                      base_path=base_path)
    else :
        raise ValueError("unspported space type.")

def rdspace_to_dict(space) : 
    """
    Builds a dict from a RDSpace.
    """    
    
    if type(space) == RDGridSpace :
        return rdgridspace_to_dict(space)
    elif type(space) == RDGraphSpace : 
        return rdgraphspace_to_dict(space)
    else :
        raise TypeError("unsupported space type.")

def load_rdspace(path, parent_units_system = UnitsSystem()):
    """
    Loads an RDSpace object from a JSON file.
    """

    f = open(path, "r", encoding="utf-8")
    d = json.load(f)
    return rdspace_from_dict(d, parent_units_system, base_path=filepath.get_base_path(path))

def save_rdspace(space, path):
    """
    Saves an RDSpace object as a JSON file.
    """

    d = rdspace_to_dict(space)
    f = open(path, "w", encoding="utf-8")
    json.dump(d, f, indent = 4)    
    f.close()
