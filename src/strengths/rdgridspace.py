import numpy as np
from strengths.typechecking import *
import strengths.value_processing as valproc
from strengths.units import *
from strengths import filepath
from strengths import text_array_rw
import json
import copy

"""
Module that contains the implementation of the Reaction Diffusion Space classes, and related functions.
"""

class RDGridSpace : 
    """
    represent a grid of cells in 3 dimensions. each cell is associated with an environment (index).
    all cells have the same volume cell_vol.
    
    :param w: grid width
    :type w: int
    :param h: grid height
    :type h: int
    :param d: grid depth
    :type d: int
    :param cell_env: cell(s) environment(s)
    :param cell_vol: cell(s) volume(s)   
    """

    def __init__(
            self, 
            w=1, 
            h=1, 
            d=1, 
            cell_env=0, 
            cell_vol=1, 
            boundary_conditions=None, 
            units_system=UnitsSystem()) : 

        """
        constructor
        """
        
        self.units_system = units_system.copy()

        self._w = int(w)
        self._h = int(h)
        self._d = int(d)
        
        if self._w <= 0 : 
            raise ValueError("system width must a strictly positive integer.")
        if self._h <= 0 : 
            raise ValueError("system height must a strictly positive integer.")        
        if self._d <= 0 : 
            raise ValueError("system depth must a strictly positive integer.")

        self.cell_env = cell_env
        self.cell_vol = cell_vol
        
        if isnone(boundary_conditions):
            boundary_conditions = {}
            
        self.set_boundary_conditions(boundary_conditions)
    
    @property
    def w(self) :         
        """
        grid width in cells (int). Cannot be set.
        """
        
        return self._w
    
    @property
    def h(self) :         
        """
        grid height in cells (int). Cannot be set.
        """

        return self._h

    @property
    def d(self) :         
        """
        grid depth in cells (int). Cannot be set.
        """

        return self._d
    
    def size(self) :
        """
        returns the number of cells in the grid (int).
        size = w*h*d
        """

        return self.w*self.h*self.d 

    def get_cell_coordinates(self, cell_index, return_type=tuple) :
        """
        Returns the 3D (x, y, z) coordinates of the cell with the given index.
        
        :param cell_index: cell index
        :type cell_index: number
        :param return_type: type of the returned coordinates. default tuple.
        :type return_type: class/type        
        :returns: index of the species at the given position.
        :rtype: int
        
        ie.
        
        .. code:: python
        
            rdspace.get_cell_index(1) # returns (1,0,0)
            
        """
        
        if not self.is_within_bounds(cell_index) :
            raise ValueError("index "+str(cell_index)+" is not within the grid bounds.")
        
        cell_index = int(cell_index)
        x = int(cell_index%self.w)
        y = int((cell_index%(self.w*self.h))/self.w)
        z = int(cell_index/(self.w*self.h)) 

        if return_type == tuple :
            return (x, y, z)
        else :
            r = return_type()
            r.x = x
            r.y = y
            r.z = z
            return r
        
    def get_cell_index(self, position) :
        """
        Returns the linear index at a given position.
        
        :param position: cell index or coordinates
        :type position: number, Coord like or tuple
            If a number is given, it is interpreted as a linear index (the function returns int(index)).
            If a tuple or a Coord like object is given, those are interperted as (x,y,z) coordinates, and the index 
            corresponding to those coordinates is returned.
        :returns: index of the species at the given position.
        :rtype: int
        
        ie.
        
        .. code:: python
        
            rdspace.get_cell_index((1,0,0))
            rdspace.get_cell_index(Coord(1,0,0))
            rdspace.get_cell_index(1)
            
        """
        if not self.is_within_bounds(position) :
            raise ValueError("position "+str(position)+" is not within the grid bounds.")
            
        if isnumber(position) : 
            return int(position)
        elif isarray(position) :
            return int(position[0]) + int(position[1])*self.w + int(position[2])*self.w*self.h
        else :
            return int(position.x)  + int(position.y)*self.w  + int(position.z)*self.w*self.h

    @property
    def cell_vol(self) :        
        """
        volume of a cell (:py:class:`UnitValue`).
        can be set from a number, a :py:class:`UnitValue` or a str.
        ie.
        
        .. code:: python
        
            rdspace.cell_vol = 1
            rdspace.cell_vol = "1 µm3"
            rdspace.cell_vol = UnitValue(1, "µm3")
        
        """
        
        return self._cell_vol

    @cell_vol.setter        
    def cell_vol(self, v) :
        self._cell_vol = valproc.process_unitvar_input(
            v, 
            self.units_system, 
            volume_units_dimensions(),
            accepts_singlevalue = True,
            accepts_dict = False,
            accepts_array = False)

    @property
    def cell_env(self) :    
        """
        environment index associated with each cell (array of int).
        can be set from a number or an array of numbers.
        if set with number, this number will be applied to all cells.
        if set with an array, the array size must match the rdspace size.
        ie. 
        
        .. code:: python
            
            rdspace.cell_env = 0
            rdspace.cell_env = [0,0,0]
        """

        if isnumber(self._cell_env) :
            return np.array([self._cell_env for i in range(self.size())], dtype=int)
        elif isarray(self._cell_env) :
            return self._cell_env

    @cell_env.setter        
    def cell_env(self, v) :
        if isnumber(v) :
            self._cell_env = int(v)
        elif isarray(v) :
            if len(v) != self.size() : 
                raise ValueError("cell_env size must match the system size.")
            self._cell_env = np.array([int(vi) for vi in v], dtype=int)
        else : 
            raise TypeError("cell_env must a number or an array.")
    
    def get_cell_vol_array(self) :
        """
        Aims to be part of the general reaction-diffusion space interface.
        Returns an array containting the volume of each cell.
        """
        
        return UnitArray([self.cell_vol.value for i in range(self.size())], self.cell_vol.units).convert(self.units_system)
    
    def get_cell_env_array(self) :
        """
        Aims to be part of the general reaction-diffusion space interface.
        Returns an array containting the environment index associated with each cell.
        """
        
        return copy.deepcopy(self.cell_env)
    
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

    def is_within_bounds(self, position) :
        """
        Tells if a given position is within the bounds of the cell grid.
        returns True if the position refers to a cell inside the grid, False otherwise.
        
        :param position: position to be evaluated
        :type position: number, tuple or Coord like
        :returns: True if the position refers to a cell inside the grid, False otherwise.
        :rtype: bool
        """
        
        if isnumber(position) :
            return (int(position)>=0 or int(position)<self.size())
        elif isarray(position) : 
            return (int(position[0])>=0 and
                    int(position[0])<self.w and
                    
                    int(position[1])>=0 and
                    int(position[1])<self.h and
                    
                    int(position[2])>=0 and
                    int(position[2])<self.d)
        else :
            return (int(position.x)>=0 and
                    int(position.x)<self.w and
                    
                    int(position.y)>=0 and
                    int(position.y)<self.h and
                    
                    int(position.z)>=0 and
                    int(position.z)<self.d)    

    def are_neighbors(self, position1, position2): 
        """
        tells if the two positions reffres to a couple of neighbor cells.
        Return True is that is the case, False otherwise.
        
        :param position1: position of the first cell
        :type position1: number, tuple or Coord like
        :param position2: position of the second cell
        :type position2: number, tuple or Coord like
        """
        
        if not self.is_within_bounds(position1) :
            raise ValueError("position1 "+str(position1)+" is not within the grid bounds.")
            
        if not self.is_within_bounds(position2) :
            raise ValueError("position2 "+str(position2)+" is not within the grid bounds.")
            
        coord1 = self.get_cell_coordinates(self.get_cell_index(position1))
        coord2 = self.get_cell_coordinates(self.get_cell_index(position2))
        
        dx = abs(coord1[0] - coord2[0])
        dy = abs(coord1[1] - coord2[1])
        dz = abs(coord1[2] - coord2[2])
        
        if self._boundary_conditions["x"] == "periodical" :
            dx = min(dx, abs(self.w-dx))
        if self._boundary_conditions["y"] == "periodical" :
            dy = min(dy, abs(self.h-dy))
        if self._boundary_conditions["z"] == "periodical" :
            dz = min(dz, abs(self.d-dz))
            
        return ((dx + dy + dz) == 1)

    def get_boundary_conditions(self) :
        """
        Returns a dict of boundary conditions associated with each axis.
        """
        
        return self._boundary_conditions.copy()
        
    def set_boundary_conditions(self, boundary_conditions) :
        """
        Sets the boundary conditons of the system space.
        """
        self._boundary_conditions = {
            "x" : "reflecting", 
            "y" : "reflecting", 
            "z" : "reflecting"}
        
        if not isdict(boundary_conditions) :
            raise TypeError("boundary_conditions must be a dict.")

        for axis in list(boundary_conditions) :
            if axis not in ["x", "y", "z"] :
                raise ValueError("axis must be \"x\", \"y\" or \"z\".")
            if boundary_conditions[axis] not in ["reflecting", "periodical"] :
                raise ValueError("conditon must be \"reflecting\" or \"periodical\".")
            self._boundary_conditions[axis] = boundary_conditions[axis]
        
def rdgridspace_from_dict(d, parent_units_system = UnitsSystem(), base_path=None) : 
    """
    Builds a RDGridSpace from a dict.
    """    

    d = valproc.process_input_dict_keys(d, [
                ["type"],
                ["w", "width"],
                ["h", "height"],
                ["d", "depth"],
                ["cell_env", "cell_environments", "cell environments", "environments", "env"],
                ["cell_volume", "cell_vol"],
                ["boundary_conditions"],
                ["units", "units_system", "units system", "u"]
                
            ]
        )
    
    da = {}
    
    if "w" in d : da["w"] = d["w"]
    if "h" in d : da["h"] = d["h"]
    if "d" in d : da["d"] = d["d"]
    if "cell_env" in d :
        cell_env = d["cell_env"]
        
        if isstr(cell_env) : 
            cell_env = filepath.get_path_with_base(cell_env, base_path)
            if filepath.have_extension(cell_env, ".npy") :
                cell_env = np.load(cell_env)
            else :
                cell_env = text_array_rw.load_1D_array_txt(cell_env, int)
        da["cell_env"] = cell_env
    if "cell_volume"         in d : da["cell_vol"] = d["cell_volume"]
    if "boundary_conditions" in d : da["boundary_conditions"] = d["boundary_conditions"]
    da["units_system"] = valproc.retrive_units_system_from_dict(
        d = d, 
        default = "inherit",
        parent_units_system = parent_units_system)    
    
    return RDGridSpace(**da)

def rdgridspace_to_dict(mg) : 
    """
    Builds a dict from a RDGridSpace.
    """    
    
    d = { "type" : "grid",
          "units" : unitssystem_to_dict(mg.units_system),
          "w" : mg.w,
          "h" : mg.h,
          "d" : mg.d,
          "cell_env" : array_to_list(mg.cell_env),
          "cell_volume" : valproc.format_unitvar_for_save(mg.cell_vol, mg.units_system),
          "boundary_conditions" : mg.get_boundary_conditions()
          }

    return d