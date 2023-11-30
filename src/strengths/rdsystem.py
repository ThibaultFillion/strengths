from strengths.rdspace import RDGridSpace, RDGraphSpace, rdspace_from_dict, rdspace_to_dict, load_rdspace, save_rdspace
from strengths.rdnetwork import RDNetwork, rdnetwork_from_dict, rdnetwork_to_dict, load_rdnetwork, save_rdnetwork
from strengths.units import *
import strengths.value_processing as valproc
from strengths import filepath
from strengths import text_array_rw
import json


def generate_species_state(species, network, space, units_system) :
    """
    Makes the default initial state for a given species. 
    The default state is generated using the species density, the space environment map, and the cell volume : 

    .. math:: 
        
        x_{s,i} = [s]_{env_i} V

    with x the system state, s the species, i the cell position, env_i the environment of cell i, and v the cell volume.
    
    :param species: species for which the state is generated. Must be part of network.
    :type species: Species
    :param network: reaction network the specied is part of.
    :type network: RDNetwork
    :param space: space for which the state must be generated
    :type space: RDSpace
    :param units_system: units system in which the state should be experssed
    :type units_system: UnitsSystem
    :returns: defalut species state in the space.
    :rtype: :py:class:`UnitsSystem`
        
    """

    state = [0 for i in range(space.size())]
    cell_vol = space.get_cell_vol_array()
    cell_env = space.get_cell_env_array()

    for i in range(space.size()) :
        cell_species_density = valproc.get_value_in_env(
            value = species.density, 
            environment = network.environments[cell_env[i]], 
            default = UnitValue(0, "molecule/µm3"))
        
        state[i] = (cell_species_density * cell_vol.get_at(i)).convert(units_system).value
    
    return UnitArray(state, Units(units_system, quantity_units_dimensions()))
        
def generate_system_state(network, space, units_system, state_dict={}) : 
    """
    Makes the initial state for a given reaction network

    :param network: reaction network the specied is part of.
    :type network: RDNetwork
    :param space: space for which the state must be generated
    :type space: RDSpace
    :param units_system: units system in which the state should be experssed
    :type units_system: UnitsSystem
    :param state_dict: dictionnary containing state to apply for the species, overriding the default species states.
        is empty by default.
    :type state_dict: dict
    :returns: full system state for the given reaction diffusion network and space.
    :rtype: :py:class:`UnitsSystem`

    """    

    state = np.array([])

    for s in network.species : 
        if s.label in list(state_dict) : 
            if type(state_dict[s.label]) != UnitArray : 
                raise TypeError("state_dicts's values must be UnitArrays.")
            if state_dict[s.label].units.dim != {"space":0, "time":0, "quantity":1}: 
                raise ValueError("state_dicts's UnitArray dimension must be a quantity.")
            if state_dict[s.label].len() != space.size() : 
                raise ValueError("state_dicts's UnitArrays length must match the syqtem size.")
            state = np.concatenate((state, state_dict[s.label].convert(units_system).value))
        else :
            species_state = generate_species_state(s, network, space, units_system)
            state = np.concatenate((state, species_state.value))

    return UnitArray(state, Units(units_system, {"space":0, "time":0, "quantity":1}))

def generate_species_chemostats(species, network, space) : 
    """
    Makes the default chemostat space map for a given species
    
    :param species: species for which the chemostat map is generated. Must be part of network.
    :type species: Species
    :param network: reaction network the specied is part of.
    :type network: RDNetwork
    :param space: space for which the chemostat map must be generated
    :type space: RDSpace
    :returns: defalut species chemostate map in the space.
    :rtype: array of int
    
    """

    chstt = [0 for i in range(space.size())]
    cell_env = space.get_cell_env_array()
    
    for i in range(space.size()) :
        cell_species_cstt = valproc.get_value_in_env(
            value = species.chstt, 
            environment = network.environments[cell_env[i]], 
            default = 0)
        chstt[i] = int(cell_species_cstt)
            
    return np.array(chstt, dtype=int)

def generate_system_chemostats(network, space, chstt_dict={}) : 
    """
    Makes the chemostat map for a given reaction network
    
    :param network: reaction network the specied is part of.
    :type network: RDNetwork
    :param space: space for which the state must be generated
    :type space: RDSpace
    :param chstt_dict: dictionnary containing chemostat map to apply for the species, overriding the default species states.
        is empty by default.
    :type chstt_dict: dict
    :returns: full system chemostate map for the given reaction diffusion network and space.
    :rtype: array of int

    """
    
    chstt = np.array([], dtype=int)

    for s in network.species : 
        if s.label in list(chstt_dict) :
            if not is_array(chstt_dict[s.label]) : 
                raise TypeError("chstt_dict's values must be arrays")
            if len(chstt_dict[s.label]) != space.size() : 
                raise ValueError("chstt_dict's array size must match the system size.")            
            chstt = np.concatenate((chstt, np.array(chstt_dict[s.label], dtype=int)))            
        else :
            species_chstt = generate_species_chemostats(s, network, space)
            chstt = np.concatenate((chstt, species_chstt))

    return chstt

class RDSystem : 
    """
    Describes a reaction diffusion system, coupling a reaction_network with a space, characterized by a system state, which corresponds
    to the distribution of the chemical species quantities inside the system space, and chemostat maps for each species in the system.
    
    reaction network and space sould no be changed. If a change of system size or number of species is required, creation of a new RDSystem is highly recomended.
    
    
    :param network: system reaction_diffusion network
    :type network: RDNetwork
    :param space: system space (default RDSpace())
    :type space: RDSpace
    :param state: system state to be used. if set to None or with a dict, the a default system state is generated from network and space (see generate_system_state) (default None). 
    :type state: array, UnitArray, None or dict    
    :param chemostats: system chemostat map to be used. if set to None or with a dict, the a default system chemostat map is generated from network and space (see generate_system_chemostats) (default None). 
    :type chemostats: array, None or dict
    :param units_system: default units system (default UnitsSystem())
    :type units_system: UnitsSystem
    """

    def __init__(self, network, space=RDGridSpace(), state=None, chemostats=None, units_system=UnitsSystem()) :
        """
        constructor
        """

        self.units_system = units_system.copy()

        self.network = network
        self.space = space

        if isnone(state) : 
            self.set_default_state()
        elif isdict(state) :
            self.set_default_state(state)
        else :
            self.state = state
        
        if isnone(chemostats) : 
            self.set_default_chemostats()
        elif isdict(chemostats) :
            self.set_default_chemostats(chemostats)
        else :
            self.chemostats = chemostats        

    @property
    def network(self) :
        """
        system reaction-diffusion network (:py:class:`RDNetwork`)
        """
        
        return self._network
        
    @network.setter
    def network(self, v) :
        if type(v) != RDNetwork :
            raise ValueError("network must be a :py:class:`RDNetwork`.")
        self._network = v
    
    @property
    def space(self) :
        """
        system space (:py:class:`RDSpace`).
        """
        
        return self._space

    @space.setter
    def space(self, v) :
        if type(v) not in [RDGridSpace, RDGraphSpace] :
            raise ValueError("space must be a :py:class:`RDGridSpace` or :py:class:`RDGraphSpace`.")
        self._space = v
    
    @property
    def state(self) :
        """
        State of the system (:py:class:`UnitArray`). It corresponds to the quantity of each species in each cell of the system.
        It can be set with an array or a UnitArray with quantity units dimensions.
        ie.
        
        .. code:: python
        
            # considering rdsystem.state_size() == 8
            
            rdsystem.state = [1,1,1,1,1,1,1,1]
            rdsystem.state = UnitArray([1,1,1,1,1,1,1,1], "molecule")
        
        """
        
        return self._state
        
    @state.setter
    def state(self, v) :
        if isarray(v) :
            self._state = v = UnitArray(v, Units(self.units_system, quantity_units_dimensions()))
        elif type(v) == UnitArray : 
            
            if v.units.dim != quantity_units_dimensions() : 
                raise ValueError("state must be in quantity units.")

            self._state = v.copy()
        else : 
            raise TypeError("invalid type for state. must be a UnitArray or an array.")
        
    @property
    def chemostats(self) :
        """
        Map of the chemostats for each species in the system (array of int).
        Can be set by a single value or an array.
        If set as an array, the array must match state_size.
        
        """
        
        return self._chemostats
        
    @chemostats.setter
    def chemostats(self, v) :        
        if not isarray(v) : 
            raise ValueError("chemostat map must be an array, a dict or None")
        v = np.array(v, dtype=int)
                
        self._chemostats = v

    @property 
    def units_system(self): 
        """
        Default units system used when value that require units are given without (:py:class:`UnitsSystem`).
        can be defined from a :py:class:`UnitsSystem` or a :py:class:`dict`.
        ie. 
        
        .. code:: python
        
            rdsystem.units_system = UnitsSystem(space="µm", time="s", quantity="molecule")
            rdsystem.units_system = {"space"="µm", "time"="s", "quantity"="molecule"}
            
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

    def set_default_state(self, override_species_state_dict={}) : 
        """
        sets the default system state (see generate_system_chemostats).
        
        :param override_species_state_dict: optionnal. species states that should override the default state.
        :type override_species_state_dict: dict
        
        """

        self._state = generate_system_state(self.network, self.space, self.network.units_system, override_species_state_dict)    
    
    def set_default_chemostats(self, override_species_chstt_dict={}) : 
        """
        Sets the default chemostat map (see generate_system_chemostats).

        :param override_species_chstt_dict: optionnal. species chemostat maps that should override the default state.
        :type override_species_chstt_dict: dict
        
        """

        self._chemostats = generate_system_chemostats(self.network, self.space, override_species_chstt_dict)

    def reset_state(self) :
        """
        Sets all species quantity to 0 in all cells.
        
        """

        self._state = np.zeros(self.space.size()*self.network.nspecies())

    def reset_chemostats(self) :
        """
        Sets all chemostat map values to False/0.
        
        """

        self._chemostats = np.zeros(self.space.size()*self.network.nspecies(), dtype=int)

    def get_cell_index(self, position) :
        """
        Returns the linear index of a cell from its 3D coordinates (x, y, z).
        
        :param position: cell index or coordinates
        :type position: int, Coord like or tuple
        :returns: index of the species at the given position.
        :rtype: int
        
        ie.
        
        .. code:: python
        
            rdsystem.get_cell_index((1,0,0))
            rdsystem.get_cell_index(Coord(1,0,0))
            rdsystem.get_cell_index(1)
        
        """

        return self.space.get_cell_index(position)

    def get_state_index(self, species, position) :
        """
        Returns the linear index of the quantity of a given species in a given cell in the system state array from its 3D coordinates (x, y, z).

        :param species: species of interest. it must be part of the reaction network.
            if given as a str, it must be a species label.
            if given as a int, it must be a species index.
            
        :type speces: Species, str or int
        :param position: cell index or coordinates
        :type position: int or Coord like
        :returns: index for the species state at the given position.
        :rtype: int
        
        ie.
        
        .. code:: python
        
            rdsystem.get_state_index(1, (1,0,0))
            rdsystem.get_state_index("A", Coord(1,0,0))
            rdsystem.get_state_index(rdsystem.network.get_species("A"), 1)

        """

        species_index = self.network.get_species_index(species)
        cell_index    = self.space.get_cell_index(position)

        return species_index * self.space.size() + cell_index
    
    def set_chemostat(self, species, position, value) : 
        """
        Sets the chemostat map value for a given species at a given position.
        
        :param species: species for which the value must be set
        :type species: Species, int or str
        :param pos: position at which the value must be set
        :type pos: number (interpreted as a linear index (int)) or a class with x, y and z properties/members.
        :param value: value to be set
        :type value: int or bool
        """
        
        state_index = self.get_state_index(species, position)
                
        self._chemostats[state_index] = int(value)

    def set_state(self, species, position, value) :
        """
        Sets the chemostat map value for a given species at a given position.
        
        :param species: species for which the value must be set
        :type species: Species, int or str
        :param pos: position at which the value must be set
        :type pos: number (interpreted as a linear index (int)) or a class with x, y and z properties/members.
        :param value: value to be set
        :type value: number or UnitValue with quantity units dimensions.
        """
        
        state_index = self.get_state_index(species, position)

        if type(value) != UnitValue :
            value = UnitValue(value, Units(self.units_system, quantity_units_dimensions()))

        self._state.set_at(state_index, value)

    def state_size(self) : 
        """
        Returns the size of the system state (int).

        size = w*h*d*nspecies
        """
        
        return self.space.size() * self.network.nspecies()

    def is_valid(self) : 
        """
        Checks is the system is valid. Error messages are stored internally, and can be retrived by calling the get_error method.
        
        :rtype: bool
        :returns: True : no problem have been found. False : the system is invalid.
        """
        
        valid = True
        self._error_msg = ""
        
        if self.space.size() * self.system.nspecies() != len(self.state) : 
            valid = False
            self._error_msg += "error : state size does not match space size * number of species."
        
        if self.space.size() * self.system.nspecies() != len(self.chemostats) : 
            valid = False
            self._error_msg += "error : chemostat map size does not match space size * number of species."
            
    def get_error(self) : 
        """
        Returns the error messages built during the previous call to the is_valid method.
        
        :rtype: str
        """
                
        return self._error_msg

    def make_dxdtf(self, units_system=UnitsSystem()) :
        """
        Returns a function dxdt = f(t, x) for systems with a system of size 1.
        """
        
        if(self.space.size() != 1) : 
            raise NotImplementedError("make_dxdtf is only implemented for systems of size 1.")
        
        N = len(self.network.species)
        M = len(2*self.network.reactions)
        
        env = self.network.environments[self.space.get_cell_env_array()[0]]
        vol = self.space.get_cell_vol_array().get_at(0).convert(units_system).value
        sl  = self.network.species_labels()
    
        reactions = []
    
        for r in self.network.reactions : 
            r1, r2 = r.split()
            reactions.append(r1)
            reactions.append(r2)
        
        #rates that already account for the volume
        k = []
        for r in reactions :
            k_r = 0
            if (r.environments is None) or (env in r.environments) :
                k_r = r.kf.convert(units_system).value
            k_r *= vol**(1-r.order())
            k.append(k_r)
        
        #stoechiometry
        sub = [list(r.ssto(sl)) for r in reactions]
        sto = [list(r.dsto(sl)) for r in reactions]
        
        chemostats = [1-self.chemostats[i] for i in range(N)]
        
        def dxdtf(t, x) :         
            rates = [k[i] for i in range(M)]
            dxdt  = [0    for i in range(N)]
    
            for r in range(M) : 
                for s in range(N) :
                    rates[r] *= x[s]**sub[r][s]
    
            for s in range(N) :
                for r in range(M) : 
                    dxdt[s] += rates[r]*sto[r][s]
                dxdt[s] *= chemostats[s]
    
            return dxdt
    
        return dxdtf
        
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
    
def rdsystem_from_dict(d, parent_units_system=UnitsSystem(), base_path=None):
    """
    Returns a RDSystem created acoording to the dictionnary d.

    :param d: dictionary representing the reaction diffusion system
    :type d: dict
    :returns: reaction diffusion system created from d.
    """
    
    d = valproc.process_input_dict_keys(d, [
                ["network", "rdnetwork"],
                ["space", "rdspace"],
                ["state"],
                ["chemostats"],
                ["units", "units_system", "units system", "u"]
            ]
        )    
    
    da = {}

    da["units_system"] = valproc.retrive_units_system_from_dict(
        d = d, 
        default = "inherit",
        parent_units_system = parent_units_system
        )
    
    if "network" in d : 
        rdn = d["network"]
    
        if isdict(rdn) :
            rdn = rdnetwork_from_dict(rdn, da["units_system"], base_path)
        elif isstr(rdn) :
            rdn = filepath.get_path_with_base(rdn, base_path)
            rdn = load_rdnetwork(rdn, da["units_system"])
        else :
            raise ValueError("network type is unexpected.")
        da["network"] = rdn
    else :
        raise ValueError("missing system network.")

    if "space" in d :
        space = d["space"]
    
        if isdict(space) :
            space = rdspace_from_dict(space, da["units_system"], base_path)
        elif isstr(space) :
            space = filepath.get_path_with_base(space, base_path)
            space = load_rdspace(space, da["units_system"])
        else :
            raise ValueError("space type is unexpected.")
        da["space"] = space        
        
    if "state" in d :
        state = d["state"]
        
        if isdict(state) : 
            state = unitarray_from_dict(state, base_path=base_path)
        da["state"] = state
    
    if "chemostats" in d : 
        chemostats = d["chemostats"]
        
        if isstr(chemostats) : 
            chemostats = filepath.get_path_with_base(chemostats, base_path)
            if filepath.have_extension(chemostats, ".npy") :
                chemostats = np.load(chemostats)
            else :
                chemostats = text_array_rw.load_1D_array_txt(chemostats, int)
        da["chemostats"] = chemostats

    return RDSystem(**da)

def rdsystem_to_dict(rds):
    """
    Returns a dictionnary describing the reaction network.
    
    :param network: reaction network to be converted
    :type network: RDNetwork
    :returns: dictionary representing network.
    """

    d = {
        "units" : unitssystem_to_dict(rds.units_system),
        "network"   : rdnetwork_to_dict(rds.network), 
        "space" : rdspace_to_dict(rds.space),
        "state" : unitarray_to_dict(rds.state),
        "chemostats" : array_to_list(rds.chemostats)
        }

    return d

def load_rdsystem(path, parent_units_system=UnitsSystem()):
    """
    Returns a RDSystem created acoording to the dictionnary d loaded from the JSON file at the given path.
    
    :param path: relative or absolute path to a json file.
    :type path: str
    :returns: reaction diffusion system built according to the JSON file.
    """

    f = open(path, "r", encoding="utf-8")
    d = json.load(f)
    f.close()
    return rdsystem_from_dict(d, parent_units_system=parent_units_system, base_path=filepath.get_base_path(path))

def save_rdsystem(rds, path):
    """
    Saves the dictionnary desctiption of the reaction diffusion system **rds** as a JSON file at the given **path**.

    :param rds: reaction diffusion system to be saved as a JSON file.
    :type rds: RDSystem
    :param path: JSON save file relative or absolute path. 
        a JSON is created at this location, or replaced if it already exists
    :type path: str
    :returns: None
    """

    d = rdsystem_to_dict(rds)
    f = open(path, "w", encoding="utf-8")
    json.dump(d, f, indent = 4)    
    f.close()
