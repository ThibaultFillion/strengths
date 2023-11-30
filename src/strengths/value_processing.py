from strengths.units import *
from strengths.typechecking import *

import copy
import string
import strengths.filepath


# processing of variables.


def is_unitarray_dict(d) :
    """
    returns true id d is a dict representing akin to represent a unitarray (a dict that possesses a "units"
    or a "value" key).
    """
    
    if not isdict(d) :
        return False
    
    return ((d.get("units", None) != None) or (d.get("value", None) != None))
    
    
def process_unitvar_input (
        v, 
        units_system, 
        units_dimensions,
        accepts_singlevalue,
        accepts_dict,
        accepts_array
        ) :
    """
    processing of function arguments expecting UnitValues, UnitArrays of dict of UnitValues/UnitArrays.

    Return a value equivalent to v, where numbers/arrays have been replaced by UnitValue/UnitArrays.
    Preexisting UnitValue/UnitArrays are checked to have the same dimensions as the one specified
    in the arguments, a ValueError is raised otherwise.
    """

    v = copy.deepcopy(v)
    
    if is_unitarray_dict(v) : 
        v = unitarray_from_dict(v)
    elif isdict(v) :
        for k in list(v) :
            if is_unitarray_dict(v[k]) :
                v[k] = unitarray_from_dict(v[k])
    
    v_out = copy.deepcopy(v)

    # if v is a dict, all value numbers are replaced by UnitValues
    if type(v) == dict :

        v_out = {}

        # raise an error if dicts are not accepted for the property
        if not accepts_dict : 
            raise ValueError("this property does not accepts dict as input value")

        for k in list(v) :
            if isnumber(v[k]) or isstr(v[k]) or type(v[k]) == UnitValue : 
                if not accepts_singlevalue : 
                    raise ValueError("this property does not accepts single values as input")
                    
                for ki in k.split(",") : 
                    v_out[ki.strip()] = UnitValue(v[k], Units(units_system, units_dimensions), convert=False)

            elif isarray(v[k]) or type(v[k]) == UnitArray :
                if not accepts_array : 
                    raise ValueError("this property does not accepts arrays as input values")
                
                for ki in k.split(",") :
                    v_out[ki.strip()] = UnitArray(v[k], Units(units_system, units_dimensions), convert=False)
                
            else :
                raise TypeError("unexpected property dict value type.")
    
    # if v is a number, it is converted to a UnitValue
    elif isnumber(v) or isstr(v) or type(v) == UnitValue :
        if not accepts_singlevalue : 
            raise ValueError("this property does not accepts single values as input")
        
        v_out = UnitValue(v, Units(units_system, units_dimensions), convert=False)
    
    # if v is an array, it is converted to a UnitArray
    elif isarray(v) or type(v) == UnitArray:
        if not accepts_array : 
            raise ValueError("this property does not accepts arrays as input value")

        v_out = UnitArray(v, Units(units_system, units_dimensions), convert=False)
    
    # if none of the conditions above is met, a ValueError is raised.
    else :
        raise TypeError("unexpected property type. must be dict, number or UnitValue.")

    return v_out

def format_unitvar_for_save(v, units_system) :
    """
    returns str(v) if v's units system differs from us, the default units system, v.value otherwise.'
    """
    
    if type(v) == UnitValue : 
        return str(v)
    elif type(v) == UnitArray : 
        return unitarray_to_dict(v)
    elif type(v) == dict : 
        d = {}
        for k in list(v) :
            if type(v[k]) == UnitValue : 
                d[k] = str(v[k])
            elif type(v[k]) == UnitArray : 
                d[k] = unitarray_to_dict(v[k])            
        return d
                    
def process_input_dict_keys(d, synonyms, policy="error") :
    """
    process a dict given in input for functions such as thing_from_dict.
    
    :param d: dict to be processed
    :type d: dict
    :param synonyms: different key strings expected as well as their synonyms
        it must be a list of striung key synonyms list. Each key with one synonym in the list will be replaced by the 
        first key in the synonym list. synonyms, should not contain any duplicated key.
        ie.
            
        synonyms = [
        ["chemostates", chstt],
        ["cell_volume", cell_vol, cell_v], 
        ...]
        
    :type synonyms: array[array[str]]
    :param policy: indicates how to manage unexpected or duplicated keys.
        expects :
            
        * "silent" : issues are ignored
        * "warning" : issues are presented through a console print
        * "error" : issues raise a ValueError exception
        
    :type policy: str
    :returns: processed dict with synonyms replaced.
    :rtype: dict
    """
    
    #first step : search for unexpected keys
    for k in list(d) : 
        key_found = False
        for s in synonyms : 
            if k in s : 
                key_found = True
                break
        if not key_found : 
            if policy == "error" :
                raise ValueError("unexcpected key \"" + k + "\" in input dictionary.")
            elif policy == "warning" :
                print("warning : unexcpected key \"" + k + "\" in input dictionary.")

    #second step : search for duplicated synonyms
    for s in synonyms : 
        synonyms_found = []
        for k in list(d) : 
            if k in s :
                synonyms_found.append(k)
        if len(synonyms_found)>1 : 
            if policy == "error" :
                raise ValueError("dictionnary key " + synonyms_found + " in input dictionary are synonyms.")
            elif policy == "warning" :
                print("warning : dictionnary key " + synonyms_found + " in input dictionary are synonyms.")
    
    #third step : replace synonyms
    for s in synonyms : 
        for k in list(d) : 
            if k in s :
                d[s[0]] = d[k]
    
    return d

def retrive_units_system_from_dict(d, default, parent_units_system) : 
    v = d.get("units", default)
    if isstr(v) :
        if v=="default" : 
            return UnitsSystem()
        elif v=="inherit" :
            return parent_units_system
        else :
            raise ValueError("Undefined string value \""+v+"\" for \"units\".")
    elif isdict(v) :
        return unitssystem_from_dict(v)
    else :
        raise TypeError("dict \"units\" key must be a dict or a string.")

def assert_string_is_a_valid_label(l) : 
    if isnone(l) : 
        return
    
    for c in l : 
        if c in string.whitespace : 
            raise ValueError("unexpected whitespace character in label.")
        if c in "+" :
            raise ValueError("unexpected character \"+\" in label.")
        if c.count("->") > 0 :
            raise ValueError("label cannot contain the \"->\" sequence.")

def get_value_in_env(value, environment, default) :
    if isdict(value) :
        if environment in list(value) : 
            return value[environment]
        elif "default" in list(value) :
            return value["default"]
        else :
            return default
    else :
        return value