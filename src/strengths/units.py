import numpy as np
import copy as cpy

from strengths.typechecking import *
from strengths import constants
import strengths.value_processing as valproc
from strengths import filepath

_units_labels_dict = {
    "space"    : ["km", "m", "dm", "cm", "mm", "dmm", "cmm", "µm", "nm", "pm", "fm"],
    "time"     : ["h", "min", "s", "ds", "cs", "ms", "µs", "ns", "ps", "fs"],
    "quantity" : ["kmol", "mol", "dmol", "cmol", "mmol", "µmol", "nmol", "pmol", "fmol", "molecule"],
    "density"  : ["kM", "M", "dM", "cM", "mM", "µM", "nM", "pM", "fM", "mM"],
    "volume"   : ["kL", "L", "mL", "µL", "nL", "pL", "fL"]
    }

_units_conversion_dict = {
    "space"    : {
        "km" : 1e3,
        "m"  : 1,
        "dm" : 1e-1,
        "cm" : 1e-2,
        "mm" : 1e-3,
        "dmm": 1e-4,
        "cmm": 1e-5,
        "µm" : 1e-6,
        "nm" : 1e-9,
        "pm" : 1e-12,
        "fm" : 1e-15
        },
    "time"     : {
        "h"  : 3600,
        "min": 60,
        "s"  : 1,
        "ds" : 1e-1,
        "cs" : 1e-2,
        "ms" : 1e-3,
        "µs" : 1e-6,
        "ns" : 1e-9,
        "ps" : 1e-12,
        "fs" : 1e-15
        },
    "quantity" : {
        "kmol" : 1e3*constants.avogadro_number(),
        "mol"  : 1*constants.avogadro_number(),
        "dmol" : 1e-1*constants.avogadro_number(),
        "cmol" : 1e-2*constants.avogadro_number(),
        "mmol" : 1e-3*constants.avogadro_number(),
        "µmol" : 1e-6*constants.avogadro_number(),
        "nmol" : 1e-9*constants.avogadro_number(),
        "pmol" : 1e-12*constants.avogadro_number(),
        "fmol" : 1e-15*constants.avogadro_number(),
        "molecule" : 1
        }
    }

_default_units_system_dict = {
    "space" : "µm",
    "time" : "s",
    "quantity" : "molecule"
    }

def compute_conversion_factor(su_src, su_dst, sdim) :
    """
    Compute the conversion factor between two units system for a given units dimension.

    :param su_src: units system assocuated with value
    :type su_src: UnitsSystem
    :param su_dst: units system in which value should be converted
    :type su_dst: UnitsSystem
    :param sdim: units dimensions associated with value.
    :type sdim: UnitsDimensions
    :returns: converted value
    :rtype: number
    """

    f = 1
    for k in su_src.keys() :
        f *= (_units_conversion_dict[k][su_src[k]]/_units_conversion_dict[k][su_dst[k]])**sdim[k]

    return f

def convert_value(value, su_src, su_dst, sdim) :
    """
    Converts a value expressed in a given units system (su_src) into another units system (su_dst).

    :param value: numerical value of the quantity to be converted
    :type value: number
    :param su_src: units system assocuated with value
    :type su_src: UnitsSystem
    :param su_dst: units system in which value should be converted
    :type su_dst: UnitsSystem
    :param sdim: units dimensions associated with value.
    :type sdim: UnitsDimensions
    :returns: converted value
    :rtype: number
    """

    return value*compute_conversion_factor(su_src, su_dst, sdim)

def unitssystem_from_dict(d) :
    """
    Creates a UnitsSystem from the dictionary d.

    :param d: dict from which the instance qhould be created
    :type d: dict
    :returns: UnitsSystem created from d
    :rtype: UnitsSystem
    """

    if type(d) != dict :
        raise ValueError("d must be a dict.")

    d = valproc.process_input_dict_keys(d, [
                ["space"],
                ["time"],
                ["quantity"]
            ]
        )

    return UnitsSystem(
        space    = d["space"],
        time     = d["time"],
        quantity = d["quantity"]
        )

def unitsdimensions_from_dict(d) :
    """
    Creates a UnitsDimensions from the dictionary d.

    :param d: dict from which the instance should be created
    :type d: dict
    :returns: UnitsDimensions created from d
    :rtype: UnitsDimensions
    """

    if type(d) != dict :
        raise ValueError("d must be a dict.")

    d = valproc.process_input_dict_keys(d, [
                ["space"],
                ["time"],
                ["quantity"]
            ]
        )

    return UnitsDimensions(
        space    = d["space"],
        time     = d["time"],
        quantity = d["quantity"]
        )

def unitssystem_to_dict(v) :
    if type(v) != UnitsSystem :
        raise ValueError("v must be a UnitsSystem.")
    return {
        "space"    : v["space"],
        "time"     : v["time"],
        "quantity" : v["quantity"]
        }

def unitsdimensions_to_dict(v) :
    if type(v) != UnitsDimensions :
        raise ValueError("v must be a UnitsDimensions.")
    return {
        "space"    : v["space"],
        "time"     : v["time"],
        "quantity" : v["quantity"]
        }

def convert_unitvalue(v, u) :
    """
    Converts the UnitValue v to the Units u.

    :param v: UnitValue to be converted
    :type v: UnitValue
    :param u: target units
    :type u: str, dict, UnitValue, Units or UnitsSystem
    :returns: converted UnitValue
    :rtype: UnitValue
    """

    su_dst = 0
    if type(u) == str :
        su_dst = parse_units(u)
        if su_dst.dim != v.units.dim :
            raise ValueError("unit conversion must happen in the same dimension.")
        su_dst = su_dst.sys
    elif type(u) == UnitValue :
        if u.units.dim != v.units.dim :
            raise ValueError("unit conversion must happen in the same dimension.")
        su_dst = u.units.sys
    elif type(u) == Units :
        if u.dim != v.units.dim :
            raise ValueError("unit conversion must happen in the same dimension. Trying to convert " + str(u.dim) + " to " + str(v.units.dim) + ".")
        su_dst = u.sys
    elif type(u) == UnitsSystem :
        su_dst = u
    elif type(u) == dict :
        su_dst = unitssystem_from_dict(u)
    else :
        raise ValueError("")

    return UnitValue(convert_value(v.value, v.units.sys, su_dst, v.units.dim), Units(su_dst, v.units.dim))

class _UnitsComponentDict :

    def keys(self) :
        return ["space", "time", "quantity"]

    def __str__(self) :
        return str(self.space) + ", " + str(self.time) + ", " + str(self.quantity)

    @property
    def space(self) :
        return self._space

    @property
    def time(self) :
        return self._time

    @property
    def quantity(self) :
        return self._quantity

    @space.setter
    def space(self, v) :
        self._check_space(v)
        self._space = v

    @time.setter
    def time(self, v) :
        self._check_time(v)
        self._time = v

    @quantity.setter
    def quantity(self, v) :
        self._check_quantity(v)
        self._quantity = v

    def __getitem__(self, k) :
        if type(k) != str :
            raise ValueError("k must be a string.")
        elif k=="space"    : return self.space
        elif k=="time"     : return self.time
        elif k=="quantity" : return self.quantity
        else : raise ValueError("k must be \"space\", \"time\" or \"quantity\".")

    def __setitem__(self, k, v) :
        if type(k) != str :
            raise ValueError("k must be a string.")
        elif k=="space"    : self.space = v
        elif k=="time"     : self.time = v
        elif k=="quantity" : self.quantity = v
        else : raise ValueError("k must be \"space\", \"time\" or \"quantity\".")

    def copy(self) :
        """
        Returns a deepcopy of the instance.

        .. code:: python

            instance.copy()

            # is equivalent to
            # import copy

            copy.deepcopy(instance)

        """

        return cpy.deepcopy(self)

    def __eq__(self, v) :
        if type(v) == dict or type(v) == type(self) :
            return (
                self.space    == v["space"] and
                self.time     == v["time"] and
                self.quantity == v["quantity"]
                )
        else :
            return None

class UnitsSystem(_UnitsComponentDict) :
    """
    descibes a system of fundamental units,
    one for time (h, min, s, etc.), one for space (m, mm, µm, etc.) and one for quantity (moleclue, mol, etc.).

    :param space: space distance units
    :type space: str
    :param time: time distance units
    :type time: str
    :param quantity: units of quantity of matter
    :type quantity: str
    """

    def __init__(self,
                 space    = _default_units_system_dict["space"],
                 time     = _default_units_system_dict["time"],
                 quantity = _default_units_system_dict["quantity"]) :
        self.space = space
        self.time = time
        self.quantity = quantity

    def _check_space(self, v) :
        if type(v) != str :
            raise ValueError("v must be a string.")
        if not v in _units_labels_dict["space"] :
            raise ValueError(v+" is not an accepted space unit.")

    def _check_time(self, v) :
        if type(v) != str :
            raise ValueError("v must be a string.")
        if not v in _units_labels_dict["time"] :
            raise ValueError(v+" is not an accepted time unit.")

    def _check_quantity(self, v) :
        if type(v) != str :
            raise ValueError("v must be a string.")
        if not v in _units_labels_dict["quantity"] :
            raise ValueError(v+" is not an accepted quantity unit.")

class UnitsDimensions(_UnitsComponentDict) :
    """
    descibes the dimensions associated with a system of units.

    :param space: space distance exponent
    :type space: str
    :param time: time distance exponent
    :type time: str
    :param quantity: exponent of quantity of matter
    :type quantity: str
    """

    def __init__(self, space=0, time=0, quantity=0) :
        self.space = space
        self.time = time
        self.quantity = quantity

    def _check_space(self, v) :
        if type(v) != int :
            raise ValueError("v must be an integer.")

    def _check_time(self, v) :
        if type(v) != int :
            raise ValueError("v must be an integer.")

    def _check_quantity(self, v) :
        if type(v) != int :
            raise ValueError("v must be an integer.")

def density_units_dimensions() :
    """
    returns the units dimensions associated with a density (space^-3, time^0, quantity^1).
    """

    return UnitsDimensions(space = -3, time = 0, quantity = 1)

def surface_units_dimensions() :
    """
    returns the units dimensions associated with a surface (space^2, time^0, quantity^0).
    """

    return UnitsDimensions(space = 2,  time = 0, quantity = 0)

def volume_units_dimensions() :
    """
    returns the units dimensions associated with a volume (space^3, time^0, quantity^0).
    """

    return UnitsDimensions(space = 3,  time = 0, quantity = 0)

def quantity_units_dimensions() :
    """
    returns the units dimensions associated with a quantity of matter (space^0, time^0, quantity^1).
    """

    return UnitsDimensions(space = 0,  time = 0, quantity = 1)

def space_units_dimensions() :
    """
    returns the units dimensions associated with a space distance (space^1, time^0, quantity^0).
    """

    return UnitsDimensions(space = 1,  time = 0, quantity = 0)

def time_units_dimensions() :
    """
    returns the units dimensions associated with a time distance (space^0, time^1, quantity^0).
    """

    return UnitsDimensions(space = 0,  time = 1, quantity = 0)

class Units :
    """
    describe the units of a variable.
    there is two ways to initialize a Units object :

    1) by specifying both sys and dim, respecitvely with UnitsSystem and UnitsDimensions objects, or dicts that can be converted to
        such objects.
        ie.

        .. code:: python

            Units(UnitSystem(space="m"), UnitsDimensions(space=1))         # OK, "m"
            Units({"space":"m", "time":"s", "quantity":"molecule"},
                      {"space":1, "time":0, "quantity":0 })                # OK, "m"
            Units(sys=UnitSystem(space="m"), dim=UnitsDimensions(space=1)) # OK, "m"
            Units(sys={"space":"m", "time":"s", "quantity":"molecule"},
                      dim={"space":1, "time":0, "quantity":0 })            # OK, "m"
            Units(1, 2)                                                    # error, 1 and 2 cannot be converted to
                                                                           # UnitsSystem and UnitsValues objects
            Units("1 s", UnitsDimensions(time=1))                          # error, dim is not None, and "1 s" cannot

    2) by specifying only sys or setting dim to None. sys must then be a characted string that can be parsed by parse_units.
        ie.

        .. code:: python

            Units("1 m")                     # OK, "m"
            Units("1 m", None)               # OK, "m"
            Units(sys="1 m")                 # OK, "m"
            Units(sys="1 m", dim=None)       # OK, "m"

            Units(1)                         # error, as dim=None, sys must be a string.
            Units(sys=UnitSystem(space="m")) # error, same reason as above

    :param sys: string representation of the units if dim is None/unspecified, UnitsSystem or equivalent dictonnary otherwise
    :type sys: str, UnitSystem or dict
    :param dim: None or UnitsDimensions or equivalent dictionnary
    :type dim: None, UnitsDimensions or dict
    """

    def __init__(self, sys, dim=None) :
        if isnone(dim) :
            u = parse_units(sys)
            self.sys = u.sys
            self.dim = u.dim
        else :
            self.sys = sys
            self.dim = dim

    @property
    def sys(self) :
        """
        units associated with space, time, etc.
        """

        return self._sys

    @property
    def dim(self) :
        """
        dimensions associated with time, space, etc.
        """

        return self._dim

    @sys.setter
    def sys(self, v) :
        if type(v) == UnitsSystem :
            self._sys = v.copy()
        elif type(v) == dict :
            self._sys = unitssystem_from_dict(v)
        else :
            raise ValueError("v must be a UnitsSystem instance.")

    @dim.setter
    def dim(self, v) :
        if type(v) == UnitsDimensions :
            self._dim = v.copy()
        elif type(v) == dict :
            self._dim = unitsdimensions_from_dict(v)
        else :
            raise ValueError("v must be a UnitsDimensions instance.")


    def __str__(self) :
        s = []
        for k in self.sys.keys() :
            if self.dim[k] != 0 :
                if self.dim[k] != 1 :
                    s.append(self.sys[k] + str(self.dim[k]))
                else :
                    s.append(self.sys[k])
        out = ""
        for i in range(len(s)) :
            out += s[i]
            if i!=len(s)-1:
                out+="."
        return out

    def invert(self) :
        """
        return inverted units.
        ie. µm2/s -> s/µm-2
        """

        invdim = UnitsDimensions()
        for k in self.dim.keys() :
            invdim[k] = -self.dim[k]
        return Units(self.sys, invdim)

    def multiply(self, u) :
        """
        return the product of the instance units by the Units u.
        ie. m2/s * m-1 -> m/s
        """

        if type(u) != Units :
            raise ValueError("u must be Units.")
        if self.sys != u.sys:
            raise ValueError("cannot multiply units of different systems.")
        sdim = UnitsDimensions()

        for k in self.dim.keys() :
            sdim[k] = self.dim[k] + u.dim[k]
        return Units(self.sys, sdim)

    def raiseto(self, e) :
        """
        return the product of the instance units raise to the power of e.
        e doesnt have to be an integer, but dimensions raised to the power of e have to.
        ie. (m/s)^3 = m3/s3
        """

        rdim = UnitsDimensions()
        for k in self.dim.keys() :
            rdim[k] = int(self.dim[k] * e)
            if self.dim[k]*e - rdim[k] != 0 :
                raise ValueError("error : when raising a UnitValue to some power, the resulting unit dimensions must be integers.")

        return Units(self.sys, rdim)

    def __eq__(self, v) :
        if type(v) == str :
            v = parse_units(v)
        if type(v) == Units :
            for k in ["space", "time", "quantity"] :
                if (self.dim[k] != v.dim[k]) :
                    return False
                if (self.dim[k] != 0) and (self.sys[k] != v.sys[k]) :
                    return False
            return True
        else :
            raise("v must be a Units or str.")

    def copy(self) :
        """
        Returns a deepcopy of the instance.

        .. code:: python

            instance.copy()

            # is equivalent to
            # import copy

            copy.deepcopy(instance)

        """

        return cpy.deepcopy(self)

def parse_units(s) :
    """
    Create Units from the character string s.
    ie. "µM/s" -> µmol/dm3/s

    syntax is simple :
    integral units exponent, positive or negative, must be appended to the units
    ie. µm2 or m-3
    each units is separated by "." or "/"
    ie. m/s/mol -> m.s-1.mol-1
    non fundamental units are automatically translated into fundamental units :
    ie. L -> dm3, µM/s -> µmol/dm3/s
    """
    if type(s) != str :
        raise "s must be a string."
    
    s = s.replace("um", "µm")
    s = s.replace("us", "µs")
    s = s.replace("umol", "µmol")
    s = s.replace("uL", "µL")
    s = s.replace("uM", "µM")
    
    s = s.strip()
    if s == "" :
        return Units(UnitsSystem(), UnitsDimensions())

    def get_unit_type(unitstr) :
        for k in _units_labels_dict.keys() :
            if unitstr in _units_labels_dict[k] :
                return k
        return None

    def get_volume_fundamental_unit(volstr) :
        if   volstr == "kL" : return "m"
        elif volstr == "L"  : return "dm"
        elif volstr == "mL" : return "cm"
        elif volstr == "µL" : return "mm"
        elif volstr == "nL" : return "dmm"
        elif volstr == "pL" : return "cmm"
        elif volstr == "fL" : return "µm"
        else : raise(ValueError("unexpected unit"))

    def get_concentration_fundamental_units(constr) :
        if   constr == "kM" : return "kmol", "dm"
        elif constr == "M"  : return "mol",  "dm"
        elif constr == "dM" : return "dmol", "dm"
        elif constr == "cM" : return "cmol", "dm"
        elif constr == "mM" : return "mmol", "dm"
        elif constr == "µM" : return "µmol", "dm"
        elif constr == "nM" : return "nmol", "dm"
        elif constr == "pM" : return "pmol", "dm"
        elif constr == "fM" : return "fmol", "dm"
        else : raise(ValueError("unexpected unit"))

    blocks = [[".", "", ""]]
    exp = False
    count = 0

    for i in s :
        if i == '.' or i == "/" :
            blocks.append([i,"",""])
            exp=False
            count+=1
        else:
            if i in ["-","0","1","2","3","4","5","6","7","8","9"] :
                exp=True
            if exp :
                blocks[count][2] += i
            else :
                blocks[count][1] += i
    for b in blocks :
        if b[2] == "" :
            b[2] = "1"
        b[2] = int(b[2])
        if b[0] == "/" :
            b[2] = -b[2]

    sys = {
        "space" : None,
        "time" : None,
        "quantity" : None
        }
    dim = {
        "space" : 0,
        "time" : 0,
        "quantity" : 0
        }

    def addunit(field, su, se) :
        if sys[field] == None or sys[field] == su:
            sys[field] =  su
            dim[field] += se
        else:
            raise Exception("incompatible units "+sys[dim]+" and "+su+".")

    for b in blocks :
        unittype = get_unit_type(b[1])
        if unittype == None :
            raise Exception("undefined unit \""+b[1]+"\".")
        else :
            if unittype == "space" :
                addunit("space", b[1], b[2])
            if unittype == "time" :
                addunit("time", b[1], b[2])
            if unittype == "quantity" :
                addunit("quantity", b[1], b[2])
            if unittype == "volume" :
                addunit("space", get_volume_fundamental_unit(b[1]), b[2]*3)
            if unittype == "density" :
                addunit("space", get_concentration_fundamental_units(b[1])[1], b[2]*-3)
                addunit("quantity", get_concentration_fundamental_units(b[1])[0], b[2])

    if sys["space"]    == None : sys["space"]    = _default_units_system_dict["space"]
    if sys["time"]     == None : sys["time"]     = _default_units_system_dict["time"]
    if sys["quantity"] == None : sys["quantity"] = _default_units_system_dict["quantity"]

    return Units(sys, dim)

class UnitValue :
    """
    Describes a physical value with its units.
    there is two ways to initialize a UnitValue :

    1) by specifying separately its numerical value and its units :
        value must then be a number, and units a Units obect or a string that can be
        parsed as such :
        ie.

        .. code:: python

            UnitValue(1, "m")                           # ok, 1 m
            UnitValue(1, Units("m"))                    # ok, 1 m
            UnitValue(1, Units(UnitsSystem(space="m"),
                         UnitsDimensions(space=1)))     # ok, 1 m

            UnitValue("1", "m"))                        # error : "1" is not a number
            UnitValue(1, []))                           # error : [] is not a Units obect nor a string

    2) by specifying only value or setting units to None. value must then be a string that can be parsed to
        a UnitValue by parse_unitvalue :
        ie.

        .. code:: python

            UnitValue("1 m")         # ok, 1 m
            UnitValue("1 m", None)   # ok, 1 m
            UnitValue(1)             # wrong, 1 is not a string
            UnitValue("a")           # wrong, "a" cannot be parsed into a UnitValue.

    :param value: if units are defined (not None), the numerical value of the UnitValue,
        otherwise, a string that can be parsed as a UnitValue.
    :type value: number or str
    :param units: units of the variable or None (default : None)
    :type units: Units, str or None
    :param convert: tells is the value, if it is already a UnitValue (or a string representing one), should be converted to the specified units (default : True)
    :type convert: bool
    """

    def __init__(self, value, units=None, convert=True) :

        if type(units) == str :
            units = Units(units)
        elif type(units) == Units :
            units = units.copy()
        elif isnone(units) : 
            pass
        else :
            raise TypeError("units must be a str, an instance of the Units class or None.")

        if isnumber(value) : 
            if type(units) == Units : # defined units
                self.value = value
                self.units = units
            else :                    # no units
                self.value = value
                self.units = Units("")
 
        elif isstr(value) or type(value) == UnitValue :
            if isstr(value) :
                value = parse_unitvalue(value)
        
            if type(units) == Units : # defined units
                if convert : #units compatibility check is done in the conversion.
                    value = value.convert(units)
                else : #only the units compatibility check must be done
                    if value.units.dim != units.dim : 
                        raise ValueError("the value "+str(value)+" and the units " + str(units)+ " have incompatble dimensions.")
                    
            self.value = value.value
            self.units = value.units  
                
        else :
            raise TypeError("unsupported value type.")

    @property
    def value(self) :
        """
        Numerical value of the pysical quantity (float).
        """

        return self._value

    @value.setter
    def value(self, v) :
        if isnumber(v) :
            self._value = float(v)
        else :
            raise ValueError("UnitValue's value must be a number.")

    @property
    def units(self) :
        """
        Units of the physical quantity (:py:class:`Units`).
        can be set with a string or an instance of :py:class:`Units`.

        """

        return self._units

    @units.setter
    def units(self, v) :
        if type(v) == str :
            v = parse_units(v)
        if type(v) == Units :
            self._units = v.copy()
        else :
            raise ValueError("UnitValue's unit must be an instance of Units or a unit string.")

    def __str__(self) :
        return str(self.value) + " " + self.units.__str__()

    def __repr__(self) :
        return "UnitValue(\""+str(self)+"\")"

    def convert(self, u) :
        """
        Return the UnitValue converted into the units u.
        u must have compatible dimensions.
        same as convert_value(v, u), where v is the UnitValue instance.
        """

        return convert_unitvalue(self, u)

    def _sum(self, v) :
        """
        returns self + v.
        """

        if type(v)==UnitValue :
            if self.units.dim == v.units.dim :
                return UnitValue(self.value + v.convert(self.units.sys).value, self.units)
            else :
                raise ValueError("addition of two UnitValue with different units dimensions is not supported.")

        elif isnumber(v) :
            return UnitValue(self.value+v, self.units)

        elif type(v)==UnitArray :
            if self.units.dim != v.units.dim :
                raise ValueError("addition of a UnitValue with a UnitArray with different units dimensions is not supported.")
            v = v.convert(self.units.sys)
            return UnitArray([self.value + v.value[i] for i in range(len(v))], self.units)

        else :
            raise ValueError("unexpected types for a sum.")

    def _product(self, v) :
        """
        returns self * v.
        """

        if type(v)==UnitValue:
            v = v.convert(self.units.sys)
            return UnitValue(self.value * v.value, self.units.multiply(v.units))

        elif isnumber(v) :
            return UnitValue(self.value * v, self.units)

        elif type(v)==UnitArray :
            v = v.convert(self.units.sys)
            return UnitArray([self.value * v.value[i] for i in range(len(v))], self.units.multiply(v.units))

        else :
            raise ValueError("unexpected types for a product.")

    def _modulo(self, mod):
        """
        returns self % mod.
        """

        if type(mod)==UnitValue:
            if self.units.dim == mod.units.dim :
                return UnitValue(self.value % mod.convert(self.units.sys).value, self.units)
            else :
                raise ValueError("modulo of UnitValues with different unit dimensions is not supported.")

        elif isnumber(mod) :
            return UnitValue(self.value % mod, self.units)

        elif type(mod)==UnitArray :
            if self.units.dim != mod.units.dim :
                raise ValueError("addition of a UnitValue with a UnitArray with different units dimensions is not supported.")
            mod = mod.convert(self.units.sys)
            return UnitArray([self.value % mod.value[i] for i in range(len(mod))], self.units)

        else :
            raise ValueError("unexpected types for a modulo.")

    def _rmodulo(self, v):
        """
        returns v % self.
        """

        if type(v)==UnitValue:
            if self.units.dim == v.units.dim :
                return UnitValue(v.convert(self.units.sys).value%self.value, self.units)
            else :
                raise ValueError("modulo of UnitValues with different unit dimensions is not supported.")
        elif isnumber(v) :
            return UnitValue(v%self.value, self.units)

        elif type(v)==UnitArray :
            if self.units.dim != v.units.dim :
                raise ValueError("addition of a UnitValue with a UnitArray with different units dimensions is not supported.")
            v = v.convert(self.units.sys)
            return UnitArray([v.value[i] % self.value for i in range(len(v))], self.units)

        else :
            raise ValueError("unexpected types for a modulo.")

    def invert(self) :
        """
        Invert the variable (returns 1/v).
        """

        return UnitValue(1/self.value, self.units.invert())

    def __add__    (self, v) :
        """
        defines self+v.
        Addition must be operated betwen UnitValue/UnitArray objects with the same units dimension.
        if v is not a UnitValue/UnitArray, is is conseidered to be in the same units than self.
        if v is a UnitArray, self is added to v to everyelement of v.
        Conversion of v is implicitely done if required.

        ie.

        .. code:: python

            UnitValue(1, "m") + 1                     # ok, == UnitValue(2, "m")
            UnitValue(1, "m") + UnitValue(1000, "mm") # ok, == UnitValue(2, "m")
            UnitValue(1, "m") + UnitArray([1,2], "m") # ok, == UnitArray([2,3], "m")
            UnitValue(1, "m") + UnitValue(1, "m/s")   # ValueError

        """

        return self._sum(v)

    def __radd__    (self, v) :
        """
        defines v+self. same as self.__add__(v).
        """

        return self._sum(v)

    def __sub__    (self, v) :
        """
        defines self-v. same as self.__add__(-v).
        """

        return self._sum(_neg(v))

    def __rsub__    (self, v) :
        """
        defines v-self. sames as -self.__sub__(v).
        """

        return _neg(self._sum(_neg(v)))

    def __mul__    (self, v) :
        """
        defines self*v.
        if v is not a UnitValue/UnitArray, is is conseidered to be a unitless factor.
        if v is a UnitArray, self multiplies every element of v.
        Conversion of v is implicitely done if required.

        ie.

        .. code:: python

            UnitValue(1, "µm") * UnitValue(1, "µm")    # ok == UnitValue(1, "µm2")
            UnitValue(1, "µm") * UnitValue(1, "µm-1")  # ok == UnitValue(1, "")
            UnitValue(1, "µm") * 2                     # ok == UnitValue(2, "µm")
            UnitValue(1, "µm") * UnitValue(1, "mol/s") # ok == UnitValue(1, "µm.mol/s")
            UnitValue(1, "µm") * UnitArray([1,2], "s") # ok == UnitArray([1,2], "µm.s")

        """

        return self._product(v)

    def __rmul__    (self, v) :
        """
        defines v*self. same as self.__mul__(v).
        """

        return self._product(v)

    def __truediv__(self, v) :
        """
        defines self/v.
        if v is not a UnitValue/UnitArray, is is conseidered to be a unitless factor.
        Conversion is implicitely done if required.

        ie.

        .. code:: python

            UnitValue(1, "µm") / UnitValue(1, "µm")    # ok == UnitValue(1, "")
            UnitValue(1, "µm") / UnitValue(1, "µm-1")  # ok == UnitValue(1, "µm2")
            UnitValue(1, "µm") / 2                     # ok == UnitValue(0.5, "µm")
            UnitValue(1, "µm") / UnitValue(1, "mol/s") # ok == UnitValue(1, "µm.s/mol")
            UnitValue(1, "µm") / UnitArray([1,2], "s") # ok == UnitArray([1,0.5], "µm/s")
        """

        return self._product(_inv(v))

    def __rtruediv__(self, v) :
        """
        defines v/self. same as self.invert()*v.
        """

        return self.invert()._product(v)

    def __mod__    (self, v) :
        """
        defines self%v.
        if v is not a UnitValue/UnitArray, self.value%v is returned.
        if v is a UnitValue or a UnitArray, v and self must have the same units dimensions, a ValueError is returned otherwise.
        if v is a UnitArray, the modulo is operated element wise.
        The operation is applied to the numeric values after conversion in the same units system.

        ie.

        .. code:: python

            UnitValue("3 µm") % 2                      #ok, == UnitValue("1 µm")
            UnitValue("3 µm") % UnitValue("2 µm")      #ok, == UnitValue("1 µm")
            UnitValue("3 µm") % UnitArray([2,3], "µm") #ok, == UnitArray([1,0], "µm")

            UnitValue("3 µm") % UnitValue(1, "s")      #wrong, different units dimensions

        """

        return self._modulo(v)

    def __rmod__    (self, v) :
        """
        defines v%self.
        type and units requirements are the same than for __mod__.

        ie.

        .. code:: python

            3 % UnitValue("2 µm")                      #ok, == UnitValue("1 µm")

        """

        return self._rmodulo(v)

    def __pow__    (self, v) :
        """
        defines self**v.
        v must be a number.
        pow is valid as long as the resluting units dimensions are integers.

        ie.

        .. code: python

            Units(1, "m")**3      # ok, 3*1=3 is still an integer
            Units(1, "m3")**(1/3) # ok, 3/3=1 is still an integer
            Units(1, "m")**(1/3)  # ValueError, 1/3 is not an integer
            Units(1, "")**(1/3)   # ok, 0/3=0 is still an integer

        """

        if isnumber(v) :
            return UnitValue(self.value**v, self.units.raiseto(v))
        else:
            raise TypeError("power v must be a number")

    def __rpow__    (self, v) :
        """
        defines v**self is an, invalid opeartion.
        raising something to the power of a UnitValue is not a supported opeartion.
        raise a ValueError.

        """

        raise NotImplementedError("raising something to a UnitValue is not a supported operation.")

    def __neg__(self) :
        """
        defines the unary operation -self.
        returns a copy of self, but with an opposed value.
        ie.

        .. code:: python

            -UnitValue("2 µm")  # ok, == UnitValue("-2 µm")
            -UnitValue("-2 µm") # ok, == UnitValue("2 µm")

        """

        return UnitValue(-self.value, self.units)

    def __pos__(self) :
        """
        defines the unary operation +self.
        returns a copy of self.
        ie.

        .. code:: python

            +UnitValue("2 µm")  # ok, == UnitValue("2 µm")
            +UnitValue("-2 µm") # ok, == UnitValue("-2 µm")

        """

        return self.copy()

    def __abs__(self) :
        """
        defines abs(self).
        simply returns UnitValue(abs(self.value), self.units).
        ie.

        .. code:: python

            abs(UnitValue("-2 µm")) # ok, == UnitValue("2 µm")
            abs(UnitValue("2 µm"))  # ok, == UnitValue("2 µm")

        """

        return UnitValue(abs(self.value), self.units)

    def __eq__(self, v) :
        """
        defines self == v or v == self.

        if v is not a UnitValue, self.value == v is returned.
        if v is a UnitValue, values are considered equal if they can be expressed
        with the same value and units through conversion.

        ie.

        .. code:: python

            UnitValue(1, "µm") == 1                     # True
            UnitValue(1, "µm") == UnitValue(1, "µm")    # True
            UnitValue(1, "µm") == UnitValue(1e-3, "mm") # True

            UnitValue(1, "µm") == UnitValue(1e-6, "mm") # False
            UnitValue(1, "µm") == UnitValue(1, "µm/s")  # False

        """

        if type(v) == UnitValue :
            if self.units.dim == v.units.dim :
                return (self.value == v.convert(self.units.sys).value)
            else :
                return False
        elif isnumber(v):
            return self.value == v
        else :
            return False

    def __neq__(self, v) :
        """
        defines self != v or v != self.
        defined as not self.__eq__(v), see the __eq__ method for more information.

        """

        return not self.__eq__(v)

    def __gt__(self, v) :
        """
        defines self > v or v < self.

        if v is not a UnitValue, self.value > v is returned.
        if v is a UnitValue, numerical value converted to the same units system
        are compared. A ValueError is raised if self and v have different units dimensions.
        a TypeError is raised if v is not a UnitValue or a number.

        ie.

        .. code:: python

            UnitValue(1, "µm") > UnitValue(1, "µm")   # False
            UnitValue(1, "µm") > UnitValue(1, "mm")   # True
            UnitValue(1, "µm") > 0.5                  # True
            UnitValue(1, "µm") > UnitValue(1, "µm/s") # ValueError

        """

        if type(v) == UnitValue :
            if self.units.dim == v.units.dim :
                return (self.value > v.convert(self.units.sys).value)
            else :
                raise ValueError("cannot compare UnitValues with different units dimensions.")
        elif isnumber(v):
            return self.value > v
        else :
            return TypeError("UnitValue can only be compared with numbers and UnitValues with the same units dimensions.")

    def __ge__(self, v) :
        """
        defines self >= v or v <= self.

        if v is not a UnitValue, self.value >= v is returned.
        if v is a UnitValue, numerical value converted to the same units system
        are compared. A ValueError is raised if self and v have different units dimensions.
        a TypeError is raised if v is not a UnitValue or a number.

        ie.

        .. code:: python

            UnitValue(1, "µm") >= UnitValue(1, "µm")   # True
            UnitValue(1, "µm") >= UnitValue(1, "mm")   # True
            UnitValue(1, "µm") >= 0.5                  # True
            UnitValue(1, "µm") >= UnitValue(1, "µm/s") # ValueError

        """

        if type(v) == UnitValue :
            if self.units.dim == v.units.dim :
                return (self.value >= v.convert(self.units.sys).value)
            else :
                raise ValueError("cannot compare UnitValues with different units dimensions.")
        elif isnumber(v):
            return self.value >= v
        else :
            return TypeError("UnitValue can only be compared with numbers and UnitValues with the same units dimensions.")


    def __lt__(self, v) :
        """
        defines self < v or v > self.

        if v is not a UnitValue, self.value < v is returned.
        if v is a UnitValue, numerical value converted to the same units system
        are compared. A ValueError is raised if self and v have different units dimensions.
        a TypeError is raised if v is not a UnitValue or a number.

        ie.

        .. code:: python

            UnitValue(1, "µm") < UnitValue(1, "µm")   # False
            UnitValue(1, "µm") < UnitValue(1, "mm")   # False
            UnitValue(1, "µm") < 0.5                  # False
            UnitValue(1, "µm") < UnitValue(1, "µm/s") # ValueError

        """

        if type(v) == UnitValue :
            if self.units.dim == v.units.dim :
                return (self.value < v.convert(self.units.sys).value)
            else :
                raise ValueError("cannot compare UnitValues with different units dimensions.")
        elif isnumber(v):
            return self.value < v
        else :
            return TypeError("UnitValue can only be compared with numbers and UnitValues with the same units dimensions.")

    def __le__(self, v) :
        """
        defines self <= v or v >= self.

        if v is not a UnitValue, self.value <= v is returned.
        if v is a UnitValue, numerical value converted to the same units system
        are compared. A ValueError is raised if self and v have different units dimensions.
        a TypeError is raised if v is not a UnitValue or a number.

        ie.

        .. code:: python

            UnitValue(1, "µm") <= UnitValue(1, "µm")   # True
            UnitValue(1, "µm") <= UnitValue(1, "mm")   # False
            UnitValue(1, "µm") <= 0.5                  # False
            UnitValue(1, "µm") <= UnitValue(1, "µm/s") # ValueError

        """

        if type(v) == UnitValue :
            if self.units.dim == v.units.dim :
                return (self.value <= v.convert(self.units.sys).value)
            else :
                raise ValueError("cannot compare UnitValues with different units dimensions.")
        elif isnumber(v):
            return self.value <= v
        else :
            return TypeError("UnitValue can only be compared with numbers and UnitValues with the same units dimensions.")


    def copy(self) :
        """
        Returns a deepcopy of the instance.

        .. code:: python

            instance.copy()

            # is equivalent to
            # import copy

            copy.deepcopy(instance)

        """

        return cpy.deepcopy(self)

def parse_unitvalue(s="") :
    """
    Create a UnitValue from the character string s.
    the value and the units must be separated by one or more whitespace :
    ie. "1e3 µmol/L.s-1" ok
    "2m" wrong
    """
    if type(s) != str :
        raise "s must be a string."
    s = s.strip()

    tok = s.split()
    if len(tok) == 0 :
        value = 0
        units = parse_units("")
    else :
        value = float(tok[0])
        us = ""
        for i in range(1, len(tok)):
            us += tok[i]
        units = parse_units(us)
    return UnitValue(value, units)

class UnitArray :
    """
    Array of values with the same units.
    this is not a conventionnal array per se, as it does not expose a proper arrya interface.
    instread, it should rather be seen as an array wrapper, proposed for conveinience, especially for
    bulk unit conversions.

    :param value: array of values
    :type value: array of number and/or UnitArray with the same unit dimensions
    :param units: units of the variable
    :type units: Units or str
    """

    def __init__(self, value, units=None, check_value=True, convert=True) :
        """
        constructor
        """

        if type(units) == str :
            units = Units(units)
        elif type(units) == Units :
            units = units.copy()
        elif isnone(units) : 
            pass
        else :
            raise TypeError("units must be a str, an instance of the Units class or None.")

        if isarray(value) : 
            if isnone(units) : # no units
                self.units = Units("")
            else :
                self.units = units
            
            self.set_value(value, check=check_value)  
         
        elif type(value) == UnitArray :
            if isnone (units): # no units
                self.units = value.units
                self.set_value(value.value, check=False)
            else :
                if convert : #units compatibility check is done in the conversion
                    value = value.convert(units)
                else : #only the units compatibility check must be done
                    if value.units.dim != units.dim : 
                        raise ValueError("the value "+str(value)+" and the units " + str(units)+ " have incompatble dimensions.")
                self.units = value.units
                self.set_value(value.value, check=False)
        else :
            raise TypeError("unsupported value type.")
            
        #type check :
        # self.units = units #units must be set first
        # self.set_value(value, check=check_value)

    @property
    def value(self) :
        """
        an array of numerical values, all associated with the same units.
        (numpy.ndarray [#numpy_ndarray]_ of number).
        Setting the attribute is equivalent to calling self.set_value with check=True.
        You may reffer to the latter for more information on accepted values for the property.
        """
        # references :
        # .. [#numpy_ndarray] Numpy Developers. numpy 1.25 API reference : numpy.ndarray. (consulted on september 04, 2023). https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html

        return self._value

    @value.setter
    def value(self, v) :
        self.set_value(v)

    @property
    def units(self) :
        """
        Units of the physical values (:py:class:`Units`).
        Can be set by a string or a :py:class:`Units` instance.
        """

        return self._units

    @units.setter
    def units(self, v) :
        if type(v) == str :
            v = parse_units(v)
        if type(v) == Units :
            self._units = v.copy()
        else :
            raise ValueError("UnitValue's unit must be an instance of Units or a unit string.")

    def __str__(self) :
        """
        str cast, mainly for printing purposes
        """

        return str(self.value) + " " + str(self.units)

    def set_value(self, v, check=True) :
        """
        Sets the value property of the UnitArray.

        :param v: new value for the value property. it must be an array.
            If check=True, array elements can be string and/or numbers and/or UnitValues.

            * numbers will be accepted as such.
            * for UnitValue objects, the value is taken after proper conversion to match self.units.
                a ValueError is raise if the UnitValue's units dimensions does not match self.units.dim.
            * strings are parsed as UnitValues then processed as such.

            However, if check=False, all elements should be numbers, as any element that is not a
            number may result in unexpected behaviours.

        :type v: array of numbers and/or UnitValues and/or str.
        :param check: tells if v should be checked or not. Can be set to False if one is sure that v is an array
            of numbers, without UnitValues or strings.
        :type check: bool
        """

        if not check :
            self._value = np.array(v, dtype=float)
        else :
            if isarray(v) :
                self._value = np.array(v)
                for i in range(len(self._value)) :
                    if type(self._value[i]) == str :
                        self._value[i] = parse_unitvalue(self._value[i])

                    if type(self._value[i]) == UnitValue :
                        if self._value[i].units.dim != self.units.dim :
                            raise ValueError("units dimensions of item "+str(i)+" does not match the UnitArray units.")

                        if self._value[i].units.sys != self.units.sys :
                            self._value[i] = self._value[i].convert(self.units)

                        self._value[i] = self._value[i].value
                    elif isnumber(self._value[i]) :
                        pass
                self._value = np.array(self._value, dtype=float)
            else :
                raise ValueError("UnitArray's value must be an array.")

    def convert(self, u) :
        """
        Returns the variable converted the units u.
        Same as convert_value(v, u), wher v is the UnitValue instance.
        """

        if type(u) == str :
            u =  parse_units(u)

        su_dst = 0
        
        if type(u) == UnitValue :
            if u.units.dim != self.units.dim :
                raise ValueError("unit conversion must happen in the same dimension.")
            su_dst = u.units.sys
        elif type(u) == Units :
            if u.dim != self.units.dim :
                raise ValueError("unit conversion must happen in the same dimension.")
            su_dst = u.sys
        elif type(u) == UnitsSystem :
            su_dst = u
        elif type(u) == dict :
            su_dst = unitssystem_from_dict(u)
        else :
            raise TypeError("conversion accepts strings, UnitValues, Units, UnitsSystems and dictionnaries only.")
        
        return UnitArray(
            convert_value(self.value, self.units.sys, su_dst, self.units.dim), 
            Units(su_dst, self.units.dim), 
            check_value=True
            )

    def get_at(self, i) :
        """
        Return the element of self.value at i as a UnitValue.
        """

        return UnitValue(self.value[i], self.units)

    def set_at(self, i, v) :
        """
        Sets the element of v.value at i as a the UnitValue v.
        units dimensions must be compatible.
        """

        if type(v) == UnitValue :
            self.value[i] = v.convert(self.units).value
        else :
            self.value[i] = v.convert(self.units).v

    def __len__(self) :
        return len(self.value)

    def _sum(self, v) :
        """
        returns self + v.
        """

        if type(v)==UnitValue :
            if self.units.dim != v.units.dim :
                raise ValueError("addition of a UnitArray with a UnitValue with different units dimensions is not supported.")
            v = v.convert(self.units.sys)
            return UnitArray([self.value[i] + v.value for i in range(len(self))], self.units)

        elif isnumber(v) :
            return UnitArray([self.value[i] + v for i in range(len(self))], self.units)

        elif type(v)==UnitArray :
            if self.units.dim != v.units.dim :
                raise ValueError("addition of a UnitArray with a UnitValue with different units dimensions is not supported.")
            if len(self) != len(v) :
                raise ValueError("operations between UnitArrays require both the have the same length.")
            v = v.convert(self.units.sys)
            return UnitArray([self.value[i] + v.value[i] for i in range(len(self))], self.units)

        else :
            raise ValueError("unexpected types for a sum.")

    def _product(self, v) :
        """
        returns self * v.
        """

        if type(v)==UnitValue :
            v = v.convert(self.units.sys)
            return UnitArray([self.value[i] * v.value for i in range(len(self))], self.units.multiply(v.units))

        elif isnumber(v) :
            return UnitArray([self.value[i] * v for i in range(len(self))], self.units)

        elif type(v)==UnitArray :
            if len(self) != len(v) :
                raise ValueError("operations between UnitArrays require both the have the same length.")
            v = v.convert(self.units.sys)
            return UnitArray([self.value[i] * v.value[i] for i in range(len(self))], self.units.multiply(v.units))

        else :
            raise ValueError("unexpected types for a product.")

    def _modulo(self, mod):
        """
        returns self % mod.
        """

        if type(mod)==UnitValue :
            if self.units.dim != mod.units.dim :
                raise ValueError("modulus mod must have the same units dimensions as self.")
            mod = mod.convert(self.units.sys)
            return UnitArray([self.value[i] % mod.value for i in range(len(self))], self.units)

        elif isnumber(mod) :
            return UnitArray([self.value[i] % mod for i in range(len(self))], self.units)

        elif type(mod)==UnitArray :
            if self.units.dim != mod.units.dim :
                raise ValueError("modulus mod must have the same units dimensions as self.")
            if len(self) != len(mod) :
                raise ValueError("operations between UnitArrays require both the have the same length.")
            mod = mod.convert(self.units.sys)
            return UnitArray([self.value[i] % mod.value[i] for i in range(len(self))], self.units)

        else :
            raise ValueError("unexpected types for a modulo.")

    def _rmodulo(self, v):
        """
        returns v % self.
        """

        if type(v)==UnitValue :
            if self.units.dim != v.units.dim :
                raise ValueError("v must have the same units dimensions as the modulus self.")
            v = v.convert(self.units.sys)
            return UnitArray([v.value % self.value[i] for i in range(len(self))], self.units)

        elif isnumber(v) :
            return UnitArray([v % self.value[i] for i in range(len(self))], self.units)

        elif type(v)==UnitArray :
            if self.units.dim != v.units.dim :
                raise ValueError("modulo operation on two UnitArrays requires both to have the same units dimensions.")
            if len(self) != len(v) :
                raise ValueError("operations between UnitArrays require both the have the same length.")
            v = v.convert(self.units.sys)
            return UnitArray([v.value[i] % self.value[i] for i in range(len(self))], self.units)

        else :
            raise ValueError("unexpected types for a modulo.")

    def invert(self) :
        """
        Invert the variable (returns 1/v).

        ie.

        .. code:: python

            UnitArray([1,2], "m").invert() # and
            UnitArray([1,0.5], "m-1")      # are equivalent

        """

        return UnitArray(1/self.value, self.units.invert())

    def __add__    (self, v) :
        """
        defines self+v.
        Addition must be operated betwen UnitValue/UnitArrays with the same units dimension.
        if v is not a UnitValue nor a UnitArray, is is conseidered to be in the same units than self.
        Conversion is implicitely done if required.

        ie.

        .. code:: python

            UnitArray([1,2], "m") + 1                     # ok, == UnitArray([3,4], "m")
            UnitArray([1,2], "m") + UnitValue(1000, "mm") # ok, == UnitArray([3,4], "m")
            UnitArray([1,2], "m") + UnitArray([3,4], "m") # ok, == UnitArray([4,6], "m")

            UnitArray([1,2], "m") + UnitValue(1, "m/s")   # ValueError

        """

        return self._sum(v)

    def __radd__    (self, v) :
        """
        defines v+self. same as self.__add__(v).
        """

        return self._sum(v)

    def __sub__    (self, v) :
        """
        defines self-v. same as self.__add__(-v).
        """

        return self._sum(_neg(v))

    def __rsub__    (self, v) :
        """
        defines v-self. sames as -self.__sub__(v).
        """

        return _neg(self._sum(_neg(v)))

    def __mul__    (self, v) :
        """
        defines self*v.
        if v is not a UnitValue no a UnitArray, is is conseidered to be a unitless factor.
        Conversion is implicitely done if required.

        ie.

        .. code:: python

            UnitArray([1,2], "m") * UnitValue(1, "m-1")   # ok, == UnitArray([1,2], "")
            UnitArray([1,2], "m") * 2                     # ok, == UnitArray([2,4], "m")
            UnitArray([1,2], "m") * UnitValue(1, "mol/s") # ok, == UnitArray([1,2], "m.mol/s")
            UnitArray([1,2], "m") * UnitArray([1,2], "s") # ok, == UnitArray([1,4], "m.s")

        """

        return self._product(v)

    def __rmul__    (self, v) :
        """
        defines v*self. same as self.__mul__(v).
        """

        return self._product(v)

    def __truediv__(self, v) :
        """
        defines self/v.
        if v is not a UnitValue nor a UnitArray, is is conseidered to be a unitless factor.
        Conversion is implicitely done if required.

        ie.

        .. code:: python

            UnitArray([1,2], "m") / 2                     # ok, == UnitArray([0.5,1], "m")
            UnitArray([1,2], "m") / UnitValue("2 m")      # ok, == UnitArray([0.5,1], "")
            UnitArray([1,2], "m") / UnitArray([1,2], "s") # ok, == UnitArray([1,1], "m/s")

        """

        return self._product(_inv(v))

    def __rtruediv__(self, v) :
        """
        defines v/self. ...
        """

        return self.invert()._product(v)

    def __mod__    (self, v) :
        """
        defines self%v.
        if v is not a UnitValue, self.value%v is returned.
        if v is a UnitValue, v and self must have the same units dimensions, a ValueError is returned otherwise.
        The operation is applied to the numeric values after conversion in the same units system.

        """

        return self._modulo(v)

    def __rmod__    (self, v) :
        """
        defines v%self.
        if v is not a UnitValue, v%self.value is returned.
        if v is a UnitValue, v and self must have the same units dimensions, a ValueError is returned otherwise.
        The operation is applied to the numeric values after conversion in the same units system.

        """

        return self._rmodulo(v)

    def __pow__    (self, v) :
        """
        defines self**v.
        pow is valid as long as the resluting units dimensions are integers.

        ie.

        .. code: python

            Units(1, "m")**3      # ok, 3*1=3 is still an integer
            Units(1, "m3")**(1/3) # ok, 3/3=1 is still an integer
            Units(1, "m")**(1/3)  # ValueError, 1/3 is not an integer
            Units(1, "")**(1/3)   # ok, 0/3=0 is still an integer

        """
        raise NotImplementedError()

        # if type(v) != UnitValue :
        #     return UnitValue(self.value**v, self.units.raiseto(v))
        # else:
        #     raise ValueError("UnitValue**UnitValue is invalid. UnitValue**int is.")

    def __rpow__    (self, v) :
        """
        defines v**self is an, invalid opeartion.
        raising a number or a UnitValue to the power  of a UnitValue is not a supported opeartion.
        raise a ValueError.

        """

        raise ValueError("UnitValue**UnitValue is invalid. UnitValue**int is.")

    def __neg__(self) :
        """
        defines the unary operation -self.
        returns a copy of self, but with an opposed value.

        """

        return UnitArray(-self.value, self.units)

    def __pos__(self) :
        """
        defines the unary operation +self.
        returns a copy of self.

        """

        return self.copy()

    def __abs__(self) :
        """
        defines abs(self).
        simply returns UnitValue(abs(self.value), self.units).

        """

        return UnitArray(abs(self.value), self.units)

    def copy(self) :
        """
        Returns a deepcopy of the instance.

        .. code:: python

            instance.copy()

            # is equivalent to
            # import copy

            copy.deepcopy(instance)

        """

        return cpy.deepcopy(self)

def unitarray_to_dict(v) :
    """
    Returns a dict representing the UnitArray v.

    :param v: unit array to be represented as a dict
    :type v: UnitArray
    :returns: dict representing v.
    :rtype: dict
    """

    return {"value" : v.value.tolist(), "units" : str(v.units)}

def unitarray_from_dict(d, base_path=None) :
    """
    Returns a the UnitArray represented by the dict d.

    :param v: dict represented a UnitArray
    :type v: dict
    :returns: UintArray represented by d.
    :rtype: UnitArray
    """

    d = valproc.process_input_dict_keys(d, [
                ["value"],
                ["units"]
            ]
        )
    if isstr(d["value"]) :
        path = filepath.get_path_with_base(d["value"], base_path)
        value = np.load(path)
        return UnitArray(value, d["units"])
    else :
        return UnitArray(d["value"], d["units"])

def _neg(v) :
    """
    returns -v.
    if v is a list, this is done element wise.
    """

    if type(v) == list :
        return [_neg(vi) for vi in v]
    else :
        return -v

def _inv(v) :
    """
    returns 1/v.
    if v is a list, this is done element wise.
    for UnitsArray or UnitsValue, the invert method is called.
    """

    if type(v) == list :
        return [_inv(vi) for vi in v]
    if type(v) == UnitValue or type(v) == UnitArray:
        return v.invert()
    else :
        return 1/v
