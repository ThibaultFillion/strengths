import numbers
import collections
import numpy as np

def isnone(x) :
    """
    Returns True if x is None, False otherwise.
    """

    return x is None

def isnumber(x):
    """
    Returns True if x is a number, False otherwise.
    """

    return isinstance(x, numbers.Number)

def isarray(x):
    """
    Returns True if x is an array (numpy.ndarray [#numpy_ndarray]_ or list or tuple), False otherwise.
    """
    # references :
    # .. [#numpy_ndarray] Numpy Developers. numpy 1.25 API reference : numpy.ndarray. (consulted on september 04, 2023). https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html
    
    return type(x) == np.ndarray or type(x) == list or type(x) == tuple

def isdict(x) :
    """
    Returns True if x is a dict, False otherwise.
    """

    return type(x) == dict

def isstr(x) :
    """
    Returns True if x is a str, False otherwise.
    """

    return type(x) == str

def array_to_list(x) :
    if type(x) == np.ndarray : 
        return x.tolist()
    else :
        return list(x)
