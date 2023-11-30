import pathlib 
from strengths.typechecking import *

def have_extension(path, ext) :
    """
    returns True if path is ends with ext, False otherwise.
    """
    
    return (path[len(path)-len(ext):len(path)] == ext)

def append_extension_if_missing(path, ext) : 
    """
    add ext to path if path doesnt already ends with ext
    """
    
    if have_extension(path, ext) : 
        return path
    else :
        return path+ext

def remove_extension_if_existing(path, ext) : 
    """
    add ext to path if path doesnt already ends with ext
    """
    
    if have_extension(path, ext) : 
        return path[0:len(path)-len(ext)]
    else :
        return path
    
def get_base_path(path) : 
    """
    returns the subdirectory path
    """
    
    p = pathlib.Path(path)
    if not p.is_absolute():
        p = p.absolute()
    return str(p.parent)
    
def get_path_with_base(path, base_path) : 
    """
    returns a path
    """
    
    if isnone(base_path) :
        return path
    elif isstr(base_path) :
        if pathlib.Path(path).is_absolute() :
            return path
        else :
            return str(pathlib.Path(base_path).joinpath(path))
    else :
        raise TypeError("base_path must be a string (path) or None.")

def get_last_element(path) : 
    """
    returns the final element of the path
    """
    
    return pathlib.Path(path).name

