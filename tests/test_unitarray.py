import sys
sys.path.append("../src/")
from strengths.units import *

def test_constructor() :
    a = UnitArray([1,2,3,4], "µm")
    b = UnitArray([1, "1 m", UnitValue(1, "nm"), 4], "µm")

def test_convert() :
    a = UnitArray([1,2,3,4], "µm")
    b = a.convert("m")
    assert list(b.value) == [1e-6, 2e-6, 3e-6, 4e-6]
