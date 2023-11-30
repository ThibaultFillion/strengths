import sys
import numpy
import ctypes
sys.path.append("../src/")
from strengths import *
from strengths.librdengine import *

# tests ###############################################


def test_get_sample_index() :
    script = RDScript(
        RDSystem(
            RDNetwork(
                species=[Species("A")],
                reactions=[]
                )
            ),
        [1,2,3,4,5]
        )
    out = RDTrajectory(UnitArray([1,2,3,4], "molecule"), UnitArray([1, 2, 3, 10], "s"), system=script.system, script=script, engine_description="", engine_option="")
    
    # same units, policy=="closest"
    
    assert out.get_sample_index("-1 s")                      == 0
    assert out.get_sample_index("0 s")                       == 0
    assert out.get_sample_index("1 s")                       == 0
    assert out.get_sample_index("1.2 s")                     == 0
    assert out.get_sample_index("1.5 s")                     == 0
    assert out.get_sample_index("1.7 s")                     == 1
    assert out.get_sample_index("2 s")                       == 1
    assert out.get_sample_index("2.41 s")                    == 1
    assert out.get_sample_index("2.49999 s")                 == 1
    assert out.get_sample_index("2.499999999999999999999 s") == 1
    assert out.get_sample_index("2.51 s")                    == 2
    assert out.get_sample_index("2.6 s")                     == 2
    assert out.get_sample_index("3 s")                       == 2
    assert out.get_sample_index("4 s")                       == 2
    assert out.get_sample_index("10.01 s")                   == 3
    assert out.get_sample_index("10000000000000 s")          == 3
    
    # same units, policy=="supeq"
    
    assert out.get_sample_index("-1 s", policy="supeq")                      == 0
    assert out.get_sample_index("0 s", policy="supeq")                       == 0
    assert out.get_sample_index("1 s", policy="supeq")                       == 0
    assert out.get_sample_index("1.2 s", policy="supeq")                     == 1
    assert out.get_sample_index("1.5 s", policy="supeq")                     == 1
    assert out.get_sample_index("1.7 s", policy="supeq")                     == 1
    assert out.get_sample_index("2 s", policy="supeq")                       == 1
    assert out.get_sample_index("2.41 s", policy="supeq")                    == 2
    assert out.get_sample_index("2.49999 s", policy="supeq")                 == 2
    assert out.get_sample_index("2.499999999999999999999 s", policy="supeq") == 2
    assert out.get_sample_index("2.51 s", policy="supeq")                    == 2
    assert out.get_sample_index("2.6 s", policy="supeq")                     == 2
    assert out.get_sample_index("3 s", policy="supeq")                       == 2
    assert out.get_sample_index("4 s", policy="supeq")                       == 3
    assert out.get_sample_index("10.01 s", policy="supeq")                   is None
    assert out.get_sample_index("10000000000000 s", policy="supeq")          is None
    
    # same units, policy=="infeq"
    
    assert out.get_sample_index("-1 s", policy="infeq")                      is None
    assert out.get_sample_index("0 s", policy="infeq")                       is None
    assert out.get_sample_index("1 s", policy="infeq")                       == 0
    assert out.get_sample_index("1.2 s", policy="infeq")                     == 0
    assert out.get_sample_index("1.5 s", policy="infeq")                     == 0
    assert out.get_sample_index("1.7 s", policy="infeq")                     == 0
    assert out.get_sample_index("2 s", policy="infeq")                       == 1
    assert out.get_sample_index("2.41 s", policy="infeq")                    == 1
    assert out.get_sample_index("2.49999 s", policy="infeq")                 == 1
    assert out.get_sample_index("2.499999999999999999999 s", policy="infeq") == 1
    assert out.get_sample_index("2.51 s", policy="infeq")                    == 1
    assert out.get_sample_index("2.6 s", policy="infeq")                     == 1
    assert out.get_sample_index("3 s", policy="infeq")                       == 2
    assert out.get_sample_index("4 s", policy="infeq")                       == 2
    assert out.get_sample_index("10.01 s", policy="infeq")                   == 3
    assert out.get_sample_index("10000000000000 s", policy="infeq")          == 3
        
    
    
    # different units, policy=="closest"
    
    assert out.get_sample_index("-1000 ms")                   == 0
    assert out.get_sample_index("0000 ms")                    == 0
    assert out.get_sample_index("1000 ms")                    == 0
    assert out.get_sample_index("1200 ms")                    == 0
    assert out.get_sample_index("1500 ms")                    == 0
    assert out.get_sample_index("1700 ms")                    == 1
    assert out.get_sample_index("2000 ms")                    == 1
    assert out.get_sample_index("2410 ms")                    == 1
    assert out.get_sample_index("2499.99 ms")                 == 1
    assert out.get_sample_index("2499.999999999999999999 ms") == 1
    assert out.get_sample_index("2510 ms")                    == 2
    assert out.get_sample_index("2600 ms")                    == 2
    assert out.get_sample_index("3000 ms")                    == 2
    assert out.get_sample_index("4000 ms")                    == 2
    assert out.get_sample_index("10010 ms")                   == 3
    assert out.get_sample_index("10000000000000000 ms")       == 3

    # different units, policy=="supeq"
    
    assert out.get_sample_index("-1000 ms", policy="supeq")                   == 0
    assert out.get_sample_index("0000 ms", policy="supeq")                    == 0
    assert out.get_sample_index("1000 ms", policy="supeq")                    == 0
    assert out.get_sample_index("1200 ms", policy="supeq")                    == 1
    assert out.get_sample_index("1500 ms", policy="supeq")                    == 1
    assert out.get_sample_index("1700 ms", policy="supeq")                    == 1
    assert out.get_sample_index("2000 ms", policy="supeq")                    == 1
    assert out.get_sample_index("2410 ms", policy="supeq")                    == 2
    assert out.get_sample_index("2499.99 ms", policy="supeq")                 == 2
    assert out.get_sample_index("2499.999999999999999999 ms", policy="supeq") == 2
    assert out.get_sample_index("2510 ms", policy="supeq")                    == 2
    assert out.get_sample_index("2600 ms", policy="supeq")                    == 2
    assert out.get_sample_index("3000 ms", policy="supeq")                    == 2
    assert out.get_sample_index("4000 ms", policy="supeq")                    == 3
    assert out.get_sample_index("10010 ms", policy="supeq")                   is None
    assert out.get_sample_index("10000000000000000 ms", policy="supeq")       is None
        
    # different units, policy=="infeq"
    
    assert out.get_sample_index("-1000 ms", policy="infeq")                   is None
    assert out.get_sample_index("0000 ms", policy="infeq")                    is None
    assert out.get_sample_index("1000 ms", policy="infeq")                    == 0
    assert out.get_sample_index("1200 ms", policy="infeq")                    == 0
    assert out.get_sample_index("1500 ms", policy="infeq")                    == 0
    assert out.get_sample_index("1700 ms", policy="infeq")                    == 0
    assert out.get_sample_index("2000 ms", policy="infeq")                    == 1
    assert out.get_sample_index("2410 ms", policy="infeq")                    == 1
    assert out.get_sample_index("2499.99 ms", policy="infeq")                 == 1
    assert out.get_sample_index("2499.999999999999999999 ms", policy="infeq") == 1
    assert out.get_sample_index("2510 ms", policy="infeq")                    == 1
    assert out.get_sample_index("2600 ms", policy="infeq")                    == 1
    assert out.get_sample_index("3000 ms", policy="infeq")                    == 2
    assert out.get_sample_index("4000 ms", policy="infeq")                    == 2
    assert out.get_sample_index("10010 ms", policy="infeq")                   == 3
    assert out.get_sample_index("10000000000000000 ms", policy="infeq")       == 3

# run all ###############################################

def run_all_tests() :
    test_get_sample_index()