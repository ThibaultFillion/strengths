import sys
sys.path.append("../src/")
from strengths import *

# tests ###############################################

def test_parse_units() : 
    q = parse_units("")
    assert str(q) == ""
    q = parse_units("µM")
    assert str(q) == "dm-3.µmol"
    q = parse_units("m/s")
    assert str(q) == "m.s-1"
    q = parse_units("m.s.mol0")
    assert str(q) == "m.s"
    q = parse_units("   m.s  ")
    assert str(q) == "m.s"

def test_parse_unitvalue() :
    q = parse_unitvalue("1 µs")
    assert str(q) == "1.0 µs"
    q = parse_unitvalue("1 ")
    assert str(q) == "1.0 "
    q = parse_unitvalue("1   µs")
    assert str(q) == "1.0 µs"
    q = parse_unitvalue("  1 µs  ")
    assert str(q) == "1.0 µs"

def test_unitvalue_conversion() :
    q1 = parse_unitvalue("1 M")
    q2 = q1.convert("µM")
    assert q2.value == 1e6
    q1 = parse_unitvalue("1 m2/s")
    q2 = q1.convert("m2/h")
    assert q2.value == 1*60*60    
    q1 = parse_unitvalue("1 mol")
    q2 = q1.convert("molecule")
    assert q2.value == constants.avogadro_number()

def test_op_addition_uv_number() :
    q1 = parse_unitvalue("1 mol")
    q2 = q1 + 1
    assert q2 == parse_unitvalue("2 mol")

def test_op_substraction_uv_number() :
    q1 = parse_unitvalue("1 mol")
    q2 = q1 - 1
    assert q2 == parse_unitvalue("0 mol")
    assert parse_unitvalue("1 m") - 0.5 == parse_unitvalue("0.5 m")
    
def test_op_multiplication_uv_number() :
    q1 = parse_unitvalue("1 mol")
    q2 = q1 * 2
    assert q2 == parse_unitvalue("2 mol")

def test_op_addition_number_uv() :
    assert 0.5 + parse_unitvalue("1 m") == parse_unitvalue("1.5 m")
    
def test_op_substraction_number_uv() :
    assert 0.5 - parse_unitvalue("1 m") == parse_unitvalue("-0.5 m")

def test_op_multiplication_number_uv() :
    assert 2 * parse_unitvalue("1 m") == parse_unitvalue("2 m")
    
def test_op_division_number_uv() :
    assert 2 / parse_unitvalue("1 m") == parse_unitvalue("2 m-1")
    assert 1 / parse_unitvalue("2 m") == parse_unitvalue("0.5 m-1")

def test_op_addition_uv_uv() :
    q1 = parse_unitvalue("1 s")
    q2 = parse_unitvalue("1 min")
    assert q1+q2 == parse_unitvalue("61 s")

def test_op_division_uv_uv() :
    assert parse_unitvalue("1 m2")/parse_unitvalue("1 m-3") == parse_unitvalue("1 m5")

def test_op_division_uv_number() :
    q1 = parse_unitvalue("1 mol")
    q2 = q1 / 2
    assert q2 == parse_unitvalue("0.5 mol")
     
def test_op_multiplication_uv_uv() :
    q1 = parse_unitvalue("1 mol")
    q2 = parse_unitvalue("1 s")
    assert q1/q2 == parse_unitvalue("1 mol/s")
    q1 = parse_unitvalue("1 molecule")
    q2 = parse_unitvalue("1 s")
    q3 = parse_unitvalue("1 m")
    assert (q1*q2*q3).convert("molecule.s.m") == parse_unitvalue("1.0 molecule.s.m")

def test_op_power_uv_number() :
    q1 = parse_unitvalue("1 s")
    q2 = parse_unitvalue("1 min")
    assert q1+q2 == parse_unitvalue("61 s")
    assert parse_unitvalue("1 m")**2 == parse_unitvalue("1 m2")
    assert parse_unitvalue("1 m")**2/parse_unitvalue("1 s") == parse_unitvalue("1 m2/s")    
    assert parse_unitvalue("1 m3")**(1/3) == parse_unitvalue("1 m")

def test_op_modulo_uv_number() :
    assert UnitValue("3 m") % 2 == UnitValue("1 m")

def test_op_modulo_number_uv() :
    assert 3 % UnitValue("2 m") == UnitValue("1 m")    

def test_op_modulo_uv_uv() :
    assert UnitValue("3 m") % UnitValue("2 m") == UnitValue("1 m")

def test_op_modulo_uv_ua() : 
    a = UnitValue("3 m") % UnitArray([2, 3], "m")
    assert (a.get_at(0)==UnitValue("1 m") and
            a.get_at(1)==UnitValue("0 m") )

def test_op_neg_ua() :
    a = -UnitArray([1,2,3], "µm")
    assert (a.get_at(0)==UnitValue("-1 µm") and
            a.get_at(1)==UnitValue("-2 µm") and
            a.get_at(2)==UnitValue("-3 µm") )

def test_ua_invert() :
    a = UnitArray([1,2,3], "µm").invert()
    assert (a.get_at(0)==UnitValue("1 µm-1") and
            a.get_at(1)==UnitValue("0.5 µm-1") and
            a.get_at(2)==UnitValue("0.33333333333333333 µm-1") )

def test_op_multiplication_uv_ua() :
    a = parse_unitvalue("1 µm") * UnitArray([1,2,3], "")
    assert (a.get_at(0)==UnitValue("1 µm") and
            a.get_at(1)==UnitValue("2 µm") and
            a.get_at(2)==UnitValue("3 µm") )

def test_op_division_uv_ua() :
    a = parse_unitvalue("1 µm") / UnitArray([1,2,3], "")
    assert (a.get_at(0)==UnitValue("1 µm") and
            a.get_at(1)==UnitValue("0.5 µm") and
            a.get_at(2)==UnitValue("0.33333333333333333333 µm") )

def test_op_addition_uv_ua() :
    a = parse_unitvalue("1 µm") + UnitArray([1,2,3], "µm")
    assert (a.get_at(0)==UnitValue("2 µm") and
            a.get_at(1)==UnitValue("3 µm") and
            a.get_at(2)==UnitValue("4 µm") )

def test_op_substraction_uv_ua() :
    a = parse_unitvalue("1 µm") - UnitArray([1,2,3], "µm")
    assert (a.get_at(0)==UnitValue("0 µm") and
            a.get_at(1)==UnitValue("-1 µm") and
            a.get_at(2)==UnitValue("-2 µm") )

def test_op_multiplication_ua_uv() :
    a = UnitArray([1,2,3], "") * UnitValue("1 µm")
    assert (a.get_at(0)==UnitValue("1 µm") and
            a.get_at(1)==UnitValue("2 µm") and
            a.get_at(2)==UnitValue("3 µm") )

def test_op_division_ua_uv() :
    a = UnitArray([1,2,3], "") / parse_unitvalue("2 µm")
    assert (a.get_at(0)==UnitValue("0.5 µm-1") and
            a.get_at(1)==UnitValue("1 µm-1") and
            a.get_at(2)==UnitValue("1.5 µm-1") )

def test_op_addition_ua_uv() :
    a = UnitArray([1,2,3], "µm") + parse_unitvalue("1 µm")
    assert (a.get_at(0)==UnitValue("2 µm") and
            a.get_at(1)==UnitValue("3 µm") and
            a.get_at(2)==UnitValue("4 µm") )

def test_op_substraction_ua_uv() :
    a = UnitArray([1,2,3], "µm") - parse_unitvalue("1 µm")
    assert (a.get_at(0)==UnitValue("0 µm") and
            a.get_at(1)==UnitValue("1 µm") and
            a.get_at(2)==UnitValue("2 µm") )

def test_op_modulo_ua_number() :
    a = UnitArray([1,2,3], "µm") % 2
    assert (a.get_at(0)==UnitValue("1 µm") and
            a.get_at(1)==UnitValue("0 µm") and
            a.get_at(2)==UnitValue("1 µm") )

def test_op_modulo_ua_uv() :
    a = UnitArray([1,2,3], "µm") % UnitValue("2 µm")
    assert (a.get_at(0)==UnitValue("1 µm") and
            a.get_at(1)==UnitValue("0 µm") and
            a.get_at(2)==UnitValue("1 µm") )

def test_op_modulo_ua_ua() :
    a = UnitArray([1,2,3], "µm") % UnitArray([1,2,3], "µm")
    assert (a.get_at(0)==UnitValue("0 µm") and
            a.get_at(1)==UnitValue("0 µm") and
            a.get_at(2)==UnitValue("0 µm") )

def test_op_addition_ua_number() :
    a = UnitArray([1,2,3], "µm") + 2
    assert (a.get_at(0)==UnitValue("3 µm") and
            a.get_at(1)==UnitValue("4 µm") and
            a.get_at(2)==UnitValue("5 µm") )

def test_op_addition_number_ua() :
    a = 2 + UnitArray([1,2,3], "µm")
    assert (a.get_at(0)==UnitValue("3 µm") and
            a.get_at(1)==UnitValue("4 µm") and
            a.get_at(2)==UnitValue("5 µm") )

def test_op_substraction_ua_number() :
    a = UnitArray([1,2,3], "µm") - 2
    assert (a.get_at(0)==UnitValue("-1 µm") and
            a.get_at(1)==UnitValue("0 µm") and
            a.get_at(2)==UnitValue("1 µm") )

def test_op_substraction_number_ua() :
    a = 2 - UnitArray([1,2,3], "µm")
    assert (a.get_at(0)==UnitValue("1 µm") and
            a.get_at(1)==UnitValue("0 µm") and
            a.get_at(2)==UnitValue("-1 µm") )

def test_op_multiplication_number_ua() :
    a = 2 * UnitArray([1,2,3], "µm")
    assert (a.get_at(0)==UnitValue("2 µm") and
            a.get_at(1)==UnitValue("4 µm") and
            a.get_at(2)==UnitValue("6 µm") )

def test_op_multiplication_ua_number() :
    a = UnitArray([1,2,3], "µm") * 2
    assert (a.get_at(0)==UnitValue("2 µm") and
            a.get_at(1)==UnitValue("4 µm") and
            a.get_at(2)==UnitValue("6 µm") )

def test_op_division_number_ua() :
    a = 2 / UnitArray([1,2,3], "µm")
    assert (a.get_at(0)==UnitValue("2 µm-1") and
            a.get_at(1)==UnitValue("1 µm-1") and
            a.get_at(2)==UnitValue("0.6666666666666666 µm-1") )

def test_op_division_ua_number() :
    a = UnitArray([1,2,3], "µm") / 2
    assert (a.get_at(0)==UnitValue("0.5 µm") and
            a.get_at(1)==UnitValue("1 µm") and
            a.get_at(2)==UnitValue("1.5 µm") )

def test_op_abs_uv() :
    assert abs(UnitValue("-2 µm")) == UnitValue("2 µm")
    assert abs(UnitValue("2 µm")) == UnitValue("2 µm")

def test_op_abs_ua() :
    a = abs(UnitArray([1,-2,3], "µm"))
    assert (a.get_at(0)==UnitValue("1 µm") and
            a.get_at(1)==UnitValue("2 µm") and
            a.get_at(2)==UnitValue("3 µm") )

def test_op_equals_uv() :
    assert UnitValue("1 µm") == 1
    assert 1 == UnitValue("1 µm")
    assert UnitValue("1 µm") == UnitValue("1 µm")
    assert UnitValue("1 µm") == UnitValue("1e-6 m")

def test_op_neq_uv() : 
    assert not(UnitValue("1 µm") != 1)
    assert not(1 != UnitValue("1 µm"))
    assert UnitValue("1 µm") != 2
    assert 2 != UnitValue("1 µm")
    assert UnitValue("1 µm") != UnitValue("2 µm")
    assert UnitValue("1 µm") != UnitValue("1 m")
    assert not(UnitValue("1 µm") == UnitValue("1e-7 m"))

def test_u_mu_substitution() : 
    assert UnitValue("1 µm") == UnitValue("1 um")
    assert UnitValue("1 µm3") == UnitValue("1 um3")
    assert UnitValue("1 µm-2/s") == UnitValue("1 um-2/s")
    assert UnitValue("1 molecule/µm1") == UnitValue("1 molecule/um")
    
# run all ###############################################

def run_all_tests() :
    test_parse_units()
    test_parse_unitvalue()
    test_unitvalue_conversion()
    test_op_addition_uv_number()
    test_op_substraction_uv_number()
    test_op_multiplication_uv_number()
    test_op_addition_number_uv()
    test_op_substraction_number_uv()
    test_op_multiplication_number_uv()
    test_op_division_number_uv()
    test_op_addition_uv_uv()
    test_op_division_uv_uv()
    test_op_division_uv_number()
    test_op_multiplication_uv_uv()
    test_op_power_uv_number()
    test_op_modulo_uv_number()
    test_op_modulo_number_uv()
    test_op_modulo_uv_uv()
    test_op_modulo_uv_ua() 
    test_op_neg_ua()
    test_ua_invert()
    test_op_multiplication_uv_ua()
    test_op_division_uv_ua()
    test_op_addition_uv_ua()
    test_op_substraction_uv_ua()
    test_op_multiplication_ua_uv()
    test_op_division_ua_uv()
    test_op_addition_ua_uv()
    test_op_substraction_ua_uv()
    test_op_modulo_ua_number()
    test_op_modulo_ua_uv()
    test_op_modulo_ua_ua()
    test_op_addition_ua_number()
    test_op_addition_number_ua()
    test_op_substraction_ua_number()
    test_op_substraction_number_ua()
    test_op_multiplication_number_ua()
    test_op_multiplication_ua_number()
    test_op_division_number_ua()
    test_op_division_ua_number()
    test_op_abs_uv()
    test_op_abs_ua()
    test_op_equals_uv()
    test_op_neq_uv()
    test_u_mu_substitution()