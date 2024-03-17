working with units
==================

*UnitsValue*, *UnitArray* and *Units*
-------------------------------------

1) Presentation

Strengths represents physical quantities with the *UnitValue* and *UnitArray* classes, associating respecitvely a numerical
value or an array of numerical values to some physical units. Units themselves are represented by the *Units* class.
Here is an example of how a *UnitValue* obsects can be declared:

.. code:: python

  a = UnitValue(5, "µM")
  a = UnitValue(5, "uM")
  a = UnitValue("5 µM")
  a = parse_unitvalue(5, "µM")
  a = UnitValue(5, Units("µM"))

All the declaration above are equivalent, setting the variable a as 5 micromolar.
Now, examples for the *UnitArray* class:

.. code:: python

  a = UnitArray([1,2,3], "µM")
  a = UnitArray([1,2,3], Units("µM"))

Both of the above declarations creates an array of values in micromolar (1 µM, 2 µM and 3 µM).

*UnitValue* and *Units* objects can be created from a string, as shown above. the UnitValue string representation
is made of a number litteral followed by a Units string representation, separated by one or more whitespaces

.. code:: python

  UnitValue("1 µm/s")        # OK
  UnitValue("1.5 µm/s")      # OK
  UnitValue("+1.3e-10 µm/s") # OK
  UnitValue("-1.3e-10 µm/s") # OK

  UnitValue("1µm/s")         # wrong. missing the white space between the value and the units.
  UnitValue("a µm/s")        # wrong. a is not a number
  UnitValue("[1, 2] µm/s")   # wrong.
  UnitValue("{'v', 1} µm/s") # wrong.

The syntax for *Units* strings is quite simple:
Those should be a list of symbols separated by "." (multiplication) or "/" (division).

.. code:: python

  parse_units("mol/µm.s")  # OK
  parse_units("mol/µm. s") # wrong, white space in the unit expression.
  parse_units("mol//µm.s") # wrong, successive "/" have no meaning.

A symbol can be immediately followed by a positive non signed or negative exponent.
If absent, the exponent is assumed to be 1. The exponent must be an integer.

.. code:: python

  parse_units("mol.µm-1.s-2")   # OK
  parse_units("mol/µm/s2")      # OK, same units as the previous one
  parse_units("mol1/µm1/s2")    # OK, same units as the previous one, but the
                                # positive exponent 1 is not required
  parse_units("mol.µm-1.5.s-2") # wrong, exponent must be integers
  parse_units("mol.µm+1.s-2")   # wrong, positive exponents must not be signed
  parse_units("mol.µm 1.s-2")   # wrong, white space in the unit expression.

Supported unit symbols are:

+-----+----+--------+------+-------+
|space|time|quantity|volume|density|
+=====+====+========+======+=======+
|     |h   |        |      |       |
+-----+----+--------+------+-------+
|km   |    |kmol    |kL    |kM     |
+-----+----+--------+------+-------+
|     |min |        |      |       |
+-----+----+--------+------+-------+
|m    |s   |mol     |L     |M      |
+-----+----+--------+------+-------+
|dm   |ds  |dmol    |      |dM     |
+-----+----+--------+------+-------+
|cm   |cs  |cmol    |      |cM     |
+-----+----+--------+------+-------+
|mm   |ms  |mmol    |mL    |mM     |
+-----+----+--------+------+-------+
|dmm  |    |        |      |       |
+-----+----+--------+------+-------+
|cmm  |    |        |      |       |
+-----+----+--------+------+-------+
|µm   |µs  |µmol    |µL    |µM     |
+-----+----+--------+------+-------+
|nm   |ns  |nmol    |nL    |nM     |
+-----+----+--------+------+-------+
|pm   |ps  |pmol    |pL    |pM     |
+-----+----+--------+------+-------+
|fm   |fs  |fmol    |fL    |fM     |
+-----+----+--------+------+-------+
|     |    |molecule|      |       |
+-----+----+--------+------+-------+

As a consequence, units are case sensitive, since m corresopnds to meters
while M corresponds to molars (mol/L).

Also, since the µ letter may be inconvenient to type
with a lot of keyboards, it can be substituted by the letter u.
Thus, um, us, umol, uL and uM will be treated as µm, µs, µmol, µL and µM.
ie.

.. code:: python

  UnitValue("1 um") == UnitValue("1 µm") # True

2) Units conversion

The units class describes units as a product of fundamental units (space, time and quantity), the units system component,
raised to some exponent, the units dimensions component/

.. math::

  \textrm{units} = \prod_i {u_i}^{e_i}

where u is the vector of fundamental units, or units system,
and v is the vector of the corresponding exponents, ot units dimensions.

*UnitValue* and *UnitArray* objects thus support unit conversion

.. code:: python

  C = UnitValue(5, "µM")
  C = C.convert("M")
  print(C) #print "5e-6 M"

which is carried out by computing and appyling a conversion factor *f*, defined as

.. math::
  f = \prod_i {\left( \frac{u_i}{v_i} \right) }^{e_i}

where *u* and *v* are units system, the conversion being from *u* to *v*.

For *UnitArrray* objects, conversion is done element wise:

.. code:: python

  C = UnitValue([1, 2, 3], "µm")
  C = C.convert("m")
  print(C) #print "[1e-6, 2e-6, 3e-6] m"

Operations on *UnitValue* and *UnitArrays*
------------------------------------------

operations involving a *UnitValue* and a *UnitArray* are performed
for every element of the *UnitArray*, returning a *UnitArray* with the same size.

operations involving two *UnitArray* are performed are performed between elemnt with the same index
and expects both terms to have the same length.

addition and substraction
^^^^^^^^^^^^^^^^^^^^^^^^^

Both terms must have the same units dimensions.
Number terms are supposed to have the same units as the *UnitValue*/*UnitArray* term:

* UnitValue + UnitValue
* UnitValue + number
* number + UnitValue
* UnitValue + UnitArray

* UnitArray + UnitValue
* UnitArray + number
* number + UnitArray
* UnitArray + UnitArray

division and multiplication
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Terms can have the any units dimensions.
Number terms are supposed to be unitless:

* UnitValue + UnitValue
* UnitValue + number
* number + UnitValue
* UnitValue + UnitArray

* UnitArray + UnitValue
* UnitArray + number
* number + UnitArray
* UnitArray + UnitArray

modulo
^^^^^^

Terms can have the any units dimensions.
Number terms are supposed to have the same units as the UnitValue/UnitArray term:

* UnitValue % UnitValue
* UnitValue % number
* number % UnitValue
* UnitValue % UnitArray

* UnitArray % UnitValue
* UnitArray % number
* number % UnitArray
* UnitArray % UnitArray

power
^^^^^

* UnitValue ** number
* UnitArray ** number

comparison
^^^^^^^^^^

* UnitValue == UnitValue
* UnitValue == number
* UnitValue != UnitValue
* UnitValue != number

For the rest, compared UnitValue object must have the same units dimensions:

* UnitValue > UnitValue
* UnitValue > number
* UnitValue >= UnitValue
* UnitValue >= number
* UnitValue < UnitValue
* UnitValue < number
* UnitValue <= UnitValue
* UnitValue <= number
