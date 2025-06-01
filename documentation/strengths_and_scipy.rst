Using Strengths with SciPy's integration ODE methods
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`SciPy <https://scipy.org/>`_ [1] proposes powerful tools for integrating 
systems of ordinary differential equations (ODE), 
such as `solve_ivp <https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html#scipy.integrate.solve_ivp>`_ [2], why are way more advanced and fast than the Euler method implemented by the ``strenghts.euler_engine()``.
To take advantage of this, you can use the ScipyRDEngine class, which is built on SciPy's ODE solvers :

.. code:: python
  
  from strengths import *
  import numpy as np
  from strengths.scipyrdengine import ScipyRDEngine
  
  system = rdsystem_from_dict({
      "network" : {
          "species" : [
              {"label" : "A", "density" : 100},
              {"label" : "B", "density" : 50},
              ],
          "reactions" : [
              {"eq" : "A -> B", "k+" : 1, "k-" : 1}
              ]
          }
      })

  out = simulate(system, t_sample=np.linspace(0, 10, 1000), engine=ScipyRDEngine())

By default, it uses the Scipy's LSODA solver [3], but it is possible to specfiy another one 

.. code:: python
  
  from strengths import *
  import numpy as np
  from strengths.scipyrdengine import ScipyRDEngine
  from scipy.integrate import Radau
  
  system = rdsystem_from_dict({
      "network" : {
          "species" : [
              {"label" : "A", "density" : 100},
              {"label" : "B", "density" : 50},
              ],
          "reactions" : [
              {"eq" : "A -> B", "k+" : 1, "k-" : 1}
              ]
          }
      })

  out = simulate(system, t_sample=np.linspace(0, 10, 1000), engine=ScipyRDEngine(Radau))

References
^^^^^^^^^^
[1] SciPy package website (https://scipy.org/)

[2] SciPy's API Reference documentation for the solve_ivp function (https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html#scipy.integrate.solve_ivp)

[3] SciPy's API Reference documentation for the LSODA class (https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.LSODA.html#scipy.integrate.LSODA)
