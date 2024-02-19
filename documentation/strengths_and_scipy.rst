Using Strengths with SciPy's integration ODE methods
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`SciPy <https://scipy.org/>`_ [1] proposes powerful tools for integrating 
systems of ordinary differential equations (ODE), 
such as `solve_ivp <https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html#scipy.integrate.solve_ivp>`_ [2], why are way more advanced and fast
than the Euler method available in the built-in engine.
We show here how to use Strengths with Scipy.

Here is a first example using Strengths's built-in engine, with the Euler method :

.. code:: python
  from strengths import *
  import numpy as np
  import matplotlib.pyplot as plt

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

  out = simulate(system, t_sample=np.linspace(0, 10, 1000), engine=euler_engine())
  plt.plot(out.t.value, out.get_trajectory("A").value, label="A")
  plt.plot(out.t.value, out.get_trajectory("B").value, label="B")
  plt.xlabel("t (s)")
  plt.ylabel("n (molecules/µm$^3$)")
  plt.show()

Now, let us do the same thing, but using `solve_ivp <https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html#scipy.integrate.solve_ivp>`_ [2] instead of simulate.
The integration function is obtained by calling the make_dxdtf method of the 
RDSystem object. 

.. code:: python
  from strengths import *
  import numpy as np
  import matplotlib.pyplot as plt
  from scipy.integrate import solve_ivp

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
  f = system.make_dxdtf()
  out = solve_ivp(f, t_span=(0, 10), y0=system.state.value, method="LSODA")
  plt.plot(out.t, out.y[0], label="A")
  plt.plot(out.t, out.y[1], label="B")
  plt.xlabel("t (s)")
  plt.ylabel("n (molecules/µm$^3$)")
  plt.show()

References
^^^^^^^^^^
[1] SciPy package website (https://scipy.org/)

[2] SciPy's API Reference documentation for the solve_ivp function (https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html#scipy.integrate.solve_ivp)
