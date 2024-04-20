import strengths.text_array_rw as strntxt
import strengths.coarsegrain as strncg
import strengths as strn
import strengths.plot as strnplt

import matplotlib.pyplot as plt
import numpy as np

tlp = strn.load_rdtrajectory("traj_eul2.json")
eul = strn.load_rdtrajectory("traj_tlp2.json")

tlpcg = strn.load_rdtrajectory("traj_cg_eul2.json")
eulcg = strn.load_rdtrajectory("traj_cg_tlp2.json")

if not( tlp.t.units == strn.Units("s") and\
        eul.t.units == strn.Units("s") and\
        tlpcg.t.units == strn.Units("s") and\
        eulcg.t.units == strn.Units("s") ) :
    raise ValueError()

if not( tlp.data.units == strn.Units("molecule") and\
        eul.data.units == strn.Units("molecule") and\
        tlpcg.data.units == strn.Units("molecule") and\
        eulcg.data.units == strn.Units("molecule") ) :
    raise ValueError()

plt.gcf().dpi = 600
plt.title("global trajectory of Y")
plt.xlabel("time (s)")
plt.ylabel("quantity (molecules)")
plt.plot(tlp.t.value, tlp.get_trajectory("Y", merge=True).value, label="tau leap")
plt.plot(eul.t.value, eul.get_trajectory("Y", merge=True).value, label="Euler")

plt.plot(tlpcg.t.value, tlpcg.get_trajectory("Y", merge=True).value, label="tau leap (coarse grained)")
plt.plot(eulcg.t.value, eulcg.get_trajectory("Y", merge=True).value, label="Euler (coarse grained)")

plt.legend(loc="best")
plt.show()

for t in ["0 s", "100 s", "1500 s"] :
    xmaxa = [max(tlp.get_state("Y", i).value) for i in range(tlp.nsamples())]
    xmax = max(xmaxa)

    sample = tlp.get_sample_index(t)
    plt.gcf().dpi = 300
    strnplt.plot_sample_state_2D(tlp,   "Y", sample, environments=["mmb", "cyt"], xmin=0, xmax=xmax)

cgmap = strntxt.load_1D_array_txt("cg_map.txt", int)

tlpcgucg = strncg.uncoarsegrain_trajectory(tlpcg, tlp.system, cgmap)

for t in ["0 s", "100 s", "1500 s"] :
    xmaxa = [max(tlpcgucg.get_state("Y", i).value) for i in range(tlpcgucg.nsamples())]
    xmax = max(xmaxa)

    sample = tlpcgucg.get_sample_index(t)
    plt.gcf().dpi = 300
    strnplt.plot_sample_state_2D(tlpcgucg,   "Y", sample, environments=["mmb", "cyt"], xmin=0, xmax=xmax)
