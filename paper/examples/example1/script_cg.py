import strengths as strn
import strengths.plot as strnplt
import strengths.text_array_rw as strntxt
import strengths.coarsegrain as strncg
import numpy as np

##################################################

system = strn.load_rdsystem("system.json")
cgmap = strntxt.load_1D_array_txt("cg_map.txt", int)
graph = strncg.coarsegrain_grid(system.space, cgmap)
system = strn.RDSystem(system.network, graph)

##################################################

output_euler = strn.simulate(
    system,
    t_sample = strn.UnitArray(np.linspace(0, 3600, 360), "s"),
    time_step = "0.5 ms",
    engine = strn.euler_engine(),
    print_progress=True
    )

output_tauleap = strn.simulate(
    system,
    t_sample = strn.UnitArray(np.linspace(0, 3600, 360), "s"),
    time_step = "0.5 ms",
    engine = strn.tauleap_engine(),
    print_progress=True
    )

strn.save_rdtrajectory(output_euler,   "traj_cg_tlp2.json")
strn.save_rdtrajectory(output_tauleap, "traj_cg_eul2.json")
