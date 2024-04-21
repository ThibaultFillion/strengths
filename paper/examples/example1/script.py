import strengths as strn
import strengths.plot as strnplt

##################################################

system = strn.load_rdsystem("system.json")
import numpy as np

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

strn.save_rdtrajectory(output_euler,   "traj_tlp2.json")
strn.save_rdtrajectory(output_tauleap, "traj_eul2.json")
