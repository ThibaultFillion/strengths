import matplotlib.pyplot as plt
import matplotlib.cm as pltcm
import matplotlib.colors as pltcol
import numpy as np

import strengths as strn
from strengths.text_array_rw import save_1D_array_txt, load_1D_array_txt

w=63
env = load_1D_array_txt("env.txt", int)
cgmap = load_1D_array_txt("cgmap.txt", int)

system = strn.rdsystem_from_dict({
    "network" : {
        "units" : {"time" : "h", "space" : "µm", "quantity" : "molecule"},
        "environments" : ["ext", "a", "b", "c"],
        "species" : [
            {"label" : "A", "density" : {"ext" : 0, "default" : 0.1}, "D" : {"ext" : 0, "default" : 80}},
            {"label" : "B", "density" : {"ext" : 0, "default" : 0.1}, "D" : {"ext" : 0, "default" : 80}}
            ],
        "reactions" : [
            {"eq" : " -> A", "k+" : 0.000105, "k-" : 0.001, "env" : ["b"]},
            {"eq" : " -> B", "k+" : 0.000105, "k-" : 0.001, "env" : ["c"]},
            
            {"eq" : " -> A", "k+" : 0.0001, "k-" : 0.001, "env" : ["a", "c"]},
            {"eq" : " -> B", "k+" : 0.0001, "k-" : 0.001, "env" : ["a", "b"]},
            
            {"eq" : "2 A + B -> 3 A", "k+" : 1},
            {"eq" : "2 B + A -> 3 B", "k+" : 1}
            ]
        },
    "space" : {
        "w" : w,
        "h" : w,
        "d" : w,
        "cell_vol" : 343000,
        "cell_env" : env
        }
    })

out = strn.simulate(system, strn.UnitArray(list(range(0, 3500, 50)), "h"), time_step="1 h", engine=strn.tauleap_engine(), print_progress=True, cgmap=cgmap)
strn.save_rdtrajectory(out, "sim_output_1.json")  

    