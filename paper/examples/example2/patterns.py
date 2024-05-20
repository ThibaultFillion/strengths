import strengths as strn
import strengths.plot as strnplt

import matplotlib.pyplot as plt
import numpy as np

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
        "w" : 400,
        "h" : 1,
        "cell_vol" : 343000,
        "cell_env" : [1 for i in range(400)]
        }
    })

out = strn.simulate(system, strn.UnitArray(list(range(0, 800, 1)), "h"), time_step="1 h", engine=strn.tauleap_engine(), print_progress=True)
strn.save_rdtrajectory(out, "pattern1D_output_1.json")

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
        "w" : 100,
        "h" : 100,
        "cell_vol" : 343000,
        "cell_env" : [1 for i in range(100*100)]
        }
    })

out = strn.simulate(system, strn.UnitArray(list(range(0, 3500, 50)), "h"), time_step="1 h", engine=strn.tauleap_engine(), print_progress=True)
strn.save_rdtrajectory(out, "pattern2D_output_1.json")