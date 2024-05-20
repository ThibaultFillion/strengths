import matplotlib.pyplot as plt
import matplotlib.cm as pltcm
import matplotlib.colors as pltcol
import numpy as np

import strengths as strn

out = strn.load_rdtrajectory("sim_output_1.json")  
w = out.system.space.w
env = np.array(out.system.space.cell_env)

A = out.data.value.reshape((out.nsamples(), 2, w*w*w))[:,0,:].reshape((out.nsamples()*w*w*w))
state_min = min(A)
state_max = max(A)
cmap = pltcm.ScalarMappable(norm=pltcol.Normalize(vmin=state_min, vmax=state_max))

for t in ["0 h", "300 h", "650 h", "3450 h"] :
    sample = out.get_sample_index(t)
    state = out.get_state("A", sample)
    assert str(state.units) == "molecule"
    state = state.value
    
    col = np.zeros(w*w*w*3)
    
    for i in range(w*w*w) :
        col_i = cmap.to_rgba(state[i])
        col[3*i]   = col_i[0] 
        col[3*i+1] = col_i[1]
        col[3*i+2] = col_i[2]
                
    fig, ax = plt.subplot_mosaic("ab", dpi=300, figsize=(5.5,4), width_ratios=[1,0.08], per_subplot_kw={"a" : {"projection" : "3d"}, "b" : {}})
    fig.suptitle("t = "+str(out.t.get_at(sample).convert("h").value)+" h")
    ax["a"].voxels(filled=env.reshape((w,w,w)).transpose((2,1,0)), facecolors=col.reshape((w,w,w,3)).transpose((2,1,0,3)))
    ax["a"].elev = 15
    ax["a"].set_xlabel("x (cell)")
    ax["a"].set_ylabel("y (cell)")
    ax["a"].set_zlabel("z (cell)")
    plt.colorbar(cmap, cax=ax["b"], label="quantity (molecule)")