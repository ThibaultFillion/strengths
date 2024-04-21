import matplotlib.pyplot as plt
import matplotlib.cm as pltcm
import matplotlib.colors as pltcol
import numpy as np
from strengths.text_array_rw import load_1D_array_txt

w=63
env = np.array(load_1D_array_txt("env.txt", int))
envcol = np.zeros(w*w*w*3)

for i in range(w*w*w):
    if env[i] == 2 :
        envcol[3*i] = 1; envcol[3*i+1] = 0; envcol[3*i+2] = 0;
    elif env[i] == 3 :
        envcol[3*i] = 0; envcol[3*i+1] = 0; envcol[3*i+2] = 1;
    else:
        envcol[3*i] = 0; envcol[3*i+1] = 1; envcol[3*i+2] = 0;
                        
fig, ax = plt.subplot_mosaic("ab", dpi=300, figsize=(5.5,4), width_ratios=[1,0.08], per_subplot_kw={"a" : {"projection" : "3d"}, "b" : {}})
fig.suptitle("System layout")
ax["a"].voxels(filled=env.reshape((w,w,w)).transpose((2,1,0)), facecolors=envcol.reshape((w,w,w,3)).transpose((2,1,0,3)))
ax["a"].elev = 15
ax["a"].set_xlabel("x (cell)")
ax["a"].set_ylabel("y (cell)")
ax["a"].set_zlabel("z (cell)")
color_map = pltcol.ListedColormap([(1,1,1), (0,1,0), (1,0,0), (0,0,1)])
plt.colorbar(pltcm.ScalarMappable(norm=pltcol.Normalize(vmin=0, vmax=1), cmap=color_map), cax=ax["b"], label="environments")
ax["b"].set_yticks([1/8, 3/8, 5/8, 7/8], ["$ext$", "$a$", "$b$", "$c$"])
