import strengths as strn
import strengths.plot as strnplt
import matplotlib.pyplot as plt
import numpy as np

out = strn.load_rdtrajectory("pattern1D_output_1.json")    
plt.gcf().dpi = 300
plt.title("Quantity of $A$ over space and time")
plt.xlabel("position (cell)")
plt.ylabel("time (h)")
plt.imshow(out.data.value.reshape(out.nsamples(), 2, 400)[:,0,:])
plt.colorbar(label="quantity (molecule)")
plt.show()

out = strn.load_rdtrajectory("pattern2D_output_1.json")    
plt.gcf().dpi = 300
strnplt.plot_sample_state_2D(out, "A", out.nsamples()-1, units_system=strn.UnitsSystem(time="h"))
plt.show()
