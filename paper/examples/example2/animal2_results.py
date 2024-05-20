import strengths as strn
import strengths.plot as strnplt
import matplotlib.pyplot as plt
import numpy as np

xmaxa = []

for l in ["detail", "0", "1", "2"] : 
    out = strn.load_rdtrajectory("animal2_output_(" + l + ").json")    
    xmaxa.append(max(out.data.value.copy().reshape((out.nsamples(), 2, out.system.space.size()))[:, out.system.network.get_species_index("A"), :].reshape(out.nsamples()*out.system.space.size())))

xmax = np.max(xmaxa)

out = strn.load_rdtrajectory("animal2_output_(detail).json")   

plt.gcf().dpi = 300
strnplt.plot_environments_2D(out.system, env_color_dict={"ext" : (1,1,1), "a" : (0,1,0), "b" : (1, 0, 0), "c" : (0, 0, 1)})
plt.show()

for i in range(out.nsamples()) :    
    plt.gcf().dpi = 150
    strnplt.plot_sample_state_2D(out, "A", i, environments=["a", "b", "c"], xmin=0, xmax=xmax, units_system=strn.UnitsSystem(time="h", quantity="molecule"))    
    
for run in range(3) : 
    out = strn.load_rdtrajectory("animal2_output_(" + str(run) + ").json")   
    
    for i in range(out.nsamples()) :    
        plt.gcf().dpi = 300
        strnplt.plot_sample_state_2D(out, "A", i, environments=["a", "b", "c"], xmin=0, xmax=xmax, units_system=strn.UnitsSystem(time="h", quantity="molecule"))     
    