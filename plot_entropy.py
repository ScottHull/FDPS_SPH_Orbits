import os
import shutil
import numpy as np
import matplotlib.pyplot as plt

from src.identify import ParticleMap
from src.combine import CombineFile
from src.animate import animate

time = 3000
number_processes = 100
path = "/scratch/shull4/target"

cf = CombineFile(num_processes=number_processes, time=time, output_path=path)
formatted_time = cf.sim_time
combined_file = cf.combine()
f = os.getcwd() + "/merged_{}.dat".format(time)
pm = ParticleMap(path=f, center=True, relative_velocity=False)
particles = pm.collect_particles(find_orbital_elements=False)
os.remove(f)

fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
ax.scatter(
    [p.distance / 1000 for p in particles if p.tag % 2 == 0],
    [p.entropy for p in particles if p.tag % 2 == 0],
    marker="+",
    color='red',
    label="Silicate"
)
ax.scatter(
    [p.distance / 1000 for p in particles if p.tag % 2 != 0],
    [p.entropy for p in particles if p.tag % 2 != 0],
    marker="+",
    color='blue',
    label="Iron"
)
ax.set_xlabel("Distance from Center (km)")
ax.set_ylabel("Entropy")
ax.grid()
ax.legend()
plt.savefig("entropy.png", format='png')

