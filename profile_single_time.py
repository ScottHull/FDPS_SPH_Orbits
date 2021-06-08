import os
import shutil
import numpy as np
import matplotlib.pyplot as plt

from src.identify import ParticleMap
from src.combine import CombineFile
from src.animate import animate

time = 3000
number_processes = 100
path = "/scratch/shull4/gi_new"

cf = CombineFile(num_processes=number_processes, time=time, output_path=path)
formatted_time = cf.sim_time
combined_file = cf.combine()
f = os.getcwd() + "/merged_{}.dat".format(time)
pm = ParticleMap(path=f, center=False, relative_velocity=False)
particles = pm.collect_particles()
pm.solve(particles=particles)
if f in os.listdir(os.getcwd()):
    os.remove(f)

planet = [p for p in particles if p.label == "PLANET"]
disk = [p for p in particles if p.label == "DISK"]
escape = [p for p in particles if p.label == "ESCAPE"]

fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
ax.scatter(
    [p.position[0] for p in planet],
    [p.position[1] for p in planet],
    marker="+",
    color="blue",
    label="PLANET"
)
ax.scatter(
    [p.position[0] for p in disk],
    [p.position[1] for p in disk],
    marker="+",
    color="green",
    label="DISK"
)
ax.scatter(
    [p.position[0] for p in escape],
    [p.position[1] for p in escape],
    marker="+",
    color="red",
    label="ESCAPE"
)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title("Time: {} sec (iteration: {})".format(formatted_time, time))
ax.grid()
ax.legend(loc="upper left")
ax.set_xlim(-1e8, 1e8)
ax.set_ylim(-1e8, 1e8)

plt.savefig("{}.png".format(time), format='png')
