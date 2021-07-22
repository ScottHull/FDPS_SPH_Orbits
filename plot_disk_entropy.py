import os
import shutil
import numpy as np
import matplotlib.pyplot as plt

from src.identify import ParticleMap
from src.combine import CombineFile
from src.animate import animate

time = 3000
number_processes = 100
path = "/scratch/shull4/gi"

cf = CombineFile(num_processes=number_processes, time=time, output_path=path)
formatted_time = cf.sim_time
combined_file = cf.combine()
f = os.getcwd() + "/merged_{}.dat".format(time)
pm = ParticleMap(path=f, center=True, relative_velocity=False)
particles = pm.collect_particles(find_orbital_elements=True)
os.remove(f)

disk_particles = [p for p in particles if p.label == "DISK"]
avg_disk_entropy = sum([p.entropy for p in disk_particles]) / len(disk_particles)

fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
ax.scatter(
    [p.distance / 1000 for p in disk_particles],
    [p.entropy for p in disk_particles],
    marker="+",
    color='black',
    label="Disk Particles"
)
ax.axhline(avg_disk_entropy, color='red', linewidth=2.0, linestyle="--", label="Avg. Disk Entropy")
ax.set_xlabel("Distance from Center (km)")
ax.set_ylabel("Disk Entropy")
ax.grid()
ax.legend()
plt.savefig("disk_entropy.png", format='png')

