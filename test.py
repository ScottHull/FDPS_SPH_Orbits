import os
import shutil
import numpy as np
from math import cos
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable

from src.identify import ParticleMap
from src.combine import CombineFile
from src.animate import animate
from src import vapor

time = 3000
number_processes = 100
path = "/Users/scotthull/Desktop/merged_3000.dat"

pm = ParticleMap(path=path, center=True, relative_velocity=False)
particles = pm.collect_particles()
pm.solve(particles=particles)

vmf = vapor.calc_vapor_mass_fraction(particles=particles)
vapor.plot_disk_entropy(particles=particles)

disk = [p for p in particles if p.label == "DISK"]

low_density_frac = sum([p.mass for p in disk if p.density < 10]) / sum([p.mass for p in disk]) * 100.0
print("low density mass frac: {}".format(low_density_frac))



fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
cmap = plt.get_cmap("viridis")
scales = [p.density for p in particles if p.label == "DISK" and p.density < 1000]
norm = plt.Normalize(min(scales), max(scales))
ax.scatter(
    [p.position[0] for p in disk],
    [p.position[1] for p in disk],
    marker="+",
    c=scales,
    label="DISK"
)
sm = ScalarMappable(norm=norm, cmap=cmap)
cbar = fig.colorbar(sm, ax=ax)
cbar.ax.set_title("Density")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title("Time: {} sec (iteration: {})".format(None, time))
ax.grid()
ax.legend(loc="upper left")
# ax.set_xlim(-1e8, 1e8)
# ax.set_ylim(-1e8, 1e8)
plt.savefig("{}.png".format(time), format='png')
