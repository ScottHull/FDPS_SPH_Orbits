import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from src.identify import ParticleMap
from src.combine import CombineFile
from src.animate import animate
from src.cross_section import cross_section_xy, sort_particles_by_closest

# time = 3000
# number_processes = 100
# path = "/scratch/shull4/gi"
# output = "/scratch/shull4/thermal_map"
#
# if os.path.exists(output):
#     shutil.rmtree(output)
# os.mkdir(output)
#
# cf_end = CombineFile(num_processes=number_processes, time=time, output_path=path)
# formatted_time_end = cf_end.sim_time
# combined_file_end = cf_end.combine()
# f = os.getcwd() + "/merged_{}.dat".format(time)
f = "/Users/scotthull/Downloads/merged_3000.dat"
pm_end = ParticleMap(path=f, center=True, relative_velocity=False)
particles = pm_end.collect_particles()
# pm_end.solve(particles=particles)
# os.remove(f)

# particles = cross_section_xy(particles=particles, min_z=-100 * 10 ** 3, max_z=100 * 10 ** 3)
particles = sort_particles_by_closest(particles=particles)

print("{} particles sampled".format(len(particles)))
iron_particles = [p for p in particles if p.tag % 2 != 0 and p.label == "PLANET"]
silicate_particles = [p for p in particles if p.tag % 2 == 0 and p.label == "PLANET"]

sns.set()

fig, ax = plt.subplots(figsize=(16, 9))
a0 = ax.scatter(
    [p.position[0] for p in silicate_particles],
    [p.position[1] for p in silicate_particles],
    c=[p.entropy for p in silicate_particles],
    cmap='Reds',
    marker='+',
    label="Silicate"
)
a1 = ax.scatter(
    [p.position[0] for p in iron_particles],
    [p.position[1] for p in iron_particles],
    c=[p.entropy for p in iron_particles],
    cmap='Blues',
    marker='*',
    label="Iron"
)
cbar1 = fig.colorbar(a0, ax=ax)
cbar1.ax.set_ylabel('Entropy (Silicate)', rotation=270)
cbar2 = fig.colorbar(a1, ax=ax)
cbar2.ax.set_ylabel('Entropy (Iron)', rotation=270)
ax.legend(loc='upper left')
ax.set_xlim(-1e7, 1e7)
ax.set_ylim(-1e7, 1e7)

plt.show()
