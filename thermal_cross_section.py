import os
import shutil
import numpy as np
import matplotlib.pyplot as plt

from src.identify import ParticleMap
from src.combine import CombineFile
from src.animate import animate

start_time = 0
end_time = 3000
interval = 20
number_processes = 100
path = "/scratch/shull4/gi"
output = "/scratch/shull4/thermal_map"

if os.path.exists(output):
    shutil.rmtree(output)
os.mkdir(output)

cf_end = CombineFile(num_processes=number_processes, time=end_time, output_path=path)
formatted_time_end = cf_end.sim_time
combined_file_end = cf_end.combine()
f = os.getcwd() + "/merged_{}.dat".format(end_time)
pm_end = ParticleMap(path=f, center=False, relative_velocity=False)
particles = pm_end.collect_particles()
pm_end.solve(particles=particles)
os.remove(f)

iron_particles = [i for i in particles if p.tag % 2 != 0]
silicate_particles = [i for i in particles if p.tag % 2 == 0]

_min, _max = np.amin(np.array([i for i in particles.entropy])), np.amax([i for i in particles.entropy])
norm = plt.Normalize(_min, _max)
fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
ax.scatter(
    [p.position[0] for p in silicate_particles],
    [p.position[1] for p in silicate_particles],
    c=[p.entropy for p in silicate_particles],
    cmap='viridis',
    marker='+', 
    norm=norm
)
plt.clim(_min, _max)
ax.scatter(
    [p.position[0] for p in iron_particles],
    [p.position[1] for p in iron_particles],
    c=[p.entropy for p in iron_particles],
    cmap='cividis',
    marker='+', 
    norm=norm
)
plt.clim(_min, _max)
ax.set_xlabel("x")
ax.set_ylabel("y")
plt.colorbar().set_label('Entropy', rotation=270)

plt.savefig("entropy_cross_section.png", format='png')
