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
output = "/scratch/shull4/thermal_map"

if os.path.exists(output):
    shutil.rmtree(output)
os.mkdir(output)

cf_end = CombineFile(num_processes=number_processes, time=time, output_path=path)
formatted_time_end = cf_end.sim_time
combined_file_end = cf_end.combine()
f = os.getcwd() + "/merged_{}.dat".format(time)
pm_end = ParticleMap(path=f, center=False, relative_velocity=False)
particles = pm_end.collect_particles()
# pm_end.solve(particles=particles)
os.remove(f)

iron_particles = [p for p in particles if p.tag % 2 != 0]
silicate_particles = [p for p in particles if p.tag % 2 == 0]

fig, ax = plt.subplots(figsize=(16, 9))
ax.scatter(
    [p.position[0] for p in silicate_particles],
    [p.position[1] for p in silicate_particles],
    c=[p.entropy for p in silicate_particles],
    cmap='viridis',
    marker='+', 
)
ax.scatter(
    [p.position[0] for p in iron_particles],
    [p.position[1] for p in iron_particles],
    c=[p.entropy for p in iron_particles],
    cmap='cividis',
    marker='+', 
)
