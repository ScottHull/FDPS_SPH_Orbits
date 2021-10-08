import os
import shutil
import numpy as np
import matplotlib.pyplot as plt

from src.identify import ParticleMap
from src.combine import CombineFile
from src.animate import animate
from src.time import get_nearest_iteration_to_time, seconds_to_hours

min_time = 0.0
max_time = 3.000030e+05
min_iteration = 0
max_iteration = 3000
number_processes = 100
sample_interval = 6
path = "sph_simulations/gi_new_eos"
inc = (max_time - min_time) / sample_interval
sample_times = []

fig, axs = plt.subplots(sample_interval, 1, figsize=(8, 8), constrained_layout=True)

for index, time in enumerate(np.arange(min_time, max_time + inc, inc)):
    closest_iteration_to_time = get_nearest_iteration_to_time(time=time, number_processes=number_processes, path=path,
                                                              min_iteration=min_iteration, max_iteration=max_iteration)
    cf = CombineFile(num_processes=number_processes, time=closest_iteration_to_time, output_path=path)
    formatted_time = cf.sim_time
    combined_file = cf.combine()
    f = os.getcwd() + "/merged_{}.dat".format(closest_iteration_to_time)
    pm = ParticleMap(path=f, center=True, relative_velocity=False)
    particles = pm.collect_particles(find_orbital_elements=False)
    # pm.solve(particles=particles)
    os.remove(f)

    ax = axs.flatten()[index]
    ax.scatter(
        [p.position[0] for p in particles if p.tag % 2 == 0],
        [p.position[1] for p in particles if p.tag % 2 == 0],
        s=4,
        marker="o",
        color='red',
        label='silicate'
    )
    ax.scatter(
        [p.position[0] for p in particles if p.tag % 2 != 0],
        [p.position[1] for p in particles if p.tag % 2 != 0],
        s=4,
        marker="o",
        color='blue',
        label='iron'
    )
    ax.set_title(seconds_to_hours(time))
    ax.set_xticks([])
    # for minor ticks
    ax.set_xticks([], minor=True)

plt.savefig("sample_planet_evolution.png", format='png')

