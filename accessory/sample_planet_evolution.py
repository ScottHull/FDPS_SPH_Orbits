#!/usr/bin/env python3
import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm

from src.identify import ParticleMap
from src.combine import CombineFile
from src.time import get_nearest_iteration_to_time, seconds_to_hours, get_all_iterations_and_times

min_time = 0.0
max_time = 3.000030e+05
min_iteration = 0
max_iteration = 3000
number_processes = 100
sample_interval = 3
path = "/home/theia/scotthull/sph_simulations/gi_new_eos"
inc = (max_time - min_time) / sample_interval
sample_times = [0, 15, 20, 30, 500, 2000, 3000]
square_scale = 2e7

all_iterations_and_times = get_all_iterations_and_times(number_processes=number_processes, path=path,
                                                        min_iteration=min_iteration, max_iteration=max_iteration)

# fig, axs = plt.subplots(sample_interval + 1, 1, figsize=(8, 16), sharex='all',
#                         gridspec_kw={"hspace": 0.0})
fig, axs = plt.subplots(len(sample_times), 1, figsize=(8, 16), sharex='all',
                        gridspec_kw={"hspace": 0.0})
fig.patch.set_facecolor('xkcd:black')


# for index, time in enumerate(np.arange(min_time, max_time + inc, inc)):
for index, time in enumerate(sample_times):
    # closest_iteration_to_time = get_nearest_iteration_to_time(time=time, sampled_times=all_iterations_and_times)
    # cf = CombineFile(num_processes=number_processes, time=closest_iteration_to_time, output_path=path)
    cf = CombineFile(num_processes=number_processes, time=time, output_path=path)
    combined_file = cf.combine()
    formatted_time = cf.sim_time
    # f = os.getcwd() + "/merged_{}.dat".format(closest_iteration_to_time)
    f = os.getcwd() + "/merged_{}.dat".format(time)
    pm = ParticleMap(path=f, center=True, relative_velocity=False)
    particles = pm.collect_particles(find_orbital_elements=False)
    # pm.solve(particles=particles)
    os.remove(f)

    ax = axs.flatten()[index]
    ax.set_facecolor('xkcd:black')
    ax.spines['left'].set_color('white')
    ax.spines['right'].set_color('white')
    ax.spines['bottom'].set_color('white')
    ax.spines['top'].set_color('white')
    ax.scatter(
        [p.position[0] for p in particles if p.tag % 2 == 0 and p.position[2] < 0],
        [p.position[1] for p in particles if p.tag % 2 == 0 and p.position[2] < 0],
        s=0.02,
        marker="o",
        color='red',
        label='silicate'
    )
    ax.scatter(
        [p.position[0] for p in particles if p.tag % 2 != 0 and p.position[2] < 0],
        [p.position[1] for p in particles if p.tag % 2 != 0 and p.position[2] < 0],
        s=0.02,
        marker="o",
        color='blue',
        label='iron'
    )
    ax.text(
        square_scale - (square_scale / 2.2),
        square_scale - (square_scale / 1.2),
        str(round(seconds_to_hours(formatted_time), 2)) + " hrs",
        c="white",
        fontsize=10
    )
    # ax.set_title(seconds_to_hours(time))
    ax.set_xticks([])
    # for minor ticks
    ax.set_xticks([], minor=True)
    ax.set_yticks([])
    # for minor ticks
    ax.set_yticks([], minor=True)
    ax.set_xlim(-square_scale, square_scale)
    ax.set_ylim(-square_scale, square_scale)
    ax.set_box_aspect(1)

    scalebar = AnchoredSizeBar(ax.transData,
                               square_scale / 5,
                               '{:.2e} km'.format(square_scale / 5),
                               loc=8,
                               pad=0.2,
                               color='white',
                               frameon=False,
                               size_vertical=1,
                               fontproperties=fm.FontProperties(size=6))

    ax.add_artist(scalebar)

plt.savefig("planet_evolution.png", format='png')
