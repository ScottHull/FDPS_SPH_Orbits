#!/usr/bin/env python3
import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
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
new_path = "/home/theia/scotthull/sph_simulations/gi_new_eos"
old_path = "/home/theia/scotthull/sph_simulations/gi_old_eos"
inc = (max_time - min_time) / sample_interval
sample_times = [0, 15, 20, 30, 80, 200, 3000]
square_scale = 2e7

all_iterations_and_times = get_all_iterations_and_times(number_processes=number_processes, path=new_path,
                                                        min_iteration=min_iteration, max_iteration=max_iteration)

nrow = len(sample_times)
ncol = 2
fig, axs = plt.subplots(nrow, ncol, figsize=(8, 16), sharex='all',
                        gridspec_kw={"hspace": 0.0, "wspace": 0.0, "top": 1. - 0.5 / (nrow + 1),
                                     "bottom": 0.5 / (nrow + 1),
                                     "left": 0.5 / (ncol + 1), "right": 1 - 0.5 / (ncol + 1)})
fig.patch.set_facecolor('xkcd:black')


def get_particles(path, number_processes, time):
    cf = CombineFile(num_processes=number_processes, time=time, output_path=path)
    combined_file = cf.combine()
    formatted_time = cf.sim_time
    # f = os.getcwd() + "/merged_{}.dat".format(closest_iteration_to_time)
    f = os.getcwd() + "/merged_{}.dat".format(time)
    pm = ParticleMap(path=f, center=True, relative_velocity=False)
    particles = pm.collect_particles(find_orbital_elements=False)
    # pm.solve(particles=particles)
    os.remove(f)
    return particles, formatted_time


def plot(fig, axs, particles, index, time):
    ax = axs.flatten()[index]
    ax.set_facecolor('xkcd:black')
    ax.spines['left'].set_color('white')
    ax.spines['right'].set_color('white')
    ax.spines['bottom'].set_color('white')
    ax.spines['top'].set_color('white')
    ax.scatter(
        [p.position[0] for p in particles if p.position[2] < 0],
        [p.position[1] for p in particles if p.position[2] < 0],
        s=0.02,
        marker="o",
        c=[p.tag for p in particles if p.position[2] < 0],
    )
    ax.text(
        square_scale - (square_scale / 1.2),
        square_scale - (square_scale / 3),
        str(round(seconds_to_hours(time), 2)) + " hrs",
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
                               6378.1 * 1000,
                               r'1 $R_{\bigoplus}$',
                               loc=8,
                               pad=0.3,
                               color='white',
                               frameon=False,
                               size_vertical=1,
                               fontproperties=fm.FontProperties(size=6))

    ax.add_artist(scalebar)
    return ax


tracked_index = 0
for index, time in enumerate(sample_times):
    new_particles, new_time = get_particles(path=new_path, number_processes=number_processes, time=time)
    old_particles, old_time = get_particles(path=old_path, number_processes=number_processes, time=time)
    plot(fig=fig, axs=axs, index=tracked_index, time=new_time, particles=new_particles)
    tracked_index += 1
    plot(fig=fig, axs=axs, index=tracked_index, time=new_time, particles=old_particles)
    tracked_index += 1

axs.flatten()[0].set_title("New EoS", c="white")
axs.flatten()[1].set_title("Old EoS", c="white")
plt.savefig("planet_evolution.png", format='png')
