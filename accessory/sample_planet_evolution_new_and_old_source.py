#!/usr/bin/env python3
import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.font_manager as fm
from matplotlib.colors import ListedColormap

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
label_text = "Entropy"
inc = (max_time - min_time) / sample_interval
sample_time = 3000
square_scale = 1e7

all_iterations_and_times = get_all_iterations_and_times(number_processes=number_processes, path=new_path,
                                                        min_iteration=min_iteration, max_iteration=max_iteration)

plt.style.use("dark_background")
ncol = 2
fig, axs = plt.subplots(1, ncol, figsize=(10, 10),
                        gridspec_kw={"hspace": 0.0, "wspace": 0.0})
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
    ax.scatter(
        [p.position[0] for p in particles if p.position[2] < 0 and p.tag == 0],
        [p.position[1] for p in particles if p.position[2] < 0 and p.tag == 0],
        s=0.02,
        marker="o",
        label="Target Silicate"
    )
    ax.scatter(
        [p.position[0] for p in particles if p.position[2] < 0 and p.tag == 1],
        [p.position[1] for p in particles if p.position[2] < 0 and p.tag == 1],
        s=0.02,
        marker="o",
        label="Target Iron"
    )
    ax.scatter(
        [p.position[0] for p in particles if p.position[2] < 0 and p.tag == 2],
        [p.position[1] for p in particles if p.position[2] < 0 and p.tag == 2],
        s=0.02,
        marker="o",
        label="Impactor Silicate"
    )
    ax.scatter(
        [p.position[0] for p in particles if p.position[2] < 0 and p.tag == 3],
        [p.position[1] for p in particles if p.position[2] < 0 and p.tag == 3],
        s=0.02,
        marker="o",
        label="Impactor Iron"
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
                               square_scale / 5,
                               '{:.2e} km'.format(square_scale / 5),
                               loc=8,
                               pad=0.3,
                               color='white',
                               frameon=False,
                               size_vertical=1,
                               fontproperties=fm.FontProperties(size=6))

    ax.add_artist(scalebar)
    return ax

new_particles, new_time = get_particles(path=new_path, number_processes=number_processes, time=sample_time)
old_particles, old_time = get_particles(path=old_path, number_processes=number_processes, time=sample_time)
ax1 = plot(fig=fig, axs=axs, index=0, time=new_time, particles=new_particles)
ax2 = plot(fig=fig, axs=axs, index=1, time=new_time, particles=old_particles)

axs.flatten()[0].set_title("New EoS")
axs.flatten()[1].set_title("Old EoS")
legend = ax1.legend(fontsize=6)
for handle in legend.legendHandles:
    handle.set_sizes([3.0])
plt.savefig("planet_evolution_source.png", format='png', dpi=200)
