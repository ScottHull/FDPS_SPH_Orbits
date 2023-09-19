#!/usr/bin/env python3
import os
import csv
import shutil
from math import pi, asin, isnan, exp
import numpy as np
import pandas as pd
from random import randint
from statistics import mean
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import multiprocessing as mp
from matplotlib.ticker import MaxNLocator
import string

from SPHDiskIdentification.src.combine import CombinedFile
from SPHDiskIdentification.src.identify import ParticleMap

plt.rcParams.update({'font.size': 14, })
plt.style.use("dark_background")
# plt.style.use('seaborn-pastel')

base_path = "/home/theia/scotthull/Paper3_SPH/gi/"
runs = [
    [base_path + f"500_mars_b073_2v_esc/500_mars_b073_2v_esc", r"Mars ($b=0.73$, $1 v_{esc}$)", 600],
    [base_path + f"500_mars_b073_2v_esc/500_mars_b073_2v_esc", r"Mars ($b=0.73$, $2 v_{esc}$)", 600],
    [base_path + f"500_mars_b073_2v_esc/500_mars_b073_2v_esc", r"Mars ($b=0.50$, $1 v_{esc}$)", 600],
]
iterations = [0, 20, 80, 100, 200]
end_iteration = 200
square_scale = 10 ** 7
mass_mars = 6.39e23  # kg
radius_mars = 3389.5 * 1000  # m

phase_path = "src/phase_data/forstSTS__vapour_curve.txt"

file_headers = ["id", "tag", "mass", "x", "y", "z", "vx", "vy", "vz", "density", "internal energy", "pressure",
                "potential energy", "entropy", "temperature"]

def get_time(f, local=True):
    formatted_time = None
    if local:  # if not reading from remote server
        with open(f, 'r') as infile:
            reader = csv.reader(infile, delimiter="\t")
            formatted_time = float(next(reader)[0])
        infile.close()
    else:
        formatted_time = float(next(f))
    return round(formatted_time * 0.000277778, 2)  # seconds -> hours


def get_end_states():
    endstates = {}
    for s, t, i in runs:
        # create the combined file
        c = CombinedFile(
            path=s,
            iteration=end_iteration,
            number_of_processes=i,
            to_fname=f"merged_{end_iteration}_{randint(1, int(1e5))}.dat"
        )
        combined_file = c.combine_to_memory()
        combined_file.columns = file_headers
        time = c.sim_time
        # create the particle map
        particle_map = ParticleMap(particles=combined_file, mass_planet=mass_mars,
                                   equatorial_radius=radius_mars)
        particles = particle_map.loop()
        endstates[t] = particles
    return endstates

endstates = get_end_states()
# fig, axs = plt.subplots(len(iterations), len(sims), figsize=figsize, sharex='all', sharey='all')
# fig, axs = plt.subplots(len(sims), len(iterations), figsize=figsize, sharex='all', sharey='all', gridspec_kw={"hspace": 0.0, "wspace": 0.0})
fig, axs = plt.subplots(len(iterations), len(runs), figsize=(10, 20), sharex='all', sharey='all')

axs = axs.flatten()
for ax in axs:
    ax.set_xlim(-square_scale, square_scale)
    ax.set_ylim(-square_scale, square_scale)
    # ax.set_xticks([], minor=False)
    # ax.set_yticks([], minor=False)
    ax.axes.set_aspect('equal')
current_index = 0
# colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
for iteration in iterations:
    for s, t, number_processes in runs:
        c = CombinedFile(
            path=s,
            iteration=iteration,
            number_of_processes=number_processes,
            to_fname=f"merged_{iteration}_{randint(1, int(1e5))}.dat"
        )
        df = c.combine_to_memory()
        formatted_time = c.sim_time
        df.columns = file_headers
        endstate = endstates[t]
        # df = df[df['z'] <= 0]  # slice simulation
        end_planet, end_disk, end_escape = endstate[endstate['label'] == "PLANET"], endstate[
            endstate['label'] == "DISK"], endstate[endstate['label'] == "ESCAPE"]
        planet, disk, escape = df[df['id'].isin(end_planet.index.tolist())].sort_values("z"), df[
            df['id'].isin(end_disk.index.tolist())].sort_values("z"), df[
                                   df['id'].isin(end_escape.index.tolist())].sort_values("z")
        for i, label in zip([planet, disk, escape], ["Planet", "Disk", "Escape"]):
            axs[current_index].scatter(
                i['x'] / 10 ** 7, i['y'] / 10 ** 7, s=0.8, marker=".", alpha=1, label=label
            )
        if current_index % len(runs) == 0:
            axs[current_index].text(square_scale - (0.75 * square_scale), -square_scale + (0.3 * square_scale),
                                    "{} hrs".format(formatted_time), fontsize=20)
        current_index += 1

legend = axs[0].legend(loc='upper right', fontsize=20)
for handle in legend.legendHandles:
    try:
        handle.set_sizes([200.0])
    except:
        pass
# plt.tight_layout()
# plt.margins(0.005, tight=True)

for index, t in enumerate([i[1] for i in runs]):
    axs[index].set_title(t, fontsize=22)
fig.tight_layout()
fig.subplots_adjust(wspace=0, hspace=0)
for ax in axs[-len(runs):-2]:
    nbins_x = len(ax.get_xticklabels())
    ax.xaxis.set_major_locator(MaxNLocator(nbins=nbins_x, prune='upper'))
for ax in [axs[i] for i in np.arange(len(runs) * 2, len(iterations) * len(runs), len(runs))]:
    nbins_y = len(ax.get_yticklabels())
    ax.yaxis.set_major_locator(MaxNLocator(nbins=nbins_y, prune='upper'))

letters = list(string.ascii_lowercase)
for index, ax in enumerate(axs):
    x1, x2, y1, y2 = ax.axis()
    x_loc = x1 + (0.02 * (x2 - x1))
    y_loc = y2 - (0.08 * (y2 - y1))
    ax.text(x_loc, y_loc, letters[index], fontweight="bold", fontsize=20)

axs[0].annotate("x ($10^4$ km)", xy=(0.0, -5.5), ha="center", fontsize=16, weight='bold')
axs[0].annotate("y ($10^4$ km)", xy=(-5.5, 0.0), va="center", rotation=90, fontsize=16, weight='bold')

plt.tight_layout()
plt.savefig("paper3_source_scenes.png", format='png', dpi=300)
