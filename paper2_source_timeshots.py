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

from src.vapor import calc_vapor_mass_fraction_from_formatted
from src.geometry import get_impact_geometry_from_formatted, get_velocity_profile_from_formatted
from src.animate import animate
from src.identify import ParticleMap
from src.combine import CombineFile

plt.rcParams.update({'font.size': 14, })
plt.style.use("dark_background")
# plt.style.use('seaborn-pastel')

base_path = "/home/theia/scotthull/"
runs = [
    [base_path + f"Paper1_SPH/gi/500_b073_new", "Canonical", 200],
    [base_path + f"Paper2_SPH/gi/500_half_earths", "Half-Earths", 200],
    # [base_path + f"Paper2_SPH/gi/mars", "Mars", 200],
]
iterations = [50, 100, 200, 300, 1800]

square_scale = 6e7 / 10 ** 7
min_normalize = 0
max_normalize = 2e11
end_iteration = 1800

phase_path = "src/phase_data/forstSTS__vapour_curve.txt"

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
        run_name = s.split("/")[-1]
        end_state_file = s + "/circularized_{}/{}.csv".format(run_name, end_iteration)
        end_state_df = pd.read_csv(end_state_file, index_col="id")
        endstates.update({t: end_state_df})
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
        run_name = s.split("/")[-1]
        p = s + "/circularized_{}".format(run_name)
        f = p + "/{}.csv".format(iteration)
        df = pd.read_csv(f)
        path = s + "/{}".format(run_name)
        file_format = "results.{}_{}_{}.dat"
        p2 = s + "/{}".format(run_name)
        base_file = file_format.format(
            str(iteration).zfill(5),
            str(number_processes).zfill(5),
            str(0).zfill(5)
        )
        formatted_time = get_time(p2 + "/" + base_file)
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
    # increase axis font size
    ax.tick_params(axis='both', which='major', labelsize=18)

axs[0].annotate("x ($10^4$ km)", xy=(0.0, -5.5), ha="center", fontsize=16, weight='bold')
axs[0].annotate("y ($10^4$ km)", xy=(-5.5, 0.0), va="center", rotation=90, fontsize=16, weight='bold')

plt.tight_layout()
plt.savefig("paper2_source_scenes.png", format='png', dpi=300)
