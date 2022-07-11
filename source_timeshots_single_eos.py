#!/usr/bin/env python3
import os
import csv
import shutil
from math import pi, asin, isnan
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

from src.vapor import calc_vapor_mass_fraction_from_formatted
from src.geometry import get_impact_geometry_from_formatted, get_velocity_profile_from_formatted
from src.animate import animate
from src.identify import ParticleMap
from src.combine import CombineFile

plt.rcParams.update({'font.size': 10, })
plt.style.use("dark_background")

base_path = "/home/theia/scotthull/Paper1_SPH/gi/"
runs = "new"
angle = "b073"
iterations = [100, 200, 500, 1800]
cutoff_densities = [5, 500, 1000, 2000]
high = True
square_scale = 6e7
min_normalize = 0
max_normalize = 2e11
end_iteration = 1800

new_phase_path = "src/phase_data/forstSTS__vapour_curve.txt"
old_phase_path = "src/phase_data/duniteN__vapour_curve.txt"

phase_path = new_phase_path
if runs == "old":
    phase_path = old_phase_path


def get_all_sims(high=True):
    fformat = "{}_{}_{}"
    tformat = "{}{}{}"
    names = []
    titles = []
    high_res_name = None
    high_res_title = None
    n = "S"
    if runs == "old":
        n = "N"
    for cd in cutoff_densities:
        output_name = fformat.format(cd, angle, runs)
        title_name = tformat.format(cd, angle, n)
        titles.append(title_name)
        names.append(output_name)
        if cd == 5 and high and runs == "new":
            high_res_name = fformat.format(cd, angle, runs) + "_high"
            high_res_title = tformat.format(cd, angle, n) + "-high"
    if high_res_name is not None and high_res_title is not None:
        names.append(high_res_name)
        titles.append(high_res_title)
    return names, titles


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


def get_end_states(angle, high):
    endstates = {}
    sims, titles = get_all_sims(high)
    for s, t in zip(sims, titles):
        end_state_file = base_path + "{}/circularized_{}/{}.csv".format(s, s, end_iteration)
        end_state_df = pd.read_csv(end_state_file, index_col="id")
        endstates.update({t: end_state_df})
    return endstates


sims, titles = get_all_sims(high)
endstates = get_end_states(angle=angle, high=high)
figsize = (24.5, 20)
if high and angle == "b073" and runs == "new":
    figsize = (24.5, 20)
# fig, axs = plt.subplots(len(iterations), len(sims), figsize=figsize, sharex='all', sharey='all')
# fig, axs = plt.subplots(len(sims), len(iterations), figsize=figsize, sharex='all', sharey='all', gridspec_kw={"hspace": 0.0, "wspace": 0.0})
fig, axs = plt.subplots(len(iterations), len(sims), figsize=figsize, sharex='all', sharey='all')

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
    for s, t in zip(sims, titles):
        cd = cutoff_densities.index(int(s.split("_")[0]))
        p = base_path + "{}/circularized_{}".format(s, s)
        f = p + "/{}.csv".format(iteration)
        df = pd.read_csv(f)
        path = base_path + "{}/{}".format(s, s)
        if "high" in s:
            number_processes = 500
        else:
            number_processes = 200
        file_format = "results.{}_{}_{}.dat"
        p2 = base_path + "{}/{}".format(s, s)
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
                i['x'], i['y'], s=0.1, marker=".", alpha=1, label=label
            )
        axs[current_index].text(square_scale - (0.7 * square_scale), -square_scale + (0.3 * square_scale),
                                "{} hrs".format(formatted_time), fontsize=16)
        current_index += 1

# legend = axs[0].legend(loc='upper left', fontsize=14)
# for handle in legend.legendHandles:
#     try:
#         handle.set_sizes([50.0])
#     except:
#         pass
# plt.tight_layout()
# plt.margins(0.005, tight=True)

for index, t in enumerate(titles):
    axs[index].set_title(t, fontsize=18)
fig.tight_layout()
fig.subplots_adjust(wspace=0, hspace=0)
for ax in axs[-5:-1]:
    nbins_x = len(ax.get_xticklabels())
    ax.xaxis.set_major_locator(MaxNLocator(nbins=nbins_x, prune='upper'))
for ax in [axs[5], axs[10], axs[15]]:
    nbins_y = len(ax.get_yticklabels())
    ax.yaxis.set_major_locator(MaxNLocator(nbins=nbins_y, prune='upper'))

plt.savefig("source_scenes_{}_{}.png".format(angle, runs), format='png', dpi=200)
