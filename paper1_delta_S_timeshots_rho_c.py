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
from matplotlib.colors import Normalize, LogNorm
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

base_path = "/home/theia/scotthull/Paper1_SPH/gi/"
runs = "new"
angle = "b073"
iterations = [1400, 1500, 1600, 1700, 1800]
cutoff_densities = [5, 500, 1000, 2000]
high = True
square_scale = 6e7 / 10 ** 7
min_normalize = 0
max_normalize = 1000
end_iteration = 1800

new_phase_path = "src/phase_data/forstSTS__vapour_curve.txt"
old_phase_path = "src/phase_data/duniteN__vapour_curve.txt"

normalizer = Normalize(min_normalize, max_normalize)
cmap = cm.get_cmap('cool')

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
        if cd == 5 and high and runs == "new" and angle == "b073":
            high_res_name = fformat.format(cd, angle, runs) + "_high"
            high_res_title = tformat.format(cd, angle, n) + "-high"
        if cd == 2000 and high and runs == "old" and angle == "b075":
            high_res_name = fformat.format(cd, angle, runs) + "_low"
            high_res_title = tformat.format(cd, angle, n) + "-low"
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
figsize = (20, 20)
if (high and angle == "b073" and runs == "new") or (high and angle == "b075" and runs == "old"):
    figsize = (24.5, 20)
fig, axs = plt.subplots(len(iterations) - 1, len(sims), figsize=figsize, sharex='all', sharey='all')

axs = axs.flatten()
for ax in axs:
    ax.set_xlim(-square_scale, square_scale)
    ax.set_ylim(-square_scale, square_scale)
    # ax.set_xticks([], minor=False)
    # ax.set_yticks([], minor=False)
    ax.axes.set_aspect('equal')
current_index = 0

prev_disk = {s: None for s in sims}
for iteration_index, iteration in enumerate(iterations):
    for s, t in zip(sims, titles):
        cd = cutoff_densities.index(int(s.split("_")[0]))
        p = base_path + "{}/circularized_{}".format(s, s)
        f = p + "/{}.csv".format(iteration)
        df = pd.read_csv(f, index_col="id")
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
        planet, disk, escape = df[df.index.isin(end_planet.index.tolist())], df[
            df.index.isin(end_disk.index.tolist())], df[df.index.isin(end_escape.index.tolist())]
        if iteration_index > 0:
            # get the difference between the entropy of this iteration and the previous iteration for each particle
            # which is indexed by the particle id
            disk['delta_S'] = None
            for i in disk.index.values:
                try:
                    disk.loc[i, 'delta_S'] = disk.loc[i]['entropy'] - prev_disk[s].loc[i]['entropy']
                    # delta_S[i] = disk.loc[i]['entropy'] - prev_disk[s].loc[i]['entropy']
                except KeyError:
                    disk.loc[i, 'delta_S'] = 0
            for i, label in zip([disk], ["Disk"]):
                axs[current_index].scatter(
                    i['x'] / 10 ** 7, i['y'] / 10 ** 7, s=0.8, marker=".", alpha=1,
                    c=cmap(normalizer(i['delta_S'])), label=label
                )
            if current_index % len(sims) == 0:
                axs[current_index].text(square_scale - (0.7 * square_scale), -square_scale + (0.40 * square_scale),
                                        "{} hrs".format(formatted_time), fontsize=20)
            # axs[current_index].text(square_scale - (0.7 * square_scale), -square_scale + (0.30 * square_scale),
            #                         "{} %".format(fraction_at_rho_c), fontsize=20)
            current_index += 1

        prev_disk[s] = df


sm = cm.ScalarMappable(norm=normalizer, cmap=cmap)
sm.set_array([])
cbaxes = inset_axes(axs[0], width="50%", height="5%", loc=1, borderpad=1.8)
cbar = plt.colorbar(sm, cax=cbaxes, orientation='horizontal')
cbar.ax.tick_params(labelsize=12)
cbar.ax.set_title(r"$\Delta$S" + "(J/kg/K)", fontsize=14)

for index, t in enumerate(titles):
    axs[index].set_title(t, fontsize=20)
fig.tight_layout()
fig.subplots_adjust(wspace=0, hspace=0)
for ax in axs[-len(sims):-2]:
    nbins_x = len(ax.get_xticklabels())
    ax.xaxis.set_major_locator(MaxNLocator(nbins=nbins_x, prune='upper'))
for ax in [axs[i] for i in np.arange(len(sims) * 2, (len(iterations) - 1) * len(sims), len(sims))]:
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

plt.savefig("{}_source_scenes_{}_{}.png".format('delta_S', angle, runs), format='png', dpi=300)
