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
plt.style.use('seaborn-colorblind')

base_path = "/home/theia/scotthull/Paper1_SPH/gi/"
runs = "new"
angle = "b073"
other = "pressure"
units = r"Pa"
# other = "temperature"
# units = r"K"
iterations = [100, 200, 500, 1800]
cutoff_densities = [5, 500, 1000, 2000]
high = True
square_scale = 6e7 / 10 ** 7
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
# fig, axs = plt.subplots(len(iterations), len(sims), figsize=figsize, sharex='all', sharey='all')
# fig, axs = plt.subplots(len(sims), len(iterations), figsize=figsize, sharex='all', sharey='all', gridspec_kw={"hspace": 0.0, "wspace": 0.0})
fig, axs = plt.subplots(len(iterations), len(sims), figsize=figsize, sharex='all', sharey='all')

axs = axs.flatten()
for ax in axs:
    ax.axes.set_aspect('equal')
current_index = 0
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
for iteration in iterations:
    for s, t in zip(sims, titles):
        cd = cutoff_densities.index(int(s.split("_")[0]))
        # color = colors[cd]
        # if "high" in s or "low" in s:
        #     color = colors[-1]
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
        disk_at_rho_cutoff = disk[disk['density'] == cutoff_densities[cd]]
        for i, label in zip([planet, disk, escape, disk_at_rho_cutoff], ["Planet", "Disk", "Escape", r"Disk At $\rho_c$"]):
            axs[current_index].scatter(
                i['density'], i['pressure'], s=0.8, marker=".", alpha=1, color='black', label=label
            )
        if current_index % len(sims) == 0:
            # label time in upper right corner
            axs[current_index].annotate(
                f"{formatted_time} hrs.", xy=(0.95, 0.95), xycoords="axes fraction",
                horizontalalignment="left", verticalalignment="top", fontweight="bold", fontsize=20
            )
        current_index += 1

legend = axs[0].legend(loc='upper right', fontsize=20)
for handle in legend.legendHandles:
    try:
        handle.set_sizes([200.0])
    except:
        pass
# plt.tight_layout()
# plt.margins(0.005, tight=True)

for index, t in enumerate(titles):
    axs[index].set_title(t, fontsize=20)
fig.tight_layout()
fig.subplots_adjust(wspace=0, hspace=0)
for ax in axs[-len(sims):-2]:
    nbins_x = len(ax.get_xticklabels())
    ax.xaxis.set_major_locator(MaxNLocator(nbins=nbins_x, prune='upper'))
for ax in [axs[i] for i in np.arange(len(sims) * 2, len(iterations) * len(sims), len(sims))]:
    nbins_y = len(ax.get_yticklabels())
    ax.yaxis.set_major_locator(MaxNLocator(nbins=nbins_y, prune='upper'))

letters = list(string.ascii_lowercase)
for index, ax in enumerate(axs):
    x1, x2, y1, y2 = ax.axis()
    x_loc = x1 + (0.02 * (x2 - x1))
    y_loc = y2 - (0.08 * (y2 - y1))
    ax.text(x_loc, y_loc, letters[index], fontweight="bold", fontsize=20)

axs[0].annotate(r"Density (kg/m$^3$)", xy=(0.0, -5.5), ha="center", fontsize=16, weight='bold')
axs[0].annotate(f"{other} ({units})", xy=(-5.5, 0.0), va="center", rotation=90, fontsize=16, weight='bold')

plt.savefig("density_vs_{}_{}_{}.png".format(other, angle, runs), format='png', dpi=300)
