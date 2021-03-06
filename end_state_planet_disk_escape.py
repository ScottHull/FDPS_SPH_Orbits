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
plt.rcParams.update({'font.size': 10,})
plt.style.use("dark_background")

angle = "b073"
cutoff_densities = [5, 500, 1000, 2000]
iteration = 1800
base_path = "/home/theia/scotthull/Paper1_SPH/gi/"
square_scale = -6e7


def get_all_sims(high=True):
    fformat = "{}_{}_{}"
    tformat = "{}{}{}"
    names = []
    titles = []
    for runs in ["new", "old"]:
        n = "S"
        if runs == "old":
            n = "N"
        for cd in cutoff_densities:
            output_name = fformat.format(cd, angle, runs)
            title_name = tformat.format(cd, angle, n)
            titles.append(title_name)
            names.append(output_name)
            if cd == 5 and high and runs == "new" and angle == "b073":
                output_name = fformat.format(cd, angle, runs) + "_high"
                names.append(output_name)
                title_name = tformat.format(cd, angle, n) + "-high"
                titles.append(title_name)
    return names, titles


fig, axs = plt.subplots(2, 4, figsize=(16, 9), sharex="all", sharey="all",
                        gridspec_kw={"hspace": 0.14, "wspace": 0.08})

# add a big axes, hide frame
fig.add_subplot(111, frameon=False)
# hide tick and tick label of the big axes
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.grid(False)
plt.xlabel("X")
plt.ylabel("Y")

axs = axs.flatten()
cur_index = 0
sims, titles = get_all_sims(high=False)
for s, t in zip(sims, titles):
    path = base_path + "{}/circularized_{}/{}.csv".format(s, s, iteration)
    df = pd.read_csv(path)
    planet, disk, escape = df[df['label'] == "PLANET"], df[df['label'] == "DISK"], df[df['label'] == "ESCAPE"]
    axs[cur_index].scatter(
        planet['x'], planet['y'], s=0.2, label="Planet"
    )
    axs[cur_index].scatter(
        disk['x'], disk['y'], s=0.2, label="Disk"
    )
    axs[cur_index].scatter(
        escape['x'], escape['y'], s=0.2, label="Escaping"
    )
    # axs[cur_index].text(square_scale - (0.5 * square_scale), -square_scale + (0.3 * square_scale),
    #                    t, fontsize=10)
    axs[cur_index].set_title(t)
    cur_index += 1

for ax in axs:
    ax.set_xlim(-square_scale, square_scale)
    ax.set_ylim(-square_scale, square_scale)
    # ax.set_xticks([], minor=False)
    # ax.set_yticks([], minor=False)
    ax.axes.set_aspect('equal')
    ax.tick_params(direction="in")

legend = axs[0].legend(loc='upper left')
for handle in legend.legendHandles:
    handle.set_sizes([14.0])

plt.savefig("{}_planet_disk_escape.png".format(angle), format='png', dpi=200)

