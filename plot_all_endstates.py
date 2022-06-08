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
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import multiprocessing as mp

from src.vapor import calc_vapor_mass_fraction_from_formatted
from src.geometry import get_impact_geometry_from_formatted, get_velocity_profile_from_formatted
from src.animate import animate
from src.identify import ParticleMap
from src.combine import CombineFile

plt.rcParams.update({'font.size': 8, })
plt.style.use("dark_background")

base_path = "/home/theia/scotthull/Paper1_SPH/gi/"
angle = "b073"
iteration = 1800
cutoff_densities = [5, 500, 1000, 2000]
square_scale = 6e7

def get_all_sims(angle, high=False):
    fformat = "{}_{}_{}"
    tformat = "{}{}{}"
    names = []
    titles = []
    for runs in ["new", "old"]:
        n = "n"
        if runs == "old":
            n = "o"
        for cd in cutoff_densities:
            output_name = fformat.format(cd, angle, runs)
            title_name = tformat.format(cd, angle, n)
            titles.append(title_name)
            names.append(output_name)
    if high:
        output_name = fformat.format(5, angle, "new") + "_high"
        names.append(output_name)
        title_name = tformat.format(5, angle, "n") + "-high"
        titles.append(title_name)
    return names, titles

fig, axs = plt.subplots(2, len(cutoff_densities), figsize=(16, 9), sharex='all', sharey='all',
                        gridspec_kw={"hspace": 0.0, "wspace": 0.0})
axs = list(axs.flatten())
index = 0
sims, titles = get_all_sims(angle, False)
for s, t in zip(sims, titles):
    path = base_path + "{}/circularized_{}".format(s, s)
    df = pd.read_csv(path + "/{}.csv".format(iteration))
    planet = df[df['label'] == "PLANET"]
    disk = df[df['label'] == "DISK"]
    escape = df[df['label'] == "ESCAPE"]
    df = df.sort_values(by=['id'])
    df = df.reset_index(drop=True)
    axs[index].scatter(planet['x'] / 1000, planet['y'] / 1000, marker='.', s=1, label="Earth")
    axs[index].scatter(disk['x'] / 1000, disk['y'] / 1000, marker='.', s=1, label="Disk")
    axs[index].scatter(escape['x'] / 1000, escape['y'] / 1000, marker='.', s=1, label="Escaping")
    axs[index].set_xlim([-square_scale / 1000, square_scale / 1000])
    axs[index].set_ylim([-square_scale / 1000, square_scale / 1000])
    x1, x2, y1, y2 = axs[index].axis()
    x_loc = x2 - (0.32 * (x2 - x1))
    y_loc = y2 - (0.1 * (y2 - y1))
    axs[index].text(x_loc, y_loc, t, fontweight="bold", fontsize=12)

    index += 1

legend = axs[0].legend(loc='upper left')
for handle in legend.legendHandles:
    try:
        handle.set_sizes([50.0])
    except:
        pass
plt.tight_layout()
for ax in axs[5:]:
    nbins_x = len(ax.get_xticklabels())
    ax.xaxis.set_major_locator(MaxNLocator(nbins=nbins_x, prune='lower'))
for ax in [axs[4]]:
    nbins_y = len(ax.get_yticklabels())
    ax.yaxis.set_major_locator(MaxNLocator(nbins=nbins_y, prune='lower'))

plt.savefig("{}_endstates.png".format(angle), format='png', dpi=300)
