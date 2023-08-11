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

base_path = "/home/theia/scotthull/Paper2_SPH/gi/"
num_processes = 600
iterations = [50, 100, 200, 500, 1000]
paths = [['500_mars', "Mars " + r"($b=0.73$)"]]
square_scale = 1e7 / 10 ** 7
min_normalize = 0
max_normalize = 1500
end_iteration = 10000

normalizer = LogNorm(min_normalize, max_normalize)
cmap = cm.get_cmap('cool')

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


def get_endstate(s, endstate_iteration):
    path = base_path + "{}/circularized_{}".format(s, s)
    return pd.read_csv(path + "/{}.csv".format(endstate_iteration))


fig, axs = plt.subplots(len(iterations), len(paths), figsize=(20, 20), sharex='all', sharey='all')

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
    for s, t in paths:
        cd = int(s.split("_")[0])
        path = base_path + "{}/{}".format(s, s)
        to_fname = "merged_{}_{}.dat".format(iteration, randint(0, 100000))
        cf = CombineFile(num_processes=num_processes, time=iteration, output_path=path, to_fname=to_fname)
        combined_file = cf.combine()
        formatted_time = round(cf.sim_time * 0.000277778, 2)
        f = os.getcwd() + "/{}".format(to_fname)
        headers = ["id", "tag", "mass", "x", "y", "z", "vx", "vy", "vz", "density", "internal energy", "pressure",
                   "potential energy", "entropy", "temperature"]
        df = pd.read_csv(f, skiprows=2, header=None, delimiter="\t", names=headers, index_col='id')
        os.remove(f)
        axs[current_index].scatter(
            df['x'] / 10 ** 7, df['y'] / 10 ** 7, s=0.8, marker=".", alpha=1, color=cmap(normalizer(df['entropy']))
        )
        if current_index % len(paths) == 0:
            axs[current_index].text(square_scale - (0.7 * square_scale), -square_scale + (0.44 * square_scale),
                                    "{} hrs".format(formatted_time), fontsize=20)
        current_index += 1


sm = cm.ScalarMappable(norm=normalizer, cmap=cmap)
sm.set_array([])
cbaxes = inset_axes(axs[0], width="50%", height="5%", loc=1, borderpad=1.8)
cbar = plt.colorbar(sm, cax=cbaxes, orientation='horizontal')
cbar.ax.tick_params(labelsize=12)
cbar.ax.set_title("Entropy" + " (J/kg/K)", fontsize=14)


for index, (s, t) in enumerate(paths):
    axs[index].set_title(t, fontsize=20)
fig.tight_layout()
fig.subplots_adjust(wspace=0, hspace=0)
for ax in axs[-len(paths):-2]:
    nbins_x = len(ax.get_xticklabels())
    ax.xaxis.set_major_locator(MaxNLocator(nbins=nbins_x, prune='upper'))
for ax in [axs[i] for i in np.arange(len(paths) * 2, len(iterations) * len(paths), len(paths))]:
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

plt.savefig("paper3_source_scenes.png", format='png', dpi=300)
