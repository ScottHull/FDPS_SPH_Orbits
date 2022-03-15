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

from src.vapor import calc_vapor_mass_fraction_from_formatted
from src.geometry import get_impact_geometry_from_formatted, get_velocity_profile_from_formatted
from src.animate import animate
from src.identify import ParticleMap
from src.combine import CombineFile

plt.style.use("dark_background")

base_path = "/home/theia/scotthull/Paper1_SPH/gi/"
runs = "new"
angle = "b073"
iterations = [85, 155, 205, 260, 300]
square_scale = 6e7
min_normalize = 0
max_normalize = 2e11

new_phase_path = "src/phase_data/forstSTS__vapour_curve.txt"
old_phase_path = "src/phase_data/duniteN__vapour_curve.txt"

phase_path = new_phase_path
if runs == "old":
    phase_path = old_phase_path

run_set = ["{}_{}_{}".format(i, angle, runs) for i in [5, 500, 1000, 2000]]
paths = [base_path + "{}/formatted_{}".format(i, i) for i in run_set]

normalizer = Normalize(min_normalize, max_normalize)
cmap = cm.get_cmap('jet')
fig, axs = plt.subplots(len(iterations), len(paths), figsize=(20, 25), sharex='all', sharey='all')
axs = axs.flatten()
for ax in axs:
    ax.set_xlim(-square_scale, square_scale)
    ax.set_ylim(-square_scale, square_scale)
    ax.set_xticks([], minor=False)
    ax.set_yticks([], minor=False)
    ax.axes.set_aspect('equal')


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

def get_name(cd):
    eos = "n"
    if runs == "old":
        eos = "o"
    return "{}{}{}".format(cd, angle, eos)

current_index = 0
for iteration in iterations:
    for p in paths:
        f = p + "/{}.csv".format(iteration)
        formatted_time = get_time(f)
        df = pd.read_csv(f, skiprows=2)
        positions = list(zip(df['x'], df['y'], df['z']))
        velocities = list(zip(df['vx'], df['vy'], df['vz']))
        spec_am = [cmap(normalizer(np.linalg.norm(np.cross(i, j)))) for i, j in zip(positions, velocities)]
        axs[current_index].scatter(
            df['x'], df['y'], s=0.1,
            color=spec_am, marker="."
        )
        axs[current_index].text(square_scale - (0.5 * square_scale), -square_scale + (0.3 * square_scale),
                                "{}\n{} hrs".format(get_name(p.split("/")[-1].split("_")[1]), formatted_time), fontsize=10)
        current_index += 1

sm = cm.ScalarMappable(norm=normalizer, cmap=cmap)
sm.set_array([])
cbaxes = inset_axes(axs[0], width="30%", height="3%", loc=2, borderpad=1.8)
cbar = plt.colorbar(sm, cax=cbaxes, orientation='horizontal')
cbar.ax.tick_params(labelsize=6)
# cbar.ax.set_title("Entropy", fontsize=6)
cbar.ax.set_title("Specific Angular Momentum ($m^2$/s)", fontsize=6)
cbar.ax.yaxis.get_offset_text().set(size=6)  # change exponent font size
cbar.ax.xaxis.get_offset_text().set(size=6)  # change exponent font size
plt.savefig("am_scenes_{}_{}.png".format(angle, runs), format='png', dpi=200)
