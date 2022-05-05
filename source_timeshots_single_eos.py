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
iterations = [500, 1000, 1400, 1800]
cutoff_densities = [5, 500, 1000, 2000]
high = True
square_scale = 6e7
min_normalize = 0
max_normalize = 2e11

new_phase_path = "src/phase_data/forstSTS__vapour_curve.txt"
old_phase_path = "src/phase_data/duniteN__vapour_curve.txt"

phase_path = new_phase_path
if runs == "old":
    phase_path = old_phase_path



def get_all_sims(angle, high=True):
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


sims, titles = get_all_sims(angle, high)
fig, axs = plt.subplots(len(iterations), len(sims), figsize=(20, 25), sharex='all', sharey='all',
                        gridspec_kw={"hspace": 0.0, "wspace": 0.0})
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
        planet, disk, escape = df[df['label'] == "PLANET"], df[df['label'] == "DISK"], df[df['label'] == "ESCAPE"]
        for i, label in zip([planet, disk, escape], ["Planet", "Disk", "Escape"]):
            axs[current_index].scatter(
                i['x'], i['y'], s=0.1, marker=".", alpha=1, label=label
            )
        axs[current_index].text(square_scale - (0.5 * square_scale), -square_scale + (0.3 * square_scale),
                                "{}\n{} hrs".format(t, formatted_time), fontsize=10)
        current_index += 1

plt.savefig("source_scenes_{}_{}.png".format(angle, runs), format='png', dpi=200)
