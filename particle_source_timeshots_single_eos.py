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

new_phase_path = "src/phase_data/forstSTS__vapour_curve.txt"
old_phase_path = "src/phase_data/duniteN__vapour_curve.txt"

phase_path = new_phase_path
if runs == "old":
    phase_path = old_phase_path

run_set = ["{}_{}_{}".format(i, angle, runs) for i in [5, 500, 1000, 2000]]
paths = [base_path + "{}/formatted_{}".format(i, i) for i in run_set]

fig, axs = plt.subplots(len(iterations), len(paths), figsize=(20, 25), sharex='all', sharey='all',
                        gridspec_kw={"hspace": 0.0, "wspace": 0.0})
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
        target_silicate = df[df['tag'] == 0]
        target_iron = df[df['tag'] == 1]
        impactor_silicate = df[df['tag'] == 2]
        impactor_iron = df[df['tag'] == 3]
        labels = ["Target Silicate", "Target Iron", "Impactor Silicate", "Impactor Iron"]
        order = [impactor_silicate, impactor_iron, target_silicate, target_iron]
        for index, i in enumerate(order):
            axs[current_index].scatter(
                df['x'], df['y'], s=0.1, marker=".", label=labels[index]
            )
        axs[current_index].text(square_scale - (0.5 * square_scale), -square_scale + (0.3 * square_scale),
                                "{}\n{} hrs".format(get_name(p.split("/")[-1].split("_")[1]), formatted_time), fontsize=10)
        current_index += 1

legend = axs[0].legend(fontsize=6)
for handle in legend.legendHandles:
    handle.set_sizes([3.0])
plt.savefig("tag_source_scenes_{}_{}.png".format(angle, runs), format='png', dpi=200)
