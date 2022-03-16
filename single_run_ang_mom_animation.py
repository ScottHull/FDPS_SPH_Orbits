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

path = base_path = "/home/theia/scotthull/Paper1_SPH/gi/5_b073_new_high/5_b073_new_high"
to_path = "{}_spec_angular_momentum".format(path.split("/")[-1])
if not os.path.exists(to_path):
    os.mkdir(to_path)
start_iteration = 0
end_iteration = 500
increment = 2
number_processes = 500
square_scale = 6e7
min_normalize = 0
max_normalize = 2e11
normalizer = Normalize(min_normalize, max_normalize)
cmap = cm.get_cmap('jet')

for iteration in np.arange(start_iteration, end_iteration + increment, increment):
    to_fname = "merged_{}_{}.dat".format(iteration, randint(0, 100000))
    cf = CombineFile(num_processes=number_processes, time=iteration, output_path=path, to_fname=to_fname)
    combined_file = cf.combine()
    formatted_time = round(cf.sim_time * 0.000277778, 2)
    f = os.getcwd() + "/{}".format(to_fname)
    df = pd.read_csv(f, skiprows=2, header=None, delimiter="\t")
    os.remove(f)
    x, y, z = df[3], df[4], df[5]
    vx, vy, vz = df[6], df[7], df[8]
    positions = list(zip(x, y, z))
    velocities = list(zip(vx, vy, vz))
    spec_am = [cmap(normalizer(np.linalg.norm(np.cross(i, j)))) for i, j in zip(positions, velocities)]

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)
    ax.set_xlim(-square_scale, square_scale)
    ax.set_ylim(-square_scale, square_scale)
    ax.set_xticks([], minor=False)
    ax.set_yticks([], minor=False)
    ax.axes.set_aspect('equal')
    ax.scatter(
        x, y, s=0.1,
        color=spec_am, marker="."
    )
    ax.set_title("5b073n-high ({} - {} hrs)".format(iteration, formatted_time))
    sm = cm.ScalarMappable(norm=normalizer, cmap=cmap)
    sm.set_array([])
    cbaxes = inset_axes(ax, width="30%", height="3%", loc=2, borderpad=1.8)
    cbar = plt.colorbar(sm, cax=cbaxes, orientation='horizontal')
    cbar.ax.tick_params(labelsize=6)
    # cbar.ax.set_title("Entropy", fontsize=6)
    cbar.ax.set_title("Specific Angular Momentum ($m^2$/s)", fontsize=6)
    cbar.ax.yaxis.get_offset_text().set(size=6)  # change exponent font size
    cbar.ax.xaxis.get_offset_text().set(size=6)  # change exponent font size
    plt.savefig(to_path + "/{}.png".format(iteration))
