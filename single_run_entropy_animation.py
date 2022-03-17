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
from src.theia import LunaToTheia

plt.style.use("dark_background")

path = "/home/theia/scotthull/Paper1_SPH/gi/5_b073_new_high/5_b073_new_high"
remote_path = "/home/theia/scotthull/FDPS_SPH_Orbits/5_b073_new_high_entropy_remote"
to_path = "{}_entropy".format(path.split("/")[-1])
if not os.path.exists(to_path):
    os.mkdir(to_path)
start_iteration = 0
end_iteration = 500
increment = 2
number_processes = 500
square_scale = 6e7
min_normalize = 0
max_normalize = 8000
normalizer = Normalize(min_normalize, max_normalize)
cmap = cm.get_cmap('jet')

s = "epsl.earth.rochester.edu"
u = "scotthull"
p = ""
theia = LunaToTheia(s, u, p)

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


def plot_iteration(args):
    iteration = args
    f = theia.get_and_combine_files_from_iteration(remote_path=path, num_processes=number_processes,
                                                   iteration=iteration, to_base_dir="/scratch/shull4")
    formatted_time = get_time(f)
    df = pd.read_csv(f, skiprows=2, header=None, delimiter="\t")
    os.remove(f)
    x, y, z = df[3], df[4], df[5]
    entropy = df[13]

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)
    ax.set_xlim(-square_scale, square_scale)
    ax.set_ylim(-square_scale, square_scale)
    ax.set_xticks([], minor=False)
    ax.set_yticks([], minor=False)
    ax.axes.set_aspect('equal')
    ax.scatter(
        x, y, s=0.1,
        color=[cmap(normalizer(i)) for i in entropy], marker="."
    )
    ax.set_title("5b073n-high ({} - {} hrs)".format(iteration, formatted_time))
    sm = cm.ScalarMappable(norm=normalizer, cmap=cmap)
    sm.set_array([])
    cbaxes = inset_axes(ax, width="30%", height="3%", loc=2, borderpad=1.8)
    cbar = plt.colorbar(sm, cax=cbaxes, orientation='horizontal')
    cbar.ax.tick_params(labelsize=6)
    # cbar.ax.set_title("Entropy", fontsize=6)
    cbar.ax.set_title("Entropy", fontsize=8)
    cbar.ax.yaxis.get_offset_text().set(size=6)  # change exponent font size
    cbar.ax.xaxis.get_offset_text().set(size=6)  # change exponent font size
    f = to_path + "/{}.png".format(iteration)
    plt.savefig(f, dpi=200)
    theia.send_file_to_theia(to_path, remote_path, "/{}.png".format(iteration))
    os.remove(to_path + "/{}.png".format(iteration))
    os.remove(f)

pool = mp.Pool(5)
pool.map(plot_iteration, [[iteration] for iteration in np.arange(start_iteration, end_iteration + increment, increment)])
pool.close()
pool.join()
