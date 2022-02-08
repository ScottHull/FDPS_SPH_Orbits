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
import src.full_suite_profiling as profile

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

def __build_scene(d):
    meta, iteration, to_path, min_normalize_parameter, max_normalize_parameter, square_scale, s, u, p, to_client_path = d
    client = LunaToTheia(s, u, p)
    # normalizer = Normalize(min_normalize_parameter, max_normalize_parameter)
    # cmap = cm.get_cmap('jet')

    num_new = len([i for i in meta.keys() if "new" in i])
    num_old = len([i for i in meta.keys() if "old" in i])
    num_rows = max([num_new, num_old])
    plt.style.use("dark_background")
    fig, axs = plt.subplots(num_rows, 2, figsize=(16, 32), sharex='all',
                            gridspec_kw={"hspace": 0.10, "wspace": 0.12})
    fig.patch.set_facecolor('xkcd:black')
    axs = axs.flatten()
    for ax in axs:
        ax.grid(alpha=0.4)
        # ax.set_xlim(-square_scale, square_scale)
        # ax.set_ylim(-square_scale, square_scale)
        # ax.set_xticks([], minor=False)
        # ax.set_yticks([], minor=False)
        ax.set_ylim(0, 1e7)
    axs[-2].set_xlabel(r"Radius from Target Center ($R_\oplus$)")
    axs[-1].set_xlabel(r"Radius from Target Center ($R_\oplus$)")
    index_new, index_old = 0, 1
    for i in meta.keys():
        n = meta[i]['name']
        p = meta[i]['path']
        if s is not None:
            f = client.get_file(client.theia_client, p, "{}.csv".format(iteration))
            formatted_time = get_time(f, local=False)
            df = client.get_df_from_theia(p, "{}.csv".format(iteration), skiprows=2)
        else:
            formatted_time = get_time(p + "/{}.csv".format(iteration))
            df = pd.read_csv(p + "/{}.csv".format(iteration), skiprows=2)
        df = df[df['z'] < 0]
        df = df[df['distance'] > (2 * 6371 * 1000)]  # only secondary impactor material
        if "new" in i:
            axs[index_new].scatter(
                df['radius'] / (6371 * 1000),
                [df['pressure'] / df['density']],
                s=2,
            )
            axs[index_new].set_title(n + " {} hrs".format(formatted_time))
            axs[index_new].set_ylabel(r"P / $\rho$")
            index_new += 2
        else:
            axs[index_old].scatter(
                df['radius'] / (6371 * 1000),
                [df['pressure'] / df['density']],
                s=2,
            )
            axs[index_old].set_title(n + " {} hrs".format(formatted_time))
            index_old += 2
    # sm = cm.ScalarMappable(norm=normalizer, cmap=cmap)
    # sm.set_array([])
    # cbaxes = inset_axes(axs[0], width="30%", height="3%", loc=2, borderpad=1.8)
    # cbar = plt.colorbar(sm, cax=cbaxes, orientation='horizontal')
    # cbar.ax.tick_params(labelsize=6)
    # # cbar.ax.set_title("Entropy", fontsize=6)
    # cbar.ax.set_title("Angular Momentum (kg $\cdot m^2$/s )", fontsize=6)
    plt.savefig(to_path + "/{}.png".format(iteration), format='png')
    if s is not None:
        client.send_file_to_theia(to_path, to_client_path, "/{}.png".format(iteration))
        os.remove(to_path + "/{}.png".format(iteration))
    client.theia_client.close()


def build_scenes(name, meta, to_path, start_iteration, end_iteration, min_normalization_param, max_normalization_param,
                 increment, s=None, u=None, p=None, to_client_path="", fill=False, proc=10, square_scale=6e7):
    # if os.path.exists(to_path):
    #     shutil.rmtree(to_path)
    if not fill:
        try:
            os.mkdir(to_path)
        except:
            pass
    pool = mp.Pool(proc)
    to_make = np.arange(start_iteration, end_iteration + increment, increment)
    if fill:
        if s is not None:
            client = LunaToTheia(s, u, p)
            ldir = client.listdir(client.theia_client, to_client_path)
            to_make = [i for i in to_make if "{}.png".format(i) not in ldir]
        else:
            to_make = [i for i in to_make if "{}.png".format(i) not in os.listdir(to_path)]
        to_make = [i for i in to_make if "{}.png".format(i) not in os.listdir(to_path)]
    pool.map(__build_scene, [[meta, iteration, to_path, min_normalization_param, max_normalization_param, square_scale,
                              s, u, p, to_client_path] for iteration in to_make])
    pool.close()
    pool.join()
    animate(
        start_time=start_iteration,
        end_time=end_iteration,
        interval=increment,
        path=to_path,
        fps=30,
        filename="{}_animated_from_formatted.mp4".format(name),
    )


base_path = "/home/theia/scotthull/Paper1_SPH/gi/"
base_path_setups = "/home/theia/scotthull/Paper1_SPH/setups/"

gi_b073_runs = {
    "5_b073_new": {
        "name": "5b073n",
        "path": base_path + "5_b073_new/formatted_5_b073_new",
        "setup": profile.get_setup_file_data(base_path_setups +
                                             "{}/setups/{}".format("setups_b073", "5_new_gi_setup_b_073.txt"))
    },
    "5_b073_old": {
        "name": "5b073o",
        "path": base_path + "5_b073_old/formatted_5_b073_old",
        "setup": profile.get_setup_file_data(base_path_setups +
                                             "{}/setups/{}".format("setups_b073", "5_old_gi_setup_b_073.txt"))
    },
    "500_b073_new": {
        "name": "500b073n",
        "path": base_path + "500_b073_new/formatted_500_b073_new",
        "setup": profile.get_setup_file_data(base_path_setups +
                                             "{}/setups/{}".format("setups_b073", "500_new_gi_setup_b_073.txt"))
    },
    "500_b073_old": {
        "name": "500b073o",
        "path": base_path + "500_b073_old/formatted_500_b073_old",
        "setup": profile.get_setup_file_data(base_path_setups +
                                             "{}/setups/{}".format("setups_b073", "500_old_gi_setup_b_073.txt"))
    },
    "1000_b073_new": {
        "name": "1000b073n",
        "path": base_path + "1000_b073_new/formatted_1000_b073_new",
        "setup": profile.get_setup_file_data(base_path_setups +
                                             "{}/setups/{}".format("setups_b073", "1000_new_gi_setup_b_073.txt"))
    },
    "1000_b073_old": {
        "name": "1000b073o",
        "path": base_path + "1000_b073_old/formatted_1000_b073_old",
        "setup": profile.get_setup_file_data(base_path_setups +
                                             "{}/setups/{}".format("setups_b073", "1000_old_gi_setup_b_073.txt"))
    },
    "2000_b073_new": {
        "name": "2000b073n",
        "path": base_path + "2000_b073_new/formatted_2000_b073_new",
        "setup": profile.get_setup_file_data(base_path_setups +
                                             "{}/setups/{}".format("setups_b073", "2000_new_gi_setup_b_073.txt"))
    },
    "2000_b073_old": {
        "name": "2000b073o",
        "path": base_path + "2000_b073_old/formatted_2000_b073_old",
        "setup": profile.get_setup_file_data(base_path_setups +
                                             "{}/setups/{}".format("setups_b073", "2000_old_gi_setup_b_073.txt"))
    },
}

s = "epsl.earth.rochester.edu"
u = "scotthull"
p = "PW"

build_scenes("gi_b073_runs", gi_b073_runs, s=s, u=u, p=p, proc=30,
                     to_client_path="/home/theia/scotthull/FDPS_SPH_Orbits/ang_mom_gi_b073_runs_scenes",
                     start_iteration=0, end_iteration=350, increment=5, to_path="pres_dens_gi_b073_runs_scenes", fill=False,
                     min_normalization_param=1e34, max_normalization_param=2.5e34, square_scale=6e7)
    