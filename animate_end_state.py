#!/usr/bin/env python3
import os
import csv
import pandas as pd
import numpy as np
from math import atan, pi, acos
from operator import contains
from random import randint
import multiprocessing as mp
import matplotlib.pyplot as plt

from src.combine import CombineFile

plt.rcParams.update({'font.size': 8, })
# plt.style.use("dark_background")
plt.style.use('seaborn-colorblind')

min_iteration = 0
max_iteration = 1800
increment = 20
angle = 'b073'
cutoff_densities = [5, 500, 1000, 2000]
number_processes = 200
square_scale = 6e7
base_path = "/home/theia/scotthull/Paper1_SPH/gi/"


def get_all_sims(angle, high=True):
    fformat = "{}_{}_{}"
    tformat = "{}{}{}"
    names = []
    titles = []
    for runs in ["new", "old"]:
        n = "s"
        if runs == "old":
            n = "n"
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


def get_com(x, y, z, mass):
    total_mass = float(np.sum(mass))
    x_center = sum(x * mass) / total_mass
    y_center = sum(y * mass) / total_mass
    z_center = sum(z * mass) / total_mass
    return np.array([x_center, y_center, z_center])


def get_end_states(angle, high):
    endstates = {}
    sims, titles = get_all_sims(angle, high)
    for s, t in zip(sims, titles):
        end_state_file = base_path + "{}/circularized_{}/{}.csv".format(s, s, max_iteration)
        end_state_df = pd.read_csv(end_state_file, index_col="id")
        endstates.update({t: end_state_df})
    return endstates


def plot_iteration(iteration, time, dfs, end_dfs, to_path):
    num_new = len([i for i in dfs.keys() if "n" in i])
    num_old = len([i for i in dfs.keys() if "o" in i])
    num_rows = max([num_new, num_old])
    plt.style.use("dark_background")
    fig, axs = plt.subplots(num_rows, 2, figsize=(16, 32), sharex='all',
                            gridspec_kw={"hspace": 0.10, "wspace": 0.12})
    fig.patch.set_facecolor('xkcd:black')
    axs = axs.flatten()
    for ax in axs:
        ax.set_xlim(-square_scale, square_scale)
        ax.set_ylim(-square_scale, square_scale)
        ax.grid(alpha=0.1)
    index_new, index_old = 0, 1
    for t, df in dfs.items():
        to_index = index_new
        if "o" in t:
            to_index = index_old
        end_df = end_dfs[t]
        planet = end_df[end_df['label'] == "PLANET"]
        disk = end_df[end_df['label'] == "DISK"]
        escape = end_df[end_df['label'] == "ESCAPE"]
        to_planet = df[df['id'].isin(planet.index.tolist())]
        to_disk = df[df['id'].isin(disk.index.tolist())]
        to_escape = df[df['id'].isin(escape.index.tolist())]
        for d, label in zip([to_planet, to_disk, to_escape], ["PLANET", "DISK", 'ESCAPE']):
            axs[to_index].scatter(
                d['x'], d['y'], marker=".", s=1, label=label
            )
            axs[to_index].set_title("{} {} hrs. ({})".format(t, time, iteration))
        if "o" in t:
            index_old += 2
        else:
            index_new += 2

    legend = axs[0].legend(loc='upper left', fontsize=8)
    for handle in legend.legendHandles:
        try:
            handle.set_sizes([30.0])
        except:
            pass

    plt.savefig(to_path + "/{}.png".format(iteration), format='png', dpi=200)


def run_proc(args):
    iteration, to_path = args
    if not os.path.exists(to_path):
        os.mkdir(to_path)
    data = {}
    high = False
    sims, titles = get_all_sims(angle, high=high)
    formatted_time = None
    endstates = get_end_states(angle=angle, high=high)
    for s, t in zip(sims, titles):
        print("{} - {}".format(iteration, t))
        path = base_path + "{}/{}".format(s, s)
        to_fname = "merged_{}_{}.dat".format(iteration, randint(0, 100000))
        cf = CombineFile(num_processes=number_processes, time=iteration, output_path=path, to_fname=to_fname)
        combined_file = cf.combine()
        formatted_time = round(cf.sim_time * 0.000277778, 2)
        f = os.getcwd() + "/{}".format(to_fname)
        headers = ["id", "tag", "mass", "x", "y", "z", "vx", "vy", "vz", "density", "internal energy", "pressure",
                   "potential energy", "entropy", "temperature"]
        df = pd.read_csv(f, skiprows=2, header=None, delimiter="\t", names=headers)
        os.remove(f)
        zero_coords = get_com(df[df['tag'] == 1]['x'], df[df['tag'] == 1]['y'],
                              df[df['tag'] == 1]['z'], df[df['tag'] == 1]['mass'])
        df['x'] -= zero_coords[0]
        df['y'] -= zero_coords[1]
        df['z'] -= zero_coords[2]

        data.update({t: df})

    plot_iteration(iteration=iteration, time=formatted_time, dfs=data, end_dfs=endstates, to_path=to_path)


def run():
    to_path = "{}_endstates".format(angle)
    pool = mp.Pool(10)
    pool.map(run_proc, [[iteration, to_path] for iteration in
                        np.arange(min_iteration, max_iteration + increment, increment)])
    pool.close()
    pool.join()


run()
