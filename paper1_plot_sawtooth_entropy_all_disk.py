#!/usr/bin/env python3
import os
import shutil
import csv
import pandas as pd
import numpy as np
import string
from random import randint
import multiprocessing as mp
import matplotlib.pyplot as plt

from src.animate import animate
from src.combine import CombineFile
from src.identify import ParticleMap
from src.report import get_sim_report, write_report_at_time, build_latex_table_from_disk_report, rows_map


plt.rcParams.update({'font.size': 16, })
# plt.style.use("dark_background")
plt.style.use('seaborn-colorblind')

base_path = "/home/theia/scotthull/Paper1_SPH/gi/"
angle = "b073"
cutoff_densities = [5]
min_density = 0
max_density = 20
min_iteration = 0
max_iteration = 1800
mid_iteration = 800
endstate_iteration = max_iteration
increment = 1
number_processes = 200

headers = ["id", "tag", "mass", "x", "y", "z", "vx", "vy", "vz", "density", "internal energy", "pressure",
           "potential energy", "entropy", "temperature"]

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
    return names, titles


def get_endstate(s, iteration):
    path = base_path + "{}/circularized_{}".format(s, s)
    df = pd.read_csv(path + "/{}.csv".format(iteration))
    return df[df['label'] == "DISK"]


# make a 3 column plot of density, entropy, internal energy
ylabels = [r'Density (kg/m$^3$)', r'Entropy (J/kg/K)', r'Internal Energy (kJ)']

names, titles = get_all_sims()
for s, t in zip(names, titles):
    if os.path.exists(f"paper1_sawtooth_{s}"):
        shutil.rmtree(f"paper1_sawtooth_{s}")
    os.mkdir(f"paper1_sawtooth_{s}")
    # get the endstate df
    endstate = get_endstate(s, endstate_iteration)
    midstate = get_endstate(s, mid_iteration)
    # get a list of the difference in entropy between the endstate and midstate
    select_particles = [i for i in endstate['id'].tolist() if endstate[endstate['id'] == i]['entropy'].values[0] - midstate[midstate['id'] == i]['entropy'].values[0] > 500]
    # get the color cycle
    time = {i: [] for i in endstate['id'].values}
    entropy = {i: [] for i in endstate['id'].values}
    internal_energy = {i: [] for i in endstate['id'].values}
    density = {i: [] for i in endstate['id'].values}
    temperature = {i: [] for i in endstate['id'].values}
    sel_time = {i: [] for i in select_particles}
    sel_entropy = {i: [] for i in select_particles}
    sel_internal_energy = {i: [] for i in select_particles}
    sel_density = {i: [] for i in select_particles}
    sel_temperature = {i: [] for i in select_particles}
    for iteration in np.arange(min_iteration, max_iteration + increment, increment):
        print(f"Working on {s} at iteration {iteration}...")
        path = base_path + "{}/{}".format(s, s)
        to_fname = "merged_{}_{}.dat".format(iteration, randint(0, 100000))
        cf = CombineFile(num_processes=number_processes, time=iteration, output_path=path, to_fname=to_fname)
        df = cf.combine_df()
        formatted_time = round(cf.sim_time * 0.000277778, 2)
        for i in endstate['id'].values:
            time[i].append(formatted_time)
            entropy[i].append(df[df['id'] == i]['entropy'].values[0])
            internal_energy[i].append(df[df['id'] == i]['internal energy'].values[0])
            density[i].append(df[df['id'] == i]['density'].values[0])
            temperature[i].append(df[df['id'] == i]['temperature'].values[0])
        for i in select_particles:
            sel_time[i].append(formatted_time)
            sel_entropy[i].append(df[df['id'] == i]['entropy'].values[0])
            sel_internal_energy[i].append(df[df['id'] == i]['internal energy'].values[0])
            sel_density[i].append(df[df['id'] == i]['density'].values[0])
            sel_temperature[i].append(df[df['id'] == i]['temperature'].values[0])

    fig, axs = plt.subplots(1, 4, figsize=(24, 6), sharex='all')
    axs = axs.flatten()
    ylabels = [r'Density (kg/m$^3$)', r'Entropy (J/kg/K)', r'Internal Energy (kJ)', r'Temperature (K)']
    for ax, y in zip(axs, ylabels):
        ax.grid()
        ax.set_xlabel("Time (hours)", fontsize=16)
        ax.set_ylabel(y, fontsize=16)

    for i in endstate['id'].values:
        axs[0].plot(time[i], density[i], alpha=0.1)
        axs[1].plot(time[i], entropy[i], alpha=0.1)
        axs[2].plot(time[i], internal_energy[i], alpha=0.1)
        axs[3].plot(time[i], temperature[i], alpha=0.1)
    for i in select_particles:
        axs[0].plot(sel_time[i], sel_density[i], alpha=1)
        axs[1].plot(sel_time[i], sel_entropy[i], alpha=1)
        axs[2].plot(sel_time[i], sel_internal_energy[i], alpha=1)
        axs[3].plot(sel_time[i], sel_temperature[i], alpha=1)

    plt.tight_layout()
    plt.savefig(f"{cutoff_densities[0]}_{angle}_sawtooth_entropy_all_disk.png", dpi=300)
    break



