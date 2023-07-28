#!/usr/bin/env python3
import os
import csv
import pandas as pd
import numpy as np
import string
from random import randint
import multiprocessing as mp
import matplotlib.pyplot as plt

from src.combine import CombineFile
from src.identify import ParticleMap
from src.report import get_sim_report, write_report_at_time, build_latex_table_from_disk_report, rows_map


plt.rcParams.update({'font.size': 16, })
# plt.style.use("dark_background")
plt.style.use('seaborn-colorblind')

base_path = "/home/theia/scotthull/Paper1_SPH/gi/"
angle = "b073"
cutoff_densities = [5]
min_iteration = 0
max_iteration = 1800
endstate_iteration = max_iteration
increment = 1
number_processes = 200

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


def get_endstate(s):
    path = base_path + "{}/circularized_{}".format(s, s)
    df = pd.read_csv(path + "/{}.csv".format(endstate_iteration))
    return df[df['label'] == "DISK"]


# make a 3 column plot of density, entropy, internal energy
fig, axs = plt.subplots(1, 3, figsize=(18, 6), sharex=True)
axs = axs.flatten()
for ax in axs:
    ax.grid()
    ax.set_xlabel("Time (hours)", fontsize=16)
ylabels = [r'Density (kg/m$^3$)', r'Entropy (J/kg/K)', r'Internal Energy (kJ)']
for ax, y in zip(axs, ylabels):
    ax.set_ylabel(y, fontsize=16)

names, titles = get_all_sims()
for s, t in zip(names, titles):
    linestyle = "-"
    if "old" in s:
        linestyle = "--"
    # get the endstate df
    endstate = get_endstate(s)
    endstate_target_particles = endstate[endstate['entropy'] > 9000]
    # get 5 random particle ids from the endstate df
    endstate = endstate_target_particles.sample(n=5)['id'].tolist()
    # get the color cycle
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    time = {i: [] for i in endstate}
    entropy = {i: [] for i in endstate}
    internal_energy = {i: [] for i in endstate}
    density = {i: [] for i in endstate}
    for iteration in np.arange(min_iteration, max_iteration + increment, increment):
        path = base_path + "{}/{}".format(s, s)
        to_fname = "merged_{}_{}.dat".format(iteration, randint(0, 100000))
        cf = CombineFile(num_processes=number_processes, time=iteration, output_path=path, to_fname=to_fname)
        combined_file = cf.combine()
        formatted_time = round(cf.sim_time * 0.000277778, 2)
        f = os.getcwd() + "/{}".format(to_fname)
        headers = ["id", "tag", "mass", "x", "y", "z", "vx", "vy", "vz", "density", "internal energy", "pressure",
                   "potential energy", "entropy", "temperature"]
        df = pd.read_csv(f, skiprows=2, header=None, delimiter="\t", names=headers)
        disk = df[df['id'].isin(endstate)]
        os.remove(f)
        for i in endstate:
            time[i].append(formatted_time)
            entropy[i].append(disk[disk['id'] == i]['entropy'].values[0])
            internal_energy[i].append(disk[disk['id'] == i]['internal energy'].values[0])
            density[i].append(disk[disk['id'] == i]['density'].values[0])

    for ax, i in zip(axs, [density, entropy, internal_energy]):
        for index, j in enumerate(endstate):
            ax.plot(time[j], i[j], linestyle=linestyle, color=colors[index])

axs[-1].plot(
    [], [], linestyle="-", label="Stewart M-ANEOS", color="black"
)
axs[-1].plot(
    [], [], linestyle="--", label="N-SPH M-ANEOS", color="black"
)

axs[0].set_ylim(bottom=0, top=20)

axs[-1].legend(fontsize=14)
plt.tight_layout()
plt.savefig("sawtooth_entropy.png", dpi=300)



