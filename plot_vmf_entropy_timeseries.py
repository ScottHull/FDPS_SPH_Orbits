#!/usr/bin/env python3
import os
import csv
import pandas as pd
import numpy as np
from random import randint
import multiprocessing as mp
import matplotlib.pyplot as plt

plt.style.use("dark_background")

base_path = "/home/theia/scotthull/Paper1_SPH/gi/"
angle = "b073"
cutoff_densities = [5, 500, 1000, 2000]
min_iteration = 0
max_iteration = 1800
increment = 50
increment_high = 100
number_processes = 200
number_processes_high = 500

new_phase_path = "src/phase_data/forstSTS__vapour_curve.txt"
old_phase_path = "src/phase_data/duniteN__vapour_curve.txt"


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


def get_all_sims(high=True):
    fformat = "{}_{}_{}"
    tformat = "{}{}{}"
    names = []
    titles = []
    high_res_name = None
    high_res_title = None
    for runs in ["new", "old"]:
        n = "n"
        if runs == "old":
            n = "o"
        for cd in cutoff_densities:
            output_name = fformat.format(cd, angle, runs)
            title_name = tformat.format(cd, angle, n)
            titles.append(title_name)
            names.append(output_name)
            if cd == 5 and high and runs == "new":
                high_res_name = fformat.format(cd, angle, runs) + "_high"
                high_res_title = tformat.format(cd, angle, n) + "-high"
    if high_res_name is not None and high_res_title is not None:
        names.append(high_res_name)
        titles.append(high_res_title)
    return names, titles


def plot_entropy_and_vmf_vs_time():
    sims, titles = get_all_sims(high=False)
    fig, axs = plt.subplots(3, 2, figsize=(16, 24), sharex="all",
                            gridspec_kw={"hspace": 0.14, "wspace": 0.14})
    axs = axs.flatten()
    axs[0].set_title("Stewart M-ANEOS")
    axs[1].set_title("GADGET M-ANEOS")
    for ax in axs:
        ax.grid(alpha=0.4)
    new_index, old_index = 0, 1
    d = {}
    for sim, title in zip(sims, titles):
        if title not in d.keys():
            d.update({title: {'TIME_HRS': [], 'MEAN_DISK_ENTROPY': [], 'DISK VMF': [], "DISK_MASS": []}})
        for iteration in np.arange(min_iteration, max_iteration + increment, increment):
            path = base_path + "{}/{}_reports/".format(sim, sim)
            df = pd.read_csv(path + "{}.csv".format(iteration))
            time, avg_disk_entropy, disk_vmf, disk_mass = df['TIME_HRS'][0], df['MEAN_DISK_ENTROPY'][0], \
                                                          df['DISK VMF'][0], df['DISK_MASS'][0]
            d[title]['TIME_HRS'].append(time)
            d[title]['MEAN_DISK_ENTROPY'].append(avg_disk_entropy)
            d[title]['DISK VMF'].append(disk_vmf)
            d[title]['DISK_MASS'].append(disk_mass)
    for i in d.keys():
        for j in d[i].keys():
            d[i][j] = [float(str(k).split(" ")[0]) for k in d[i][j]]
    for sim in d.keys():
        to_index = new_index
        if "o" in sim:
            to_index = old_index
        axs[to_index].plot(
            d[sim]['TIME_HRS'], d[sim]['MEAN_DISK_ENTROPY'], linewidth=2.0, label=sim
        )
        axs[to_index + 2].plot(
            d[sim]['TIME_HRS'], d[sim]['DISK VMF'], linewidth=2.0, label=sim
        )
        axs[to_index + 4].plot(
            d[sim]['TIME_HRS'], d[sim]['DISK_MASS'], linewidth=2.0, label=sim
        )
    for ax in [axs[0], axs[1]]:
        ax.set_ylabel("Avg. Disk Entropy")
    for ax in [axs[2], axs[3]]:
        ax.set_ylabel("Disk VMF (%)")
    for ax in [axs[4], axs[5]]:
        ax.set_ylabel(r"Disk Mass ($M_{L}$)")
        ax.set_xlabel("Time (hrs)")
    axs[0].legend(loc='lower left')
    axs[1].legend(loc='lower left')
    plt.savefig("{}_disk_entropy_and_vmf.png".format(angle), format='png', dpi=200)


plot_entropy_and_vmf_vs_time()
