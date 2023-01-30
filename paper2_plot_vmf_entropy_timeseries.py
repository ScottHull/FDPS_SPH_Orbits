#!/usr/bin/env python3
import os
import csv
import pandas as pd
import numpy as np
import string
from random import randint
import multiprocessing as mp
import matplotlib.pyplot as plt

from src.report import rows_map

plt.rcParams.update({'font.size': 14, })
# plt.style.use("dark_background")
plt.style.use('seaborn-colorblind')

base_path = "/home/theia/scotthull/"
runs = [
    [base_path + f"Paper1_SPH/gi/500_b073_new/", "Canonical", 200],
    [base_path + f"Paper2_SPH/gi/500_half_earths/", "Half-Earths", 200],
]

min_iteration = 0
max_iteration = 1800
increment = 50
increment_high = 100
increment_low = 20

# new_phase_path = "src/phase_data/forstSTS__vapour_curve.txt"
# old_phase_path = "src/phase_data/duniteN__vapour_curve.txt"
phase_path = "src/phase_data/forstSTS__vapour_curve.txt"


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

def plot_entropy_and_vmf_vs_time():
    sims, titles = [i[0] for i in runs], [i[1] for i in runs]
    fig, axs = plt.subplots(2, 3, figsize=(18, 9), sharex="all")
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    axs = axs.flatten()
    d = {}
    d_wo_circ = {}
    for sim, title in zip(sims, titles):
        s = sim.split("/")[-2]
        if title not in d.keys():
            d.update({title: {'TIME_HRS': [], 'MEAN_DISK_ENTROPY_W_CIRC': [], 'MEAN_DISK_TEMPERATURE': [],
                              'DISK_VMF_W_CIRC': [], "DISK_MASS": [],
                              "DISK_ANGULAR_MOMENTUM": [], 'DISK_THEIA_MASS_FRACTION': []}})
            d_wo_circ.update({title: {'TIME_HRS': [], 'MEAN_DISK_ENTROPY_WITHOUT_CIRC': [], 'DISK_VMF_WITHOUT_CIRC': []}})
        for iteration in np.arange(min_iteration, max_iteration + increment, increment):
            try:
                if not "Paper1" in sim:
                    path = sim + "{}_reports/".format(s)
                else:
                    path = sim + "{}_reports2/".format(s)
                df = pd.read_csv(path + "{}.csv".format(iteration))
                time, avg_disk_entropy, disk_temperature, disk_vmf, disk_mass, disk_am, theia_mass_frac = df['TIME_HRS'][0], df['MEAN_DISK_ENTROPY_W_CIRC'][0], \
                                                                       df['MEAN_DISK_TEMPERATURE'][0], df['DISK_VMF_W_CIRC'][0], df['DISK_MASS'][0], \
                                                                       df['DISK_ANGULAR_MOMENTUM'][0], df['DISK_THEIA_MASS_FRACTION'][0]
                d[title]['TIME_HRS'].append(time)
                d[title]['MEAN_DISK_ENTROPY_W_CIRC'].append(avg_disk_entropy)
                d[title]['DISK_VMF_W_CIRC'].append(disk_vmf)
                d[title]['DISK_MASS'].append(disk_mass)
                d[title]['DISK_ANGULAR_MOMENTUM'].append(disk_am)
                d[title]['MEAN_DISK_TEMPERATURE'].append(disk_temperature)
                d[title]['DISK_THEIA_MASS_FRACTION'].append(theia_mass_frac)
                d_wo_circ[title]['TIME_HRS'].append(time)
                d_wo_circ[title]['MEAN_DISK_ENTROPY_WITHOUT_CIRC'].append(df['MEAN_DISK_ENTROPY_WITHOUT_CIRC'][0])
                d_wo_circ[title]['DISK_VMF_WITHOUT_CIRC'].append(df['DISK_VMF_WITHOUT_CIRC'][0])
            except Exception as e:
                print(e)
                pass
    for i in d.keys():
        for j in d[i].keys():
            d[i][j] = [float(str(k).split(" ")[0]) for k in d[i][j]]
    for index, sim in enumerate(list(d.keys())):
        color = colors[index]
        if "high" in sim or "low" in sim:
            linestyle = "dotted"
        for index, h in enumerate(list(d[sim].keys())[1:]):
            axs[index].plot(
                d[sim]['TIME_HRS'], d[sim][h], linewidth=2.0, color=color
            )
            axs[index].set_ylabel(rows_map[h][1:-1], fontsize=16)
            if h == "MEAN_DISK_ENTROPY_W_CIRC":
                axs[index].plot(
                    d_wo_circ[sim]['TIME_HRS'], d_wo_circ[sim]['MEAN_DISK_ENTROPY_WITHOUT_CIRC'], linewidth=2.0,
                    linestyle="--", color=color
                )
                axs[index].set_ylabel("Avg. Disk Entropy (J/K)", fontsize=16)
            elif h == "DISK_VMF_W_CIRC":
                axs[index].plot(
                    d_wo_circ[sim]['TIME_HRS'], d_wo_circ[sim]['DISK_VMF_WITHOUT_CIRC'], linewidth=2.0,
                    linestyle="--", color=color
                )
    for index, title in enumerate(titles):
        axs[0].plot([], [], linewidth=2.0, linestyle="-", color=colors[index], label=title)
    letters = list(string.ascii_lowercase)
    for index, ax in enumerate(axs):
        x1, x2, y1, y2 = ax.axis()
        x_loc = x1 + (0.02 * (x2 - x1))
        y_loc = y2 - (0.08 * (y2 - y1))
        ax.grid(alpha=0.4)
        ax.text(x_loc, y_loc, letters[index], fontweight="bold")
    for ax in axs[-3:]:
        ax.set_xlabel("Time (hrs)", fontsize=16)

    plt.tight_layout()
    legend = fig.legend(loc=7, fontsize=16)
    for line in legend.get_lines():  # increase line widths in legend
        try:
            line.set_linewidth(4.0)
        except:
            pass
    for handle in legend.legendHandles:  # increase marker sizes in legend
        try:
            handle.set_sizes([120.0])
        except:
            pass
    fig.subplots_adjust(right=0.62)
    plt.savefig("paper2_disk_entropy_and_vmf.png", format='png', dpi=300)


plot_entropy_and_vmf_vs_time()
