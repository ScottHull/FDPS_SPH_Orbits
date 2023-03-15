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

base_path = "/home/theia/scotthull/Paper1_SPH/gi/"
angle = "b073"
cutoff_densities = [5, 100, 500, 1000, 2000]
additional = [
    "5_b073_new_h_max",
    "5b073-h-max",
]
min_iteration = 0
max_iteration = 1800
increment = 50
increment_high = 100
increment_low = 20
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


def get_all_sims(angle, high=True):
    fformat = "{}_{}_{}"
    tformat = "{}{}{}"
    names = []
    titles = []
    for runs in ["new", "old"]:
        high_res_name = None
        high_res_title = None
        n = "S"
        if runs == "old":
            n = "N"
        for cd in cutoff_densities:
            output_name = fformat.format(cd, angle, runs)
            title_name = tformat.format(cd, angle, n)
            titles.append(title_name)
            names.append(output_name)
            if cd == 5 and high and runs == "new" and angle == "b073":
                high_res_name = fformat.format(cd, angle, runs) + "_high"
                high_res_title = tformat.format(cd, angle, n) + "-high"
            if cd == 2000 and high and runs == "old" and angle == "b075":
                high_res_name = fformat.format(cd, angle, runs) + "_low"
                high_res_title = tformat.format(cd, angle, n) + "-low"
        if high_res_name is not None and high_res_title is not None:
            names.append(high_res_name)
            titles.append(high_res_title)
    for i in additional:
        names.append(i[0])
        titles.append(i[1])
    return names, titles



def plot_entropy_and_vmf_vs_time():
    sims, titles = get_all_sims(high=True, angle=angle)
    fig, axs = plt.subplots(2, 3, figsize=(18, 9), sharex="all")
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    axs = axs.flatten()
    d = {}
    for sim, title in zip(sims, titles):
        if title not in d.keys():
            d.update({title: {'TIME_HRS': [], 'MEAN_DISK_ENTROPY_W_CIRC': [], 'MEAN_DISK_TEMPERATURE': [],
                              'DISK_VMF_W_CIRC': [], "DISK_MASS": [],
                              "DISK_ANGULAR_MOMENTUM": [], 'DISK_THEIA_MASS_FRACTION': []}})
        inc = increment
        if "high" in sim:
            inc = increment_high
        elif "low" in sim:
            inc = increment_low
        for iteration in np.arange(min_iteration, max_iteration + inc, inc):
            print(sim, title, iteration)
            try:
                path = base_path + "{}/{}_reports/".format(sim, sim)
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
            except Exception as e:
                print(e)
                pass
    for i in d.keys():
        for j in d[i].keys():
            d[i][j] = [float(str(k).split(" ")[0]) for k in d[i][j]]
    for sim in list(d.keys()):
        cutoff_density = int(sim.split("b")[0])
        color = colors[cutoff_densities.index(cutoff_density)]
        linestyle = "-"
        if "N" in sim:
            linestyle = "--"
        if "high" in sim or "low" in sim:
            linestyle = "dotted"
        if "h-max" in sim:
            linestyle = 'dashdot'
        for index, h in enumerate(list(d[sim].keys())[1:]):
            axs[index].plot(
                d[sim]['TIME_HRS'], d[sim][h], linewidth=2.0, linestyle=linestyle, color=color,
            )
            axs[index].set_ylabel(rows_map[h][1:-1], fontsize=16)
    for c in cutoff_densities:
        axs[0].scatter(
            [], [], marker="s", s=80, label=r"$\rho_c$ = {} kg/m$^3$".format(c)
        )
    axs[0].plot(
        [], [], c='black', linewidth=2.0, linestyle="-", label="Stewart M-ANEOS"
    )
    axs[0].plot(
        [], [], c='black', linewidth=2.0, linestyle="--", label="N-SPH M-ANEOS"
    )
    if angle == "b073":
        axs[0].plot(
            [], [], c=colors[cutoff_densities.index(5)], linewidth=2.0, linestyle="dotted", label="5b073S-high"
        )
    if angle == "b075":
        axs[0].plot(
            [], [], c=colors[cutoff_densities.index(2000)], linewidth=2.0, linestyle="dotted", label="2000b075N-low"
        )
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
    fig.subplots_adjust(right=0.82)
    plt.savefig("{}_disk_entropy_and_vmf.png".format(angle), format='png', dpi=200)


plot_entropy_and_vmf_vs_time()
