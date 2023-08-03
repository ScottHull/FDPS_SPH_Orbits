##!/usr/bin/env python3
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
cutoff_densities = [5, 500, 1000, 2000]
min_iteration = 0
max_iteration = 1800
increment = 50
increment_high = 100
increment_low = 20
number_processes = 200
number_processes_high = 500
earth_mass = 5.972e24  # kg


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
    return names, titles


times = []
disk_mass = []
disk_energy = []
total_energy = []

sims, titles = get_all_sims(high=True, angle=angle)
# make a 2 column 1 row plot
fig, axs = plt.subplots(1, 3, figsize=(18, 6))
axs = axs.flatten()

for ax in axs:
    ax.grid()
    ax.set_xlabel("Time (hrs.)")
axs[0].set_ylabel(r"Disk Mass ($M_{\oplus}$)")
axs[1].set_ylabel(r"Disk Energy (kJ/kg)")
axs[2].set_ylabel(r"Total Energy (kJ/kg)")

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
d = {}
for sim, title in zip(sims, titles):
    cutoff_density = int(title.split("b")[0])
    color = colors[cutoff_densities.index(cutoff_density)]
    linestyle = "-"
    if "N" in sim:
        linestyle = "--"
    if "high" in sim or "low" in sim:
        linestyle = "dotted"
    for iteration in np.arange(min_iteration, max_iteration + increment, increment):
        path = base_path + "{}/circularized_{}/".format(sim, sim)
        report_path = base_path + "{}/{}_reports/".format(sim, sim)
        report_df = pd.read_csv(report_path + "{}.csv".format(iteration))
        time = report_df['TIME_HRS'][0]
        df = pd.read_csv(path + "/{}.csv".format(iteration))
        df['velocity'] = np.sqrt(df['vx']**2 + df['vy']**2 + df['vz']**2)
        disk = df[df['label'] == 'DISK']
        specific_energy_total = sum(df['internal_energy']) + sum(df['potential_energy']) + sum((0.5 * df['velocity']**2)) / 1000
        specific_energy_disk = sum(disk['internal_energy']) + sum(disk['potential_energy']) + sum((0.5 * disk['velocity']**2)) / 1000

        times.append(time)
        disk_mass.append(report_df['DISK_MASS'][0])
        disk_energy.append(specific_energy_disk)
        total_energy.append(specific_energy_total)


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
plt.savefig("energy_conservation_{}.png".format(angle), dpi=300)
