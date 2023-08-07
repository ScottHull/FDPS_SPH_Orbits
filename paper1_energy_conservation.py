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

sims, titles = get_all_sims(high=True, angle=angle)
# make a 2 column 1 row plot
fig, axs = plt.subplots(1, 3, figsize=(20, 6))
axs = axs.flatten()

for ax in axs:
    ax.grid()
    ax.set_xlabel("Time (hrs.)")
    # make all axis tick font size larger
    ax.tick_params(axis='both', which='major', labelsize=14)
axs[0].set_ylabel(r"Disk Mass ($M_{\oplus}$)", fontsize=16)
axs[1].set_ylabel(r"Disk Energy (kJ/kg)", fontsize=16)
axs[2].set_ylabel(r"Total Energy (kJ/kg)", fontsize=16)

times_dict = {}
total_potential_energy = {}
total_internal_energy = {}
total_kinetic_energy = {}

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
d = {}
for sim, title in zip(sims, titles):
    times = []
    disk_mass = []
    disk_energy = []
    total_energy = []
    times_dict.update({title: []})
    total_potential_energy.update({title: []})
    total_internal_energy.update({title: []})
    total_kinetic_energy.update({title: []})
    cutoff_density = int(title.split("b")[0])
    color = colors[cutoff_densities.index(cutoff_density)]
    linestyle = "-"
    to_increment = increment
    if "high" in sim:
        to_increment = increment_high
    if "low" in sim:
        to_increment = increment_low
    if "N" in sim:
        linestyle = "--"
    if "high" in sim or "low" in sim:
        linestyle = "dotted"
    for iteration in np.arange(min_iteration, max_iteration + to_increment, to_increment):
        path = base_path + "{}/circularized_{}/".format(sim, sim)
        report_path = base_path + "{}/{}_reports/".format(sim, sim)
        if os.path.exists(report_path + "{}.csv".format(iteration)):
            report_df = pd.read_csv(report_path + "{}.csv".format(iteration))
            time = report_df['TIME_HRS'][0]
            df = pd.read_csv(path + "/{}.csv".format(iteration))
            df['velocity'] = np.sqrt(df['vx']**2 + df['vy']**2 + df['vz']**2)
            disk = df[df['label'] == 'DISK']
            # specific_internal_energy = sum(df['internal_energy'] / df['mass'])
            # specific_potential_energy = sum(df['potential_energy'] / df['mass'])
            # specific_kinetic_energy = sum(0.5 * df['velocity']**2)
            # disk_specific_internal_energy = sum(disk['internal_energy'] / disk['mass'])
            # disk_specific_potential_energy = sum(disk['potential_energy'] / disk['mass'])
            # disk_specific_kinetic_energy = sum(0.5 * disk['velocity']**2)
            # specific_energy_total = (specific_internal_energy + specific_potential_energy + specific_kinetic_energy) / 1000
            # specific_energy_disk = (disk_specific_internal_energy + disk_specific_potential_energy + disk_specific_kinetic_energy) / 1000
            # times_dict[title].append(time)
            specific_internal_energy = sum(df['internal_energy'])
            specific_potential_energy = sum(df['potential_energy'])
            specific_kinetic_energy = sum(0.5 * df['mass'] * df['velocity']**2)
            disk_specific_internal_energy = sum(disk['internal_energy'])
            disk_specific_potential_energy = sum(disk['potential_energy'])
            disk_specific_kinetic_energy = sum(0.5 * disk['mass'] * disk['velocity']**2)
            specific_energy_total = (specific_internal_energy + specific_potential_energy + specific_kinetic_energy) / 1000
            specific_energy_disk = (disk_specific_internal_energy + disk_specific_potential_energy + disk_specific_kinetic_energy) / 1000
            times_dict[title].append(time)

            times.append(time)
            disk_mass.append(float(report_df['DISK_MASS'][0].split(" ")[0]))
            disk_energy.append(specific_energy_disk)
            total_energy.append(specific_energy_total)
            total_potential_energy[title].append(specific_potential_energy / 1000)
            total_internal_energy[title].append(specific_internal_energy / 1000)
            total_kinetic_energy[title].append(specific_kinetic_energy / 1000)

    axs[0].plot(times, disk_mass, color=color, linestyle=linestyle)
    axs[1].plot(times, disk_energy, color=color, linestyle=linestyle)
    axs[2].plot(times, total_energy, color=color, linestyle=linestyle)

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
    ax.text(x_loc, y_loc, letters[index], fontweight="bold")

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


# make a 3 column 1 row plot
fig, axs = plt.subplots(1, 3, figsize=(20, 6))
axs = axs.flatten()
for ax in axs:
    ax.grid()
    ax.set_xlabel("Time (hrs.)")
    # make all axis tick font size larger
    ax.tick_params(axis='both', which='major', labelsize=14)
axs[0].set_ylabel(r"Potential Energy (kJ/kg)", fontsize=16)
axs[1].set_ylabel(r"Internal Energy (kJ/kg)", fontsize=16)
axs[2].set_ylabel(r"Kinetic Energy (kJ/kg)", fontsize=16)

for sim, title in zip(sims, titles):
    cutoff_density = int(title.split("b")[0])
    color = colors[cutoff_densities.index(cutoff_density)]
    linestyle = "-"
    to_increment = increment
    if "high" in sim:
        to_increment = increment_high
    if "low" in sim:
        to_increment = increment_low
    if "N" in sim:
        linestyle = "--"
    if "high" in sim or "low" in sim:
        linestyle = "dotted"
    axs[0].plot(
        times_dict[title], total_potential_energy[title], color=color, linestyle=linestyle
    )
    axs[1].plot(
        times_dict[title], total_internal_energy[title], color=color, linestyle=linestyle
    )
    axs[2].plot(
        times_dict[title], total_kinetic_energy[title], color=color, linestyle=linestyle
    )

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
plt.savefig("energy_conservation_components_{}.png".format(angle), dpi=300)


