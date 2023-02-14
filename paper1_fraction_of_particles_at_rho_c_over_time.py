#!/usr/bin/env python3
import os
import csv
import json
import shutil
from math import pi, asin, isnan, exp
import numpy as np
import pandas as pd
from random import randint
from statistics import mean
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import multiprocessing as mp
from matplotlib.ticker import MaxNLocator
import string

from src.vapor import calc_vapor_mass_fraction_from_formatted
from src.geometry import get_impact_geometry_from_formatted, get_velocity_profile_from_formatted
from src.animate import animate
from src.identify import ParticleMap
from src.combine import CombineFile

plt.rcParams.update({'font.size': 14, })
# plt.style.use("dark_background")
plt.style.use('seaborn-colorblind')

base_path = "/home/theia/scotthull/Paper1_SPH/gi/"
cutoff_densities = [5, 500, 1000, 2000]
angles = ['b073', 'b075']
min_iteration = 0
max_iteration = 1800
end_iteration = 1800
increment = 50
number_processes = 200
number_processes_high = 500
new_run = False

new_phase_path = "src/phase_data/forstSTS__vapour_curve.txt"
old_phase_path = "src/phase_data/duniteN__vapour_curve.txt"

output_columns = [
            'id', 'tag', 'mass', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'density', 'internal_energy', 'pressure',
            'potential_energy', 'entropy', 'temperature'
        ]

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


def get_all_sims():
    fformat = "{}_{}_{}"
    tformat = "{}{}{}"
    names = []
    titles = []
    for angle in angles:
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
                if cd == 5 and runs == "new" and angle == "b073":
                    high_res_name = fformat.format(cd, angle, runs) + "_high"
                    high_res_title = tformat.format(cd, angle, n) + "-high"
                if cd == 2000 and runs == "old" and angle == "b075":
                    high_res_name = fformat.format(cd, angle, runs) + "_low"
                    high_res_title = tformat.format(cd, angle, n) + "-low"
            if high_res_name is not None and high_res_title is not None:
                names.append(high_res_name)
                titles.append(high_res_title)
    return names, titles


def get_end_states(s, end_iteration):
    end_state_file = base_path + "{}/circularized_{}/{}.csv".format(s, s, end_iteration)
    return pd.read_csv(end_state_file, index_col="id")


if new_run:
    data = {}
    sims, titles = get_all_sims()
    # colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    for s, t in zip(sims, titles):
        data[t] = {}
        num_proc = number_processes
        if "high" in s:
            num_proc = number_processes_high
        endstate = get_end_states(s, end_iteration)
        for iteration in np.arange(min_iteration, max_iteration + increment, increment):
            cd = cutoff_densities.index(int(s.split("_")[0]))
            path = base_path + "{}/{}".format(s, s)
            to_fname = "merged_{}_{}.dat".format(iteration, randint(0, 100000))
            cf = CombineFile(num_processes=num_proc, time=iteration, output_path=path, to_fname=to_fname)
            combined_file = cf.combine()
            formatted_time = round(cf.sim_time * 0.000277778, 2)
            data[t][formatted_time] = 0.0
            f = os.getcwd() + "/{}".format(to_fname)
            df = pd.read_csv(to_fname, skiprows=2, header=None, delimiter="\t")
            os.remove(to_fname)
            df.columns = output_columns
            endstate_disk = endstate[endstate['label'] == "DISK"]
            disk_particles = df[df['id'].isin(endstate_disk.index.values)]
            fraction_of_particles_at_rho_cutoff = len(disk_particles[disk_particles['density'] ==
                                                                     int(s.split("_")[0])]) / len(disk_particles)
            data[t][formatted_time] = fraction_of_particles_at_rho_cutoff

    # write data to file
    if os.path.exists("paper1_rho_c_fraction_data.json"):
        os.remove("paper1_rho_c_fraction_data.json")
    with open("paper1_rho_c_fraction_data.json", "w") as outfile:
        outfile.write(str({k: v for k, v in data.items()}))
else:
    # read data from file
    with open("paper1_rho_c_fraction_data.json", "r") as infile:
        data = eval(infile.read())

fig, axs = plt.subplots(1, 2, figsize=(16, 9), sharex="all", sharey="all")
axs = axs.flatten()
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

for t in data.keys():
    cd = int(cutoff_densities.index(int(t.split("b")[0])))
    linestyle = "-"
    index = 0
    if "N" in t:
        linestyle = "--"
    if "high" in t or "low" in t:
        linestyle = "dotted"
    if "b075" in t:
        index = 1
    axs[index].plot(
        list(data[t].keys()), list(data[t].values()), linewidth=2.0, color=colors[cd], linestyle=linestyle
    )

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
axs[0].plot(
    [], [], linewidth=2.0, color='black',
    linestyle="dotted", label="5b073S-high or\n2000b075N-low"
)

for ax in axs:
    ax.set_xlabel("Time (hours)")
    ax.grid(alpha=0.8)

axs[0].set_ylabel(r"Particle Fraction at $\rho_c$")
axs[0].set_title(r"b = 0.73")
axs[1].set_title(r"b = 0.75")

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
fig.subplots_adjust(right=0.80)

plt.savefig("particle_fraction_at_rho_cutoff.png", dpi=300)
