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
import seaborn as sns

from src.animate import animate
from src.combine import CombineFile
from src.identify import ParticleMap
from src.report import get_sim_report, write_report_at_time, build_latex_table_from_disk_report, rows_map

base_path = "/home/theia/scotthull/Paper1_SPH/gi/"
to_path = "sawtooth_gi_visualization"
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
square_scale = 6e7 / 10 ** 7
num_samples = 10

headers = ["id", "tag", "mass", "x", "y", "z", "vx", "vy", "vz", "density", "internal energy", "pressure",
           "potential energy", "entropy", "temperature"]

# use dark background
# plt.style.use("dark_background")

if os.path.exists(to_path):
    shutil.rmtree(to_path)
os.mkdir(to_path)

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


def get_smth_length(m, rho, xi=1.2, n=3):
    """
    h = xi * (m_j / rho_j)^(1/n)
    where xi is a coefficient, m_j is particle mass, and n is the resolution of the simulation.
    Increasing resolution decreases particle mass and therefore decreases smoothing length.
    Decreasing cutoff density increases smoothing length.
    """
    return xi * ((m / rho) ** (1 / n))


def get_endstate(s, iteration, only_disk=True):
    path = base_path + "{}/circularized_{}".format(s, s)
    df = pd.read_csv(path + "/{}.csv".format(iteration))
    if only_disk:
        return df[df['label'] == "DISK"]
    else:
        return df


# make a 3 column plot of density, entropy, internal energy
ylabels = [r'Density (kg/m$^3$)', r'Entropy (J/kg/K)', r'Internal Energy (kJ)']

# get a spectrum of colors in the jet colormap that equals the number of samples
colors = sns.color_palette('husl', n_colors=len(range(num_samples)))  # a list of RGB tuples

names, titles = get_all_sims()

# time = {j: {i: [] for i in endstate['id'].values} for j in titles}
# entropy = {i: [] for i in endstate['id'].values}
# internal_energy = {i: [] for i in endstate['id'].values}
# density = {i: [] for i in endstate['id'].values}
# temperature = {i: [] for i in endstate['id'].values}
# sel_time = {i: [] for i in select_particles}
# sel_entropy = {i: [] for i in select_particles}
# sel_internal_energy = {i: [] for i in select_particles}
# sel_density = {i: [] for i in select_particles}
# sel_temperature = {i: [] for i in select_particles}

for s, t in zip(names, titles):
    if os.path.exists(f"paper1_sawtooth_{s}"):
        shutil.rmtree(f"paper1_sawtooth_{s}")
    os.mkdir(f"paper1_sawtooth_{s}")
    # get the endstate df
    endstate = get_endstate(s, endstate_iteration)
    midstate = get_endstate(s, mid_iteration, only_disk=False)
    # get a list of the difference in entropy between the endstate and midstate
    select_particles = [i for i in endstate['id'].tolist() if endstate[endstate['id'] == i]['entropy'].values - midstate[midstate['id'] == i]['entropy'].values > 500]
    # get 5 random particles from the select_particles list
    select_particles = np.random.choice(select_particles, num_samples)
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
    sel_neighbors = {i: [] for i in select_particles}
    for iteration in np.arange(min_iteration, max_iteration + increment, increment):
        print(f"Working on {s} at iteration {iteration}...")
        path = base_path + "{}/{}".format(s, s)
        to_fname = "merged_{}_{}.dat".format(iteration, randint(0, 100000))
        cf = CombineFile(num_processes=number_processes, time=iteration, output_path=path, to_fname=to_fname)
        df = cf.combine_df()
        formatted_time = round(cf.sim_time * 0.000277778, 2)
        df.columns = headers
        # for i in endstate['id'].values:
        #     time[i].append(formatted_time)
        #     entropy[i].append(df[df['id'] == i]['entropy'].values[0])
        #     internal_energy[i].append(df[df['id'] == i]['internal energy'].values[0])
        #     density[i].append(df[df['id'] == i]['density'].values[0])
        #     temperature[i].append(df[df['id'] == i]['temperature'].values[0])
        for i in select_particles:
            sel_time[i].append(formatted_time)
            sel_entropy[i].append(df[df['id'] == i]['entropy'].values[0])
            sel_internal_energy[i].append(df[df['id'] == i]['internal energy'].values[0])
            sel_density[i].append(df[df['id'] == i]['density'].values[0])
            sel_temperature[i].append(df[df['id'] == i]['temperature'].values[0])

            smth = get_smth_length(df[df['id'] == i]['mass'].values[0], df[df['id'] == i]['density'].values[0])
            # get the number of particles in the dataframe whose distance from the sample particle is less than the smoothing length
            neighbors = len(df[(df['x'] - df[df['id'] == i]['x'].values[0]) ** 2 +
                               (df['y'] - df[df['id'] == i]['y'].values[0]) ** 2 +
                               (df['z'] - df[df['id'] == i]['z'].values[0]) ** 2 <= smth ** 2]) - 1
            sel_neighbors[i].append(neighbors)

        if iteration % 20 == 0:
            fig = plt.figure(figsize=(10, 10))
            ax = fig.add_subplot(111)
            ax.set_xlim(-square_scale, square_scale)
            ax.set_ylim(-square_scale, square_scale)
            ax.set_xticks([], minor=False)
            ax.set_yticks([], minor=False)
            ax.scatter(df['x'] / 10 ** 7, df['y'] / 10 ** 7, s=0.8, marker=".")
            for index, i in enumerate(select_particles):
                ax.scatter(
                    df[df['id'] == i]['x'] / 10 ** 7, df[df['id'] == i]['y'] / 10 ** 7, s=200, marker=".", alpha=1,
                    color=colors[index]
                )
            ax.set_title(f"{formatted_time} hrs.")
            plt.tight_layout()
            plt.savefig(f"paper1_sawtooth_{s}/{iteration}.png", dpi=200)

    plt.rcParams.update({'font.size': 16, })
    # plt.style.use("dark_background")
    plt.style.use('seaborn-colorblind')

    fig, axs = plt.subplots(1, 5, figsize=(30, 6), sharex='all')
    axs = axs.flatten()
    ylabels = [r'Density (kg/m$^3$)', r'Entropy (J/kg/K)', r'Internal Energy (kJ/kg)', r'Temperature (K)', 'Neighbors']
    for ax, y in zip(axs, ylabels):
        ax.grid()
        ax.set_xlabel("Time (hours)", fontsize=16)
        ax.set_ylabel(y, fontsize=16)

    # for i in endstate['id'].values:
    #     axs[0].plot(time[i], density[i], alpha=0.05)
    #     axs[1].plot(time[i], entropy[i], alpha=0.05)
    #     axs[2].plot(time[i], internal_energy[i], alpha=0.05)
    #     axs[3].plot(time[i], temperature[i], alpha=0.05)
    for index, i in enumerate(select_particles):
        axs[0].plot(sel_time[i], sel_density[i], alpha=1, color=colors[index])
        axs[1].plot(sel_time[i], sel_entropy[i], alpha=1, color=colors[index])
        axs[2].plot(sel_time[i], sel_internal_energy[i], alpha=1, color=colors[index])
        axs[3].plot(sel_time[i], sel_temperature[i], alpha=1, color=colors[index])
        axs[4].plot(sel_time[i], sel_neighbors[i], alpha=1, color=colors[index])
    axs[0].set_ylim(0, 50)
    axs[4].set_ylim(0, 10)

    plt.tight_layout()
    plt.savefig(f"{s}_sawtooth_entropy_all_disk.png", dpi=300)

    animate(
        start_time=min_iteration,
        end_time=max_iteration,
        interval=increment,
        path=f"paper1_sawtooth_{s}",
        fps=10,
        filename=f"paper1_sawtooth_{s}.mp4"
    )



