#!/usr/bin/env python3
import os
import shutil
import csv
import pandas as pd
import numpy as np
from math import pi
import string
from random import randint
import multiprocessing as mp
from matplotlib.colors import Normalize
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from src.animate import animate
from src.combine import CombineFile

base_path = "/home/theia/scotthull/Paper1_SPH/gi/"
angle = "b073"
paths = [
    "5_b073_new", "500_b073_new",
]

cutoff_densities = [5, 500]
density_ranges = [[0, 20], [495, 520]]
min_iteration = 0
max_iteration = 1800
endstate_iteration = max_iteration
increment = 1
number_processes = 200
norm_min = 5
norm_max = 500

normalizer = Normalize(norm_min, norm_max)
cmap = cm.get_cmap("jet")

plt.rcParams.update({'font.size': 16, })
# plt.style.use("dark_background")
plt.style.use('seaborn-colorblind')

def get_smth_length(m, rho, xi=1.2, n=3):
    """
    h = xi * (m_j / rho_j)^(1/n)
    where xi is a coefficient, m_j is particle mass, and n is the resolution of the simulation.
    Increasing resolution decreases particle mass and therefore decreases smoothing length.
    Decreasing cutoff density increases smoothing length.
    """
    return xi * ((m / rho) ** (1 / n))

def wendland_c6(r, h, support_radius=2.5, n=3):
    h *= support_radius  # add support radius
    sigma = 1 / (h ** n)
    s = abs(r) / h
    W = (1365 * sigma) / (64 * pi)
    if s < 1:
        W *= ((1 - s) ** 8) * (1 + (8 * s) + (25 * (s ** 2)) + (32 * (s ** 3)))
    else:
        W *= 0
    return W

def get_endstate(s):
    path = base_path + "{}/circularized_{}".format(s, s)
    df = pd.read_csv(path + "/{}.csv".format(endstate_iteration))
    return df[df['label'] == "DISK"]

disk_samples = {}

for index, p in enumerate(paths):
    if os.path.exists(f"shocks_{p}"):
        shutil.rmtree(f"shocks_{p}")
    os.mkdir(f"shocks_{p}")
    disk_samples.update({p: {
        'times': [],
        'samples': [],
    }})
    endstate = get_endstate(p)
    endstate_sample = endstate[endstate['density'] == cutoff_densities[index]]
    endstate_sample = endstate_sample[endstate_sample['entropy'] > 9000].sample(10)
    # append another sample with entropy < 9000
    endstate_sample = endstate_sample.append(endstate[endstate['entropy'] < 9000].sample(10))
    cd = cutoff_densities[index]
    prev_df = None
    prev_df_disk_sample = None
    for iteration in np.arange(min_iteration, max_iteration + increment, increment):
        path = base_path + "{}/{}".format(p, p)
        to_fname = "merged_{}_{}.dat".format(iteration, randint(0, 100000))
        cf = CombineFile(num_processes=number_processes, time=iteration, output_path=path, to_fname=to_fname)
        combined_file = cf.combine()
        formatted_time = round(cf.sim_time * 0.000277778, 2)
        disk_samples[p]['times'].append(formatted_time)
        f = os.getcwd() + "/{}".format(to_fname)
        headers = ["id", "tag", "mass", "x", "y", "z", "vx", "vy", "vz", "density", "internal energy", "pressure",
                   "potential energy", "entropy", "temperature"]
        df = pd.read_csv(f, skiprows=2, header=None, delimiter="\t", names=headers, index_col='id')
        os.remove(f)
        disk = df[df.index.isin(endstate['id'])]
        disk_sample = df[df.index.isin(endstate_sample['id'])]
        # calculate the smoothing length for each particle in the sample group
        disk_sample['smth'] = get_smth_length(disk_sample['mass'], disk_sample['density'])
        # get the number of particles with distances less than the smoothing length
        for i, row in disk_sample.iterrows():
            disk_sample.loc[i, 'neighbors'] = len(df[(df['x'] - row['x']) ** 2 + (df['y'] - row['y']) ** 2 +
                                                     (df['z'] - row['z']) ** 2 <= row['smth'] ** 2]) - 1
            # disk_sample.loc[i, 'neighbors'] = len(disk_sample[(disk_sample['x'] - row['x']) ** 2 +
            #                                                   (disk_sample['y'] - row['y']) ** 2 +
            #                                                   (disk_sample['z'] - row['z']) ** 2 <= row['smth'] ** 2])
        # calculate the change in entropy for each particle in the sample group
        if prev_df is not None:
            disk['entropy_change'] = disk['entropy'] - prev_df['entropy']
            disk_sample['entropy_change'] = disk_sample['entropy'] - prev_df_disk_sample['entropy']
            disk_samples[p]['samples'].append(disk_sample)

        fig, axs = plt.subplots(1, 4, figsize=(24, 6), sharex='all')
        axs = axs.flatten()
        for var, ylabel, ax in zip(['entropy', 'internal energy', 'temperature'],
                                   [r'Entropy (J/kg/K)', r'Internal Energy (kJ)', r'Temperature (K)'],
                                   axs[:-1]):
            ax.grid()
            ax.set_ylabel(ylabel)
            ax.set_xlabel(r'Density (kg/m$^3$)')
            ax.set_xlim(density_ranges[index][0], density_ranges[index][1])
            ax.set_title(f"{formatted_time} hrs.")
            if prev_df is None:
                ax.scatter(
                    disk['density'], disk[var], s=2, alpha=1, color='k'
                )
            else:
                ax.scatter(
                    disk['density'], disk[var], s=2, alpha=1, color=cmap(normalizer(disk['entropy_change']))
                )
            ax.scatter(
                disk_sample['density'], disk_sample[var], s=60, alpha=1, color='red'
            )
        axs[-1].scatter(
            disk_sample['density'], disk_sample['neighbors'], s=60, alpha=1, color='red'
        )
        axs[-1].set_ylabel(r'Neighbors')
        axs[-1].set_xlabel(r'Density (kg/m$^3$)')
        axs[-1].set_xlim(density_ranges[index][0], density_ranges[index][1])
        axs[-1].set_title(f"{formatted_time} hrs.")
        axs[-1].grid()
        plt.tight_layout()
        plt.savefig(f"shocks_{p}/{iteration}.png")

        prev_df = disk
        prev_df_disk_sample = disk_sample

    animate(
        start_time=min_iteration,
        end_time=max_iteration,
        interval=increment,
        path=f"shocks_{p}",
        fps=80,
        filename=f"{p}_shocks.mp4",
    )

    fig, axs = plt.subplots(1, 5, figsize=(30, 6), sharex='all')
    axs = axs.flatten()
    for ax in axs:
        ax.grid()
        ax.set_xlabel("Time (hours)", fontsize=16)
    ylabels = [r'Density (kg/m$^3$)', r'Entropy (J/kg/K)', r'Internal Energy (kJ)', r'Temperature (K)', r'Neighbors']
    y = ['density', 'entropy', 'internal energy', 'temperature', 'neighbors']
    for ax, ylabel, y in zip(axs, ylabels, y):
        ax.set_ylabel(ylabel, fontsize=16)
        ax.grid()
        ax.set_xlim(0, 50)
        ax.set_xlabel("Time (hrs.)", fontsize=16)
        for particle in disk_samples[p]['samples'][0].index:
            ax.plot(
                disk_samples[p]['times'], [disk_samples[p]['samples'][i][particle][y] for i in range(len(disk_samples[p]['times']))], alpha=1
            )
    plt.tight_layout()
    plt.savefig(f"{p}_shocks.png", dpi=300)

