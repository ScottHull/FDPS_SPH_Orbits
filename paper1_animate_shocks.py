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
        disk = df[df['id'].isin(endstate)]
        disk_sample = df[df['id'].isin(endstate_sample['id'])]
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
            disk_sample['entropy_change'] = disk_sample['entropy'] - prev_df['entropy']
            disk_samples[p]['samples'].append(disk_sample)