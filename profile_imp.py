#!/usr/bin/env python3
import os
import csv
import pandas as pd
import numpy as np
from math import atan, pi, acos
from operator import contains
from random import randint
import multiprocessing as mp
import matplotlib.pyplot as plt

from src.combine import CombineFile


cutoff_densities = [5, 500, 1000, 2000]
number_processes = 200
base_path = "/home/theia/scotthull/Paper1_SPH/tar-imp/"

for i in ['new', 'old']:
    eos = "Stewart M-ANEOS"
    if i == "old":
        eos = "N-SPH M-ANEOS"
    fig, axs = plt.subplots(2, 2, figsize=(10, 10))
    axs = axs.flatten()
    for cd in cutoff_densities:
        name = "{}_{}".format(cd, i)
        f = base_path + "{}/imp.dat".format(name)
        headers = ["id", "tag", "mass", "x", "y", "z", "vx", "vy", "vz", "density", "internal energy", "pressure",
                   "potential energy", "entropy", "temperature"]
        df = pd.read_csv(f, skiprows=2, header=None, delimiter="\t", names=headers)
        df['radius'] = [((i ** 2 + j ** 2 + k ** 2) ** 0.5) / 1000.0 for i, j, k in zip(df['x'], df['y'], df['z'])]
        axs[0].scatter(
            df['radius'], df['density'], marker='.', s=2, label=name
        )
        axs[1].scatter(
            df['radius'], df['pressure'], marker='.', s=2, label=name
        )
        axs[2].scatter(
            df['radius'], df['internal energy'], marker='.', s=2, label=name
        )
        axs[3].scatter(
            df['radius'], df['entropy'], marker='.', s=2, label=name
        )
        print(name, max(df['radius']), min(df['internal energy']), min(df['temperature']), min(df['pressure']))
        if cd == 5 and i == "new":
            name = "{}_{}_high".format(cd, i)
            f = base_path + "{}/imp.dat".format(name)
            headers = ["id", "tag", "mass", "x", "y", "z", "vx", "vy", "vz", "density", "internal energy", "pressure",
                       "potential energy", "entropy", "temperature"]
            df = pd.read_csv(f, skiprows=2, header=None, delimiter="\t", names=headers)
            df['radius'] = [((i ** 2 + j ** 2 + k ** 2) ** 0.5) / 1000.0 for i, j, k in zip(df['x'], df['y'], df['z'])]
            axs[0].scatter(
                df['radius'], df['density'], marker='.', s=2, label=name
            )
            axs[1].scatter(
                df['radius'], df['pressure'], marker='.', s=2, label=name
            )
            axs[2].scatter(
                df['radius'], df['internal energy'], marker='.', s=2, label=name
            )
            axs[3].scatter(
                df['radius'], df['entropy'], marker='.', s=2, label=name
            )
            print(name, max(df['radius']), min(df['internal energy']), min(df['temperature']), min(df['pressure']))

    for ax in axs:
        ax.grid(alpha=0.4)
        ax.set_xlabel("Radius (km)")
        ax.set_title(eos)

    axs[0].set_ylabel("Density (kg/m^3)")
    axs[1].set_ylabel("Pressure (Pa)")
    axs[2].set_ylabel("Internal Energy (J)")
    axs[3].set_ylabel("Entropy (J/kg/K)")
    axs[0].legend(loc='upper left')
    plt.savefig("{}_impactor_profile.png".format(eos), format='png', dpi=200)
