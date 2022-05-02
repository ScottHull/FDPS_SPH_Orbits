#!/usr/bin/env python3
import os
import csv
import pandas as pd
import numpy as np
from random import randint
import string
import multiprocessing as mp
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerLine2D

from src.report import rows_map

plt.rcParams.update({'font.size': 14, })
# plt.style.use("dark_background")
plt.style.use('seaborn-colorblind')

# base_path = "/home/theia/scotthull/Paper1_SPH/gi/"
base_path = "/Users/scotthull/Desktop/"
end_iteration = 1800
cutoff_densities = [5, 500, 1000, 2000]


secondary_impact_times = {
    '5b073n': {
        'iteration': 220,
        'time': 6.11,
        'primary_impact_iteration': 15,
        'primary_impact_time': 0.42,
        'characteristic_tail_iteration': None,
        'characteristic_tail_time': None,
    },
    '500b073n': {
        'iteration': 255,
        'time': 6.94,
        'primary_impact_iteration': 15,
        'primary_impact_time': 0.42,
        'characteristic_tail_iteration': None,
        'characteristic_tail_time': None,
    },
    '1000b073n': {
        'iteration': 235,
        'time': 6.53,
        'primary_impact_iteration': 15,
        'primary_impact_time': 0.42,
        'characteristic_tail_iteration': None,
        'characteristic_tail_time': None,
    },
    '2000b073n': {
        'iteration': 205,
        'time': 5.69,
        'primary_impact_iteration': 15,
        'primary_impact_time': 0.42,
        'characteristic_tail_iteration': None,
        'characteristic_tail_time': None,
    },
    '5b073n-high': {
        'iteration': 220,
        'time': 6.11,
        'primary_impact_iteration': 15,
        'primary_impact_time': 0.42,
        'characteristic_tail_iteration': None,
        'characteristic_tail_time': None,
    },
    '5b073o': {
        'iteration': 235,
        'time': 6.53,
        'primary_impact_iteration': 15,
        'primary_impact_time': 0.42,
        'characteristic_tail_iteration': None,
        'characteristic_tail_time': None,
    },
    '500b073o': {
        'iteration': 260,
        'time': 7.22,
        'primary_impact_iteration': 15,
        'primary_impact_time': 0.42,
        'characteristic_tail_iteration': None,
        'characteristic_tail_time': None,
    },
    '1000b073o': {
        'iteration': 265,
        'time': 7.36,
        'primary_impact_iteration': 15,
        'primary_impact_time': 0.42,
        'characteristic_tail_iteration': None,
        'characteristic_tail_time': None,
    },
    '2000b073o': {
        'iteration': 245,
        'time': 6.81,
        'primary_impact_iteration': 15,
        'primary_impact_time': 0.42,
        'characteristic_tail_iteration': None,
        'characteristic_tail_time': None,
    },
}


def plot_secondary_impact_angle_with_struct():
    r_dot_v_angle_path = "b073_secondary_impact_angles_angle_between_r_and_v.csv"
    si_and_tail_path = "b073_secondary_impact_struc_data.csv"
    r_dot_v_angle_df = pd.read_csv(base_path + r_dot_v_angle_path, index_col="Unnamed: 0")
    si_and_tail_df = pd.read_csv(base_path + si_and_tail_path, index_col="Unnamed: 0")
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    runs = si_and_tail_df.keys()

    fig, axs = plt.subplots(2, 2, figsize=(16, 9), sharex="all", sharey='all')
    axs = axs.flatten()
    for ax in axs:
        ax.grid(alpha=0.4)
    headers = ["% DISK FROM SI", "% DISK FROM TAIL", "% DISK FROM INSIDE PLANET", "% DISK FROM OUTSIDE PLANET (not tail/si)"]

    to_index = 0
    for h in headers:
        for r in runs:
            marker = "o"
            if "o" in r:
                marker = "^"
            cd = int(r.split("b")[0])
            color = colors[cutoff_densities.index(cd)]

            angle = r_dot_v_angle_df[r][secondary_impact_times[r]['iteration']]
            si_and_tail_disk_fracs = si_and_tail_df[r]
            axs[to_index].scatter(180 - angle, si_and_tail_disk_fracs[h], s=80, marker=marker, color=color)
            axs[to_index].set_xlabel(r"$\phi$")
            axs[to_index].set_ylabel("% Disk")
        axs[to_index].set_title(h)
        to_index += 1
    for cd in cutoff_densities:
        axs[0].scatter([], [], s=80, marker="s", color=colors[cutoff_densities.index(cd)], label=r"$\rho_c$ = {} kg/m$^3$".format(cd))
    for i in ["Stewart M-ANEOS", "GADGET M-ANEOS"]:
        marker = "o"
        if "GADGET" in i:
            marker = "^"
        axs[0].scatter([], [], s=80, marker=marker, color='black', label=i)
    axs[0].legend(loc='upper left')
    plt.savefig("b073_secondary_impact_angle_with_struct.png", format='png', dpi=200)
    # plt.show()

plot_secondary_impact_angle_with_struct()

