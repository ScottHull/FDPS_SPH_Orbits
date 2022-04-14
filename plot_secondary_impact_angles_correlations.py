#!/usr/bin/env python3
import os
import csv
import pandas as pd
import numpy as np
from random import randint
import multiprocessing as mp
import matplotlib.pyplot as plt

from src.report import rows_map

plt.rcParams.update({'font.size': 18, })
plt.style.use("dark_background")

base_path = "/home/theia/scotthull/Paper1_SPH/gi/"
angle = "b073"
end_iteration = 1800
cutoff_densities = [5, 500, 1000, 2000]

secondary_impact_times = {
    '5b073n': {
        'iteration': 220,
        'time': 6.11,
    },
    '500b073n': {
        'iteration': 255,
        'time': 7.08,
    },
    '1000b073n': {
        'iteration': 235,
        'time': 6.53,
    },
    '2000b073n': {
        'iteration': 205,
        'time': 5.69,
    },
    '5b073o': {
        'iteration': 235,
        'time': 6.53,
    },
    '500b073o': {
        'iteration': 260,
        'time': 7.22,
    },
    '1000b073o': {
        'iteration': 265,
        'time': 7.36,
    },
    '2000b073o': {
        'iteration': 245,
        'time': 6.81,
    },
}

def plot_impact_angles_vs_time():
    df = pd.read_csv("/Users/scotthull/Desktop/b073_secondary_impact_angles.csv", index_col="Iteration")
    df = df[df.index < 300]
    fig, axs = plt.subplots(1, 2, figsize=(16, 9), sharex="all", sharey='all')
    axs = axs.flatten()
    new_index, old_index = 0, 1
    for run in df.keys():
        to_index = new_index
        if "o" in run:
            to_index = old_index
        axs[to_index].plot(
            df.index, df[run], linewidth=2.0, label=run
        )
    for ax in axs:
        ax.grid(alpha=0.4)
        ax.set_xlabel("Iteration")
        ax.legend()
    axs[0].set_ylabel("Secondary Impact Angle (deg.)")
    axs[0].set_title("Stewart M-ANEOS")
    axs[1].set_title("GADGET M-ANEOS")
    plt.show()


def plot_moment_of_impact_vs_angle():
    df = pd.read_csv("/Users/scotthull/Desktop/b073_secondary_impact_angles.csv", index_col="Iteration")
    fig, axs = plt.subplots(1, 2, figsize=(16, 9), sharex="all", sharey='all')
    axs = axs.flatten()
    new_index, old_index = 0, 1
    for run in df.keys():
        to_index = new_index
        if "o" in run:
            to_index = old_index
        iteration, time = secondary_impact_times[run]['iteration'], secondary_impact_times[run]['time']
        angle = df[run][iteration]
        secondary_impact_times[run].update({'angle': angle})
        axs[to_index].scatter(
            time, angle, s=80, label=run
        )
    for ax in axs:
        ax.grid(alpha=0.4)
        ax.set_xlabel("Time (hrs)")
        ax.legend()
    axs[0].set_ylabel("Secondary Impact Angle (deg.)")
    axs[0].set_title("Stewart M-ANEOS")
    axs[1].set_title("GADGET M-ANEOS")
    plt.show()


# plot_impact_angles_vs_time()
# plot_moment_of_impact_vs_angle()

def get_all_sims(high=True):
    fformat = "{}_{}_{}"
    tformat = "{}{}{}"
    names = []
    titles = []
    for runs in ["new", "old"]:
        n = "n"
        if runs == "old":
            n = "o"
        for cd in cutoff_densities:
            output_name = fformat.format(cd, angle, runs)
            title_name = tformat.format(cd, angle, n)
            titles.append(title_name)
            names.append(output_name)
            if cd == 5 and high and runs == "new":
                output_name = fformat.format(cd, angle, runs) + "_high"
                names.append(output_name)
                title_name = tformat.format(cd, angle, n) + "-high"
                titles.append(title_name)
    return names, titles


def plot_vs_disk_property(r_dot_v: bool):
    angles_path = "{}_secondary_impact_angles.csv".format(angle)
    r_dot_v_path = "{}_secondary_impact_angles_r_dot_v.csv".format(angle)
    angles_df = pd.read_csv(angles_path, index_col="Unnamed: 0")
    r_dot_v_df = pd.read_csv(r_dot_v_path, index_col="Unnamed: 0")
    x = angles_df
    x_label = "Impact Angle (deg.)"
    if r_dot_v:
        x = r_dot_v_df
        x_label = "r $\cdot$ v"
    sims, titles = get_all_sims(high=False)
    points = ["DISK_MASS", "DISK_ANGULAR_MOMENTUM", "MEAN_DISK_ENTROPY_W_CIRC", "DISK_VMF_WITH_CIRC"]
    fig, axs = plt.subplots(2, 2, figsize=(16, 9), gridspec_kw={"hspace": 0.14, "wspace": 0.16})
    axs = axs.flatten()
    for s, t in zip(sims, titles):
        impact_point = secondary_impact_times[t]
        report_path = base_path + "{}/{}_reports/{}.csv".format(s, s, end_iteration)
        report = pd.read_csv(report_path)
        for index, p in enumerate(points):
            axs[index].scatter(
                x[t][impact_point["iteration"]], float(report[p][0].split(" ")[0]), s=80, label=t
            )
            axs[index].set_ylabel(rows_map[p][1:-1])
    for ax in axs:
        ax.grid(alpha=0.4)
        ax.legend()
    for ax in axs[-2:]:
        ax.set_xlabel(x_label)
    f = "angles"
    if r_dot_v:
        f = "r_dot_v"
    plt.savefig("{}_{}_secondary_impact_vs_disk_property.png".format(angle, f), format='png', dpi=200)

plot_vs_disk_property(r_dot_v=False)
plot_vs_disk_property(r_dot_v=True)
