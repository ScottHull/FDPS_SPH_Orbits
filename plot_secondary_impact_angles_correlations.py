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

base_path = "/home/theia/scotthull/Paper1_SPH/gi/"
angle = "b073"
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

def change_legend_marker(handle, original):
    """ Change the marker style of the legend handles """
    handle.update_from(original)
    handle.set_marker('s')

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
    # x_label = "Impact Angle (deg.)"
    x_label = r"$\theta$"
    if r_dot_v:
        x = r_dot_v_df
        x_label = "r $\cdot$ v"
    sims, titles = get_all_sims(high=False)
    points = ["DISK_MASS", "DISK_ANGULAR_MOMENTUM"]
    fig, axs = plt.subplots(1, 2, figsize=(16, 9), gridspec_kw={"hspace": 0.20, "wspace": 0.20})
    axs = axs.flatten()
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    for s, t in zip(sims, titles):
        cutoff_density = int(s.split("_")[0])
        color = colors[cutoff_densities.index(cutoff_density)]
        scatter_point = "o"
        if "old" in s:
            scatter_point = "^"
        impact_point = secondary_impact_times[t]
        report_path = base_path + "{}/{}_reports/{}.csv".format(s, s, end_iteration)
        report = pd.read_csv(report_path)
        for index, p in enumerate(points):
            time = impact_point["time"] - impact_point['primary_impact_time']  # hrs after primary impact
            plot_x, plot_y = x[t][impact_point["iteration"]], float(str(report[p][0]).split(" ")[0])
            axs[index].scatter(
                plot_x, plot_y, color=color, marker=scatter_point, s=80
            )
            if "new" in s:
                axs[index].scatter(
                    [], [], s=80, marker="s", label=str(cutoff_density) + " $kg/m^3$"
                )
            axs[index].set_ylabel(rows_map[p][1:-1])
            axs[index].scatter(
                [], [], s=80, marker="^", label="GADGET M-ANEOS"
            )
            axs[index].scatter(
                [], [], s=80, marker="o", label="Stewart M-ANEOS"
            )
    for ax in axs:
        ax.grid(alpha=0.4)
        ax.legend()
    for ax in axs[-2:]:
        ax.set_xlabel(x_label)
    f = "angles"
    if r_dot_v:
        f = "r_dot_v"
    plt.savefig("{}_{}_secondary_impact_vs_disk_property.png".format(angle, f), format='png', dpi=200)

def plot_vs_disk_property_all():
    angles_path = "{}_secondary_impact_angles.csv".format(angle)
    r_dot_v_path = "{}_secondary_impact_angles_r_dot_v.csv".format(angle)
    r_dot_v_angle = "{}_secondary_impact_angles_angle_between_r_and_v.csv".format(angle)
    angles_df = pd.read_csv(angles_path, index_col="Unnamed: 0")
    r_dot_v_df = pd.read_csv(r_dot_v_path, index_col="Unnamed: 0")
    r_dot_v_angle_df = pd.read_csv(r_dot_v_angle, index_col="Unnamed: 0")
    # x_label_angle = "Impact Angle (deg.)"
    x_label_angle = r"$\theta$"
    # x_label_r_dot_v = "$r \cdot v$"
    x_label_r_dot_v = r"$\phi$"
    sims, titles = get_all_sims(high=False)
    points = ["DISK_MASS", "DISK_ANGULAR_MOMENTUM"]
    fig, axs = plt.subplots(2, 4, figsize=(27, 9), gridspec_kw={"hspace": 0.24, "wspace": 0.20})
    axs = axs.flatten()
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    for s, t in zip(sims, titles):
        cutoff_density = int(s.split("_")[0])
        color = colors[cutoff_densities.index(cutoff_density)]
        scatter_point = "o"
        if "old" in s:
            scatter_point = "^"
        if "high" in s:
            scatter_point = "*"
        impact_point = secondary_impact_times[t]
        impact_angle = angles_df[t][impact_point['iteration']]
        impact_r_dot_v = r_dot_v_df[t][impact_point['iteration']]
        impact_r_v_angle = 180 - r_dot_v_angle_df[t][impact_point['iteration']]
        report_path = base_path + "{}/{}_reports/{}.csv".format(s, s, end_iteration)
        report = pd.read_csv(report_path)
        time = impact_point["time"] - impact_point['primary_impact_time']  # hrs after primary impact
        # axs[0].scatter(
        #     impact_angle, time, color=color, marker=scatter_point, s=80
        # )
        # axs[1].scatter(
        #     impact_angle, float(str(report[points[0]][0]).split(" ")[0]), color=color, marker=scatter_point, s=80
        # )
        # axs[2].scatter(
        #     impact_angle, float(str(report[points[1]][0]).split(" ")[0]), color=color, marker=scatter_point, s=80
        # )
        # axs[3].scatter(
        #     impact_r_dot_v, time, color=color, marker=scatter_point, s=80
        # )
        # axs[4].scatter(
        #     impact_r_dot_v, float(str(report[points[0]][0]).split(" ")[0]), color=color, marker=scatter_point, s=80
        # )
        # axs[5].scatter(
        #     impact_r_dot_v, float(str(report[points[1]][0]).split(" ")[0]), color=color, marker=scatter_point, s=80
        # )
        # axs[6].scatter(
        #     time, float(str(report[points[0]][0]).split(" ")[0]), color=color, marker=scatter_point, s=80
        # )
        # axs[7].scatter(
        #     time, float(str(report[points[1]][0]).split(" ")[0]), color=color, marker=scatter_point, s=80
        # )
        axs[0].scatter(
            impact_angle, time, color=color, marker=scatter_point, s=80
        )
        axs[1].scatter(
            impact_angle, float(str(report[points[0]][0]).split(" ")[0]), color=color, marker=scatter_point, s=80
        )
        axs[2].scatter(
            impact_angle, float(str(report[points[1]][0]).split(" ")[0]), color=color, marker=scatter_point, s=80
        )
        axs[3].scatter(
            impact_r_v_angle, time, color=color, marker=scatter_point, s=80
        )
        axs[4].scatter(
            impact_r_v_angle, float(str(report[points[0]][0]).split(" ")[0]), color=color, marker=scatter_point, s=80
        )
        axs[5].scatter(
            impact_r_v_angle, float(str(report[points[1]][0]).split(" ")[0]), color=color, marker=scatter_point, s=80
        )
        axs[6].scatter(
            time, float(str(report[points[0]][0]).split(" ")[0]), color=color, marker=scatter_point, s=80
        )
        axs[7].scatter(
            time, float(str(report[points[1]][0]).split(" ")[0]), color=color, marker=scatter_point, s=80
        )
        
        for ax in axs[0:3]:
            ax.set_xlabel(x_label_angle)
        for ax in axs[3:6]:
            ax.set_xlabel(x_label_r_dot_v)
        for ax in [axs[6], axs[7]]:
            ax.set_xlabel("Time After Primary Impact (hrs)")

        for index in [1, 4, 6]:
            axs[index].set_ylabel(rows_map[points[0]][1:-1])
        for index in [2, 5, 7]:
            axs[index].set_ylabel(rows_map[points[1]][1:-1])
        for index in [0, 3]:
            axs[index].set_ylabel("Time After Primary Impact (hrs)")

    for c in cutoff_densities:
        axs[0].scatter(
            [], [], s=80, marker="s", color=colors[cutoff_densities.index(c)], label="{} kg/m$^3$".format(c)
        )
    axs[0].scatter(
        [], [], s=80, marker="^", color='black', label="GADGET M-ANEOS"
    )
    axs[0].scatter(
        [], [], s=80, marker="o", color='black', label="Stewart M-ANEOS"
    )
    if "high" in ",".join(str(i) for i in sims) and angle == "b073":
        axs[0].scatter(
            [], [], s=80, marker="*", color=colors[cutoff_densities.index(5)], label="5b073n-high"
        )

    axs[0].legend(fontsize=8)
    letters = list(string.ascii_lowercase)
    for index, ax in enumerate(axs):
        x1, x2, y1, y2 = ax.axis()
        x_loc = x1 + (0.02 * (x2 - x1))
        y_loc = y2 - (0.08 * (y2 - y1))
        ax.grid(alpha=0.4)
        ax.text(x_loc, y_loc, letters[index], fontweight="bold")
    plt.savefig("{}_secondary_impact_vs_disk_property.png".format(angle), format='png', dpi=200)

# plot_vs_disk_property(r_dot_v=False)
# plot_vs_disk_property(r_dot_v=True)
plot_vs_disk_property_all()
