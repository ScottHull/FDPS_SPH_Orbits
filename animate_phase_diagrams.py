#!/usr/bin/env python3
import os
import pandas as pd
import numpy as np
import string
import matplotlib.pyplot as plt

from src.animate import animate

plt.rcParams.update({'font.size': 14, })
# plt.style.use("dark_background")
plt.style.use('seaborn-colorblind')

runs = [('500_b073_new', '500b073S')]
base_path = "/home/theia/scotthull/Paper1_SPH/gi/"
new_phase_path = "src/phase_data/forstSTS__vapour_curve.txt"
old_phase_path = "src/phase_data/duniteN__vapour_curve.txt"
animate_phase_curves_path = os.path.join(base_path, "animate_phase_curves")

new_phase_df = pd.read_fwf(new_phase_path, skiprows=1,
                           names=["temperature", "density_sol_liq", "density_vap", "pressure",
                                  "entropy_sol_liq", "entropy_vap"])
old_phase_df = pd.read_fwf(old_phase_path, skiprows=1,
                           names=["temperature", "density_sol_liq", "density_vap", "pressure",
                                  "entropy_sol_liq", "entropy_vap"])

for p in [animate_phase_curves_path]:
    if not os.path.exists(p):
        os.mkdir(p)


fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

critical_point_new = max(new_phase_df['temperature'])
critical_point_old = max(old_phase_df['temperature'])



def shade_plot(s, ax):
    phase_df = new_phase_df
    cp = critical_point_new
    if "old" in s:
        phase_df = old_phase_df
        cp = critical_point_old
        # ax.text((3000, 5000), "100% Liquid", fontsize=8)
        # ax.text((3000, 5000), "100% Liquid", fontsize=8)
    ax.plot(
        phase_df['entropy_sol_liq'],
        phase_df['temperature'],
        linewidth=2.0,
        color='black'
    )
    ax.plot(
        phase_df['entropy_vap'],
        phase_df['temperature'],
        linewidth=2.0,
        color='black'
    )
    ax.fill_between(
        x=phase_df['entropy_vap'],
        y1=phase_df['temperature'],
        y2=cp,
        color=colors[-1],
        alpha=0.2,
        # label="100% Vapor"
    )
    ax.fill_between(
        x=phase_df['entropy_sol_liq'],
        y1=phase_df['temperature'],
        y2=cp,
        color=colors[-2],
        alpha=0.2,
        # label="100% Liquid"
    )
    ax.fill_between(
        x=sorted(list(phase_df['entropy_sol_liq']) + list(new_phase_df['entropy_vap'])),
        y1=cp,
        y2=1e10,
        color=colors[-3],
        alpha=0.2,
        # label="Supercritical"
    )
    ax.fill_between(
        x=phase_df['entropy_sol_liq'],
        y1=phase_df['temperature'],
        color=colors[-4],
        edgecolor="none",
        alpha=0.2,
        # label="Mixed"
    )
    ax.fill_between(
        x=phase_df['entropy_vap'],
        y1=phase_df['temperature'],
        color=colors[-4],
        edgecolor="none",
        alpha=0.2,
    )

    for phase, c in [("100% Vapor", colors[-1]), ("100% Liquid", colors[-2]), ("Supercritical", colors[-3]),
                     ("Liquid-Vapor\nMixture", colors[-4])]:
        ax.fill_between(
            x=[],
            y1=0,
            y2=0,
            color=c,
            alpha=0.2,
            label=phase
        )


def plot_phase_diagrams(ax_no_circ, ax_circ, s, t, iteration):
    marker = "."
    to_path = base_path + "{}/circularized_{}/{}.csv".format(s, s, iteration)
    df = pd.read_csv(to_path)
    disk = df[df['label'] == "DISK"]
    disk_filtered = disk[disk['circ_entropy_delta'] < 5000]
    temp, entropy, entropy_circ = disk_filtered['temperature'], disk_filtered['entropy'], disk_filtered['entropy'] + disk_filtered['circ_entropy_delta']
    cd = int(s.split("_")[0])
    color = colors[0]
    if "old" in s:
        phase_curve = old_phase_df
    to_index = 0
    if "b073" in s and "new" in s:
        to_index = 0
    elif "b073" in s and "old" in s:
        to_index = 1
    elif "b075" in s and "new" in s:
        to_index = 2
    else:
        to_index = 3
    label = r"$\rho_c = {}$ kg/m$^3$".format(cd)
    if "high" in s:
        label = t
    ax_no_circ.scatter(
        temp, entropy, s=1, marker=marker, alpha=0.6, color=color
    )
    ax_circ.scatter(
        temp, entropy_circ, s=1, marker=marker, alpha=0.6, color=color
    )

for run, t in runs:
    path = base_path + f"{run}/circularized_{run}"
    iterations = []
    for iteration in os.listdir(path):
        iteration = int(iteration.split(".")[0])
        iterations.append(iteration)
        fig, axs = plt.subplots(1, 2, figsize=(16, 9))
        axs = axs.flatten()
        for ax in axs:
            ax.grid(alpha=0.4)
            ax.set_xlim(1000, 12000)
            ax.set_ylim(0, 12500)
            shade_plot(s="new", ax=ax)
        axs[0].set_label("No Circularization")
        axs[1].set_label("Circularization")

        plot_phase_diagrams(axs[0], axs[1], run, t, iteration)

        plt.tight_layout()

        fig.subplots_adjust(right=0.80)

        plt.savefig(animate_phase_curves_path + f"/{iteration}.png", format='png', dpi=200)

    animate(
        start_time=min(iterations),
        end_time=max(iterations),
        interval=int(iterations[1] - iterations[0]),
        path=animate_phase_curves_path,
        fps=30,
        filename="{}_animate_phase_curves.mp4".format(t),
    )
