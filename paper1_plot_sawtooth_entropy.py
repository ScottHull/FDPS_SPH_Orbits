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

from src.animate import animate
from src.combine import CombineFile
from src.identify import ParticleMap
from src.report import get_sim_report, write_report_at_time, build_latex_table_from_disk_report, rows_map


plt.rcParams.update({'font.size': 16, })
# plt.style.use("dark_background")
plt.style.use('seaborn-colorblind')

base_path = "/home/theia/scotthull/Paper1_SPH/gi/"
angle = "b073"
cutoff_densities = [5]
min_density = 0
max_density = 20
min_iteration = 0
max_iteration = 1800
endstate_iteration = max_iteration
increment = 1
number_processes = 200

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


def get_endstate(s):
    path = base_path + "{}/circularized_{}".format(s, s)
    df = pd.read_csv(path + "/{}.csv".format(endstate_iteration))
    return df[df['label'] == "DISK"]


# make a 3 column plot of density, entropy, internal energy
ylabels = [r'Density (kg/m$^3$)', r'Entropy (J/kg/K)', r'Internal Energy (kJ)']

time = {}
entropy = {}
internal_energy = {}
density = {}
temperature = {}
names, titles = get_all_sims()
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
for s, t in zip(names, titles):
    if os.path.exists(f"paper1_sawtooth_{s}"):
        shutil.rmtree(f"paper1_sawtooth_{s}")
    os.mkdir(f"paper1_sawtooth_{s}")
    # get the endstate df
    endstate = get_endstate(s)
    endstate_whole_disk = endstate['id'].tolist()
    endstate_target_particles = endstate[endstate['entropy'] > 9000]
    endstate_target_particles = endstate_target_particles[endstate_target_particles['density'] == cutoff_densities[0]]
    # get 5 random particle ids from the endstate df
    endstate = endstate_target_particles.sample(n=5)['id'].tolist()
    # get the color cycle
    time.update({s: {i: [] for i in endstate}})
    entropy.update({s: {i: [] for i in endstate}})
    internal_energy.update({s: {i: [] for i in endstate}})
    density.update({s: {i: [] for i in endstate}})
    temperature.update({s: {i: [] for i in endstate}})
    prev_entropy = {}
    for iteration in np.arange(min_iteration, max_iteration + increment, increment):
        path = base_path + "{}/{}".format(s, s)
        to_fname = "merged_{}_{}.dat".format(iteration, randint(0, 100000))
        cf = CombineFile(num_processes=number_processes, time=iteration, output_path=path, to_fname=to_fname)
        combined_file = cf.combine()
        formatted_time = round(cf.sim_time * 0.000277778, 2)
        f = os.getcwd() + "/{}".format(to_fname)
        headers = ["id", "tag", "mass", "x", "y", "z", "vx", "vy", "vz", "density", "internal energy", "pressure",
                   "potential energy", "entropy", "temperature"]
        df = pd.read_csv(f, skiprows=2, header=None, delimiter="\t", names=headers)
        disk = df[df['id'].isin(endstate)]
        whole_disk = df[df['id'].isin(endstate_whole_disk)]
        os.remove(f)
        for i in endstate:
            time[s][i].append(formatted_time)
            entropy[s][i].append(disk[disk['id'] == i]['entropy'].values[0])
            internal_energy[s][i].append(disk[disk['id'] == i]['internal energy'].values[0])
            density[s][i].append(disk[disk['id'] == i]['density'].values[0])
            temperature[s][i].append(disk[disk['id'] == i]['temperature'].values[0])

    #     # make a 3 column plot of density, entropy, internal energy
    #     fig2, axs2 = plt.subplots(1, 4, figsize=(24, 6), sharex='all')
    #     axs2 = axs2.flatten()
    #     for ax in axs2:
    #         ax.grid()
    #         ax.set_xlabel(r"Density (kg/m$^3$", fontsize=16)
    #         ax.set_title(f"{formatted_time} hrs.")
    #         ax.set_xlim(min_density, max_density)
    #     axs2[0].set_ylabel(r'Entropy (J/kg/K)', fontsize=16)
    #     axs2[1].set_ylabel(r'Internal Energy (kJ)', fontsize=16)
    #     axs2[2].set_ylabel(r'Temperature (K)', fontsize=16)
    #     axs2[3].set_ylabel(r'Delta Entropy (J/kg/K)', fontsize=16)
    #     axs2[0].set_ylim(0, 15000)
    #     axs2[1].set_ylim(0, 100000)
    #     axs2[2].set_ylim(0, 15000)
    #     axs2[3].set_ylim(0, 200)
    #     # axs2[1].set_ylim(0, 5e7)
    #     axs2[0].scatter(
    #         whole_disk['density'], whole_disk['entropy'], marker=".", s=5
    #     )
    #     axs2[0].scatter(
    #         disk['density'], disk['entropy'], marker=".", s=80, c='r'
    #     )
    #     axs2[1].scatter(
    #         whole_disk['density'], whole_disk['internal energy'] / 1000, marker=".", s=5
    #     )
    #     axs2[1].scatter(
    #         disk['density'], disk['internal energy'] / 1000, marker=".", s=80, c='r'
    #     )
    #     # annotate the average y value in the upper right corner
    #     axs2[0].annotate(
    #         "Average Entropy:\n{:.2f}".format(np.mean(whole_disk['entropy'])),
    #         xy=(0.65, 0.9),
    #         xycoords='axes fraction',
    #         fontsize=10,
    #     )
    #     axs2[1].annotate(
    #         "Average\nInternal Energy: \n{:.2f}".format(np.mean(whole_disk['internal energy']) / 1000),
    #         xy=(0.65, 0.9),
    #         xycoords='axes fraction',
    #         fontsize=10,
    #     )
    #     axs2[2].scatter(
    #         whole_disk['density'], whole_disk['temperature'], marker=".", s=5
    #     )
    #     axs2[2].scatter(
    #         disk['density'], disk['temperature'], marker=".", s=80, c='r'
    #     )
    #     axs2[2].annotate(
    #         "Average\nTemperature:\n{:.2f}".format(np.mean(whole_disk['temperature'])),
    #         xy=(0.50, 0.9),
    #         xycoords='axes fraction',
    #         fontsize=10,
    #     )
    #     try:
    #         if iteration > min_iteration:
    #             delta_s_wd = {i: prev_entropy[i] - s_curr for i, s_curr in zip(whole_disk['id'], whole_disk['entropy'])}
    #             delta_s_wd_at_rhoc = {i: prev_entropy[i] - s_curr for i, s_curr, density in zip(whole_disk['id'], whole_disk['entropy'], whole_disk['density'])}
    #             delta_s_d = {i: prev_entropy[i] - s_curr for i, s_curr in zip(disk['id'], disk['entropy'])}
    #
    #             axs2[3].scatter(
    #                 whole_disk['density'], delta_s_wd.values(), marker=".", s=5
    #             )
    #             axs2[3].scatter(
    #                 disk['density'], delta_s_d.values(), marker=".", s=80, c='r'
    #             )
    #             axs2[3].annotate(
    #                 "Average Delta Entropy at\n" + r"$\rho_c$:" + "{:.2f}".format(np.mean(list(delta_s_wd_at_rhoc.values()))),
    #                 xy=(0.50, 0.9),
    #                 xycoords='axes fraction',
    #                 fontsize=10,
    #             )
    #     except:
    #         pass
    #     # axs2[-1].scatter(
    #     #     [], [], label="Stewart M-ANEOS"
    #     # )
    #     # axs2[-1].scatter(
    #     #     [], [], label="N-SPH M-ANEOS"
    #     # )
    #     # axs2[-1].legend(fontsize=14, loc='lower right')
    #     plt.tight_layout()
    #     plt.savefig("paper1_sawtooth_{}/{}.png".format(s, iteration), dpi=200)
    #
    #     prev_entropy = dict(zip(whole_disk['id'], whole_disk['entropy']))
    #
    # animate(
    #     start_time=min_iteration,
    #     end_time=max_iteration,
    #     interval=increment,
    #     path="paper1_sawtooth_{}".format(s),
    #     fps=80,
    #     filename="{}_{}_sawtooth.mp4".format(s, cutoff_densities[0]),
    # )
    # break


fig, axs = plt.subplots(1, 4, figsize=(24, 6), sharex='all')
axs = axs.flatten()
for ax in axs:
    ax.grid()
    ax.set_xlabel("Time (hours)", fontsize=16)
ylabels = [r'Density (kg/m$^3$)', r'Entropy (J/kg/K)', r'Internal Energy (kJ)', r'Temperature (K)']
for ax, y in zip(axs, ylabels):
    ax.set_ylabel(y, fontsize=16)

for s in density.keys():
    linestyle = "-"
    if "old" in s:
        linestyle = "--"
    for ax, i in zip(axs, [density[s], entropy[s], internal_energy[s], temperature[s]]):
        for index, j in enumerate(density[s].keys()):
            ax.plot(time[s][j], i[j], linestyle=linestyle, color=colors[index])

axs[0].set_ylim(bottom=min_density, top=max_density)

# axs[-1].plot(
#     [], [], linestyle="-", label="Stewart M-ANEOS", color="black"
# )
# axs[-1].plot(
#     [], [], linestyle="--", label="N-SPH M-ANEOS", color="black"
# )
# axs[-1].legend(fontsize=14, loc='lower right')
plt.tight_layout()
plt.savefig(f"{cutoff_densities[0]}_{angle}_sawtooth_entropy.png", dpi=300)



