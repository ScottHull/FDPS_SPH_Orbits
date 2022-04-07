#!/usr/bin/env python3
import os
import csv
import pandas as pd
import numpy as np
from math import sin, cos, tan, pi
from operator import contains
from random import randint
import multiprocessing as mp
import matplotlib.pyplot as plt

plt.style.use("dark_background")

base_path = "/home/theia/scotthull/Paper1_SPH/gi/"
angle = "b073"
cutoff_densities = [5, 500, 1000, 2000]
min_iteration = 50
max_iteration = 600
increment = 50
square_scale = 6e7


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


def get_com(x, y, z, mass):
    total_mass = float(np.sum(mass))
    x_center = sum(x * mass) / total_mass
    y_center = sum(y * mass) / total_mass
    z_center = sum(z * mass) / total_mass
    return x_center, y_center, z_center


def get_impactor_com_particles(output_name):
    output_path = base_path + output_name + "/circularized_{}".format(output_name)
    df = pd.read_csv(output_path + "/{}.csv".format(min_iteration))
    impactor_iron = df[df['tag'] == 3]
    impactor_iron = impactor_iron[impactor_iron['radius'] > 1e7]
    return impactor_iron['id']


def get_angle(target_com, impactor_com):
    x_offset = impactor_com[0] - target_com[0]
    y_offset = impactor_com[1] - target_com[1]
    return tan(x_offset / y_offset) * (180 / pi)


def __build_scene(iteration, dfs, sims, titles, impact_angles, target_coms, impactor_coms, to_save_path):
    num_new = len([i for i in sims if "new" in i])
    num_old = len([i for i in sims if "old" in i])
    num_rows = max([num_new, num_old])
    plt.style.use("dark_background")
    fig, axs = plt.subplots(num_rows, 2, figsize=(16, 32), sharex='all',
                            gridspec_kw={"hspace": 0.10, "wspace": 0.12})
    fig.patch.set_facecolor('xkcd:black')
    axs = axs.flatten()
    for ax in axs:
        ax.set_xlim(-square_scale, square_scale)
        ax.set_ylim(-square_scale, square_scale)
        ax.set_xticks([], minor=False)
        ax.set_yticks([], minor=False)
    index_new, index_old = 0, 1
    for s, t in zip(sims, titles):
        df, impact_angle, target_com, impactor_com = dfs[t][-1], impact_angles[t][-1], \
                                                     target_coms[t][-1], impactor_coms[t][-1]
        target_material = df[df['tag'] <= 1]
        impactor_material = df[df['tag'] > 1]
        to_index = index_new
        if "old" in s:
            to_index = index_old
        axs[to_index].scatter(
            target_material['x'], target_material['y'], s=1, label="Target Material"
        )
        axs[to_index].scatter(
            impactor_material['x'], impactor_material['y'], s=1, label="Impactor Material"
        )
        axs[to_index].scatter(
            impactor_com[0], impactor_com[1], s=60, c='pink', marker="*", label="Impactor COM"
        )
        axs[to_index].scatter(
            target_com[0], target_com[1], s=60, c='green', marker="*", label="Target COM"
        )
        axs[to_index].plot(
            [target_com[0], impactor_com[0]], [target_com[1], impactor_com[1]], linewidth=2.0, color="white"
        )
        axs[to_index].plot(
            [target_com[0], target_com[0]], [target_com[1], impactor_com[1]], linewidth=2.0, color="white"
        )
        axs[to_index].plot(
            [impactor_com[0], target_com[0]], [impactor_com[1], impactor_com[1]], linewidth=2.0, color="white"
        )
        axs[to_index].set_title("{} - {} deg.".format(iteration, round(impact_angle, 2)))
        if "old" in s:
            index_old += 2
        else:
            index_new += 2
    axs[0].legend(loc='upper left')
    plt.savefig(to_save_path + "/{}.png".format(iteration), format='png', dpi=200)


def get_impact_angles():
    sims, titles = get_all_sims(high=False)
    impact_angles = {}
    target_coms = {}
    impactor_coms = {}
    imp_ids = {}
    for output_name, title in zip(sims, titles):
        imp_ids.update({title: get_impactor_com_particles(output_name)})
    for iteration in np.arange(min_iteration, max_iteration + increment, increment):
        dfs = {}
        for output_name, title in zip(sims, titles):
            if title not in impact_angles.keys():
                impact_angles.update({title: []})
                impactor_coms.update({title: []})
                target_coms.update({title: []})
            if title not in dfs.keys():
                dfs.update({title: []})
            output_path = base_path + output_name + "/circularized_{}".format(output_name)
            df = pd.read_csv(output_path + "/{}.csv".format(iteration))
            target_iron = df[df['tag'] == 1]
            impactor_iron = df[df['tag'] == 3]
            # impactor_iron = impactor_iron[impactor_iron['id'] in imp_ids[title]]
            impactor_iron = impactor_iron[impactor_iron['id'].isin(imp_ids[title].tolist())]
            target_com = get_com(target_iron['x'], target_iron['y'], target_iron['z'], target_iron['mass'])
            impactor_com = get_com(impactor_iron['x'], impactor_iron['y'], impactor_iron['z'], impactor_iron['mass'])
            imp_angle = get_angle(target_com, impactor_com)
            dfs[title].append(df)
            impact_angles[title].append(imp_angle)
            target_coms[title].append(target_com)
            impactor_coms[title].append(impactor_com)
        to_save_path = "{}_secondary_impact_angles_scences".format(angle)
        if not os.path.exists(to_save_path):
            os.mkdir(to_save_path)
        __build_scene(iteration=iteration, dfs=dfs, sims=sims, titles=titles, impact_angles=impact_angles,
                      target_coms=target_coms, impactor_coms=impactor_coms, to_save_path=to_save_path)

    df = pd.DataFrame(impact_angles, index=np.arange(min_iteration, max_iteration + increment, increment))
    df.to_csv("{}_secondary_impact_angles.csv".format(angle))

get_impact_angles()
