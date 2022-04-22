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

plt.rcParams.update({'font.size': 8, })
# plt.style.use("dark_background")
plt.style.use('seaborn-colorblind')

min_iteration = 50
iteration = 190
angle = 'b073'
cutoff_densities = [5, 500, 1000, 2000]
number_processes = 200
vsquare_scale = 5e7
hsquare_scale = 4e7
base_path = "/home/theia/scotthull/Paper1_SPH/gi/"


secondary_impact_lims = {
    '5b073n': {
        'si_min_x': 0.4e7,
        'si_max_x': -0.5e7,
        'si_min_y': -2.3e7,
        'si_max_y': -1.6e7,
        'tail_min_x': -2.7e7,
        'tail_max_x': -1.1e7,
        'tail_min_y': -4.1e7,
        'tail_max_y': -2.3e7,
    },
    '500b073n': {
        'si_min_x': -1.1e7,
        'si_max_x': -1.85e7,
        'si_min_y': -2.05e7,
        'si_max_y': -2.7e7,
        'tail_min_x': -1.3e7,
        'tail_max_x': -2.2e7,
        'tail_min_y': -2.7e7,
        'tail_max_y': -4.2e7,
    },
    '1000b073n': {
        'si_min_x': -0.95e7,
        'si_max_x': -1.65e7,
        'si_min_y': -1.8e7,
        'si_max_y': -2.5e7,
        'tail_min_x': -1.2e7,
        'tail_max_x': -2.3e7,
        'tail_min_y': -2.5e7,
        'tail_max_y': -4.2e7,
    },
    '2000b073n': {
        'si_min_x': -0.45e7,
        'si_max_x': -1e7,
        'si_min_y': -1.3e7,
        'si_max_y': -2e7,
        'tail_min_x': -2.6e7,
        'tail_max_x': -3e7,
        'tail_min_y': -2e7,
        'tail_max_y': -4e7,
    },
    '5b0b73o': {
        'si_min_x': -0.8e7,
        'si_max_x': -1.55e7,
        'si_min_y': -1.85e7,
        'si_max_y': -2.5e7,
        'tail_min_x': -1e7,
        'tail_max_x': -2.4e7,
        'tail_min_y': -2.5e7,
        'tail_max_y': -4.1e7,
    },
    '500b0b73o': {
        'si_min_x': -1.1e7,
        'si_max_x': -1.8e7,
        'si_min_y': -2.1e7,
        'si_max_y': -3e7,
        'tail_min_x': -1.1e7,
        'tail_max_x': -3e7,
        'tail_min_y': -3e7,
        'tail_max_y': -4.5e7,
    },
    '1000b0b73o': {
        'si_min_x': -1.05e7,
        'si_max_x': -1.9e7,
        'si_min_y': -2.2e7,
        'si_max_y': -3.15e7,
        'tail_min_x': -1.5e7,
        'tail_max_x': -2e7,
        'tail_min_y': -2.5e7,
        'tail_max_y': -4.1e7,
    },
    '2000b0b73o': {
        'si_min_x': -1.05e7,
        'si_max_x': -1.6e7,
        'si_min_y': -2e7,
        'si_max_y': -2.6e7,
        'tail_min_x': -1.1e7,
        'tail_max_x': -1.05e7,
        'tail_min_y': -2.6e7,
        'tail_max_y': -4.1e7,
    },
}


def get_all_sims(angle, high=True):
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
    if high:
        output_name = fformat.format(5, angle, "new") + "_high"
        names.append(output_name)
        title_name = tformat.format(5, angle, "n") + "-high"
        titles.append(title_name)
    return names, titles

def __build_scene(iteration, times, to_save_path, dfs, sims, titles, target_coms, impactor_coms):
    num_new = len([i for i in sims if "new" in i])
    num_old = len([i for i in sims if "old" in i])
    num_rows = max([num_new, num_old])
    plt.style.use("dark_background")
    fig, axs = plt.subplots(num_rows, 2, figsize=(16, 32), sharex='all',
                            gridspec_kw={"hspace": 0.10, "wspace": 0.12})
    fig.patch.set_facecolor('xkcd:black')
    axs = axs.flatten()
    for ax in axs:
        ax.set_xlim(-hsquare_scale, 0.5 * hsquare_scale)
        ax.set_ylim(-vsquare_scale, 2000)
        ax.grid(alpha=0.4)
    index_new, index_old = 0, 1
    for s, t in zip(sims, titles):
        df, target_com, impactor_com, time = dfs[t][-1], target_coms[t][-1], impactor_coms[t][-1], times[t][-1]
        # df = df[df['z'] < 0]
        target_silicate = df[df['tag'] == 0]
        target_iron = df[df['tag'] == 1]
        impactor_silicate = df[df['tag'] == 2]
        impactor_iron = df[df['tag'] == 3]
        t1 = ["Target Silicate", "Target Iron", "Impactor Silicate", "Impactor Iron"]
        t2 = [target_silicate, target_iron, impactor_silicate, impactor_iron]

        to_index = index_new
        if "old" in s:
            to_index = index_old
        for a, b, in zip(t1, t2):
            axs[to_index].scatter(
                b['x'], b['y'], s=1, marker='.', label=a
            )

        axs[to_index].scatter(
            impactor_com[0], impactor_com[1], s=60, c='pink', marker="*", label="Impactor COM"
        )
        axs[to_index].scatter(
            target_com[0], target_com[1], s=60, c='green', marker="*", label="Target COM"
        )
        axs[to_index].plot(
            [target_com[0], impactor_com[0]], [target_com[1], impactor_com[1]], linewidth=2.0, color="aqua"
        )
        axs[to_index].plot(
            [target_com[0], target_com[0]], [target_com[1], impactor_com[1]], linewidth=2.0, color="aqua"
        )
        axs[to_index].plot(
            [impactor_com[0], target_com[0]], [impactor_com[1], impactor_com[1]], linewidth=2.0, color="aqua"
        )

        axs[to_index].set_title("{} {} hrs. ({})".format(t, time, iteration))
        if "old" in s:
            index_old += 2
        else:
            index_new += 2
    axs[0].legend(loc='upper left')
    plt.savefig(to_save_path + "/{}.png".format(iteration), format='png', dpi=200)


def __plot_secondary(iteration, times, to_save_path, dfs, sims, titles, target_coms, impactor_coms, sis, tails, not_classifieds):
    num_new = len([i for i in sims if "new" in i])
    num_old = len([i for i in sims if "old" in i])
    num_rows = max([num_new, num_old])
    plt.style.use("dark_background")
    fig, axs = plt.subplots(num_rows, 2, figsize=(16, 32), sharex='all',
                            gridspec_kw={"hspace": 0.10, "wspace": 0.12})
    fig.patch.set_facecolor('xkcd:black')
    axs = axs.flatten()
    for ax in axs:
        ax.set_xlim(-hsquare_scale, 0.5 * hsquare_scale)
        ax.set_ylim(-vsquare_scale, 2000)
        ax.grid(alpha=0.4)
    index_new, index_old = 0, 1
    for s, t in zip(sims, titles):
        df, target_com, impactor_com, time, si, tail, not_classified = dfs[t][-1], target_coms[t][-1], \
                                                                       impactor_coms[t][-1], times[t][-1], sis[t][-1], \
                                                                       tails[t][-1], not_classifieds[t][-1]
        # df = df[df['z'] < 0]
        # target_silicate = df[df['tag'] == 0]
        # target_iron = df[df['tag'] == 1]
        # impactor_silicate = df[df['tag'] == 2]
        # impactor_iron = df[df['tag'] == 3]
        t1 = ["Secondary Impactor", "Debris Tail", "Other"]
        t2 = [si, tail, not_classified]

        to_index = index_new
        if "old" in s:
            to_index = index_old
        for a, b, in zip(t1, t2):
            axs[to_index].scatter(
                b['x'], b['y'], s=1, marker='.', label=a
            )

        axs[to_index].scatter(
            impactor_com[0], impactor_com[1], s=60, c='pink', marker="*", label="Impactor COM"
        )
        axs[to_index].scatter(
            target_com[0], target_com[1], s=60, c='green', marker="*", label="Target COM"
        )
        axs[to_index].plot(
            [target_com[0], impactor_com[0]], [target_com[1], impactor_com[1]], linewidth=2.0, color="aqua"
        )
        axs[to_index].plot(
            [target_com[0], target_com[0]], [target_com[1], impactor_com[1]], linewidth=2.0, color="aqua"
        )
        axs[to_index].plot(
            [impactor_com[0], target_com[0]], [impactor_com[1], impactor_com[1]], linewidth=2.0, color="aqua"
        )

        axs[to_index].set_title("{} {} hrs. ({})".format(t, time, iteration))
        if "old" in s:
            index_old += 2
        else:
            index_new += 2
    axs[0].legend(loc='upper left')
    plt.savefig(to_save_path + "/{}.png".format(iteration), format='png', dpi=200)


def get_com(x, y, z, mass):
    total_mass = float(np.sum(mass))
    x_center = sum(x * mass) / total_mass
    y_center = sum(y * mass) / total_mass
    z_center = sum(z * mass) / total_mass
    return np.array([x_center, y_center, z_center])


def get_secondary_imp_and_tail_particles(title, df):
    lims = secondary_impact_lims[title]
    x_lims_si = [lims['si_min_x'], lims['si_max_x']]
    y_lims_si = [lims['si_min_y'], lims['si_max_y']]
    x_lims_tail = [lims['tail_min_x'], lims['tail_max_x']]
    y_lims_tail = [lims['tail_min_y'], lims['tail_max_y']]
    min_x_si, max_x_si = min(x_lims_si), max(x_lims_si)
    min_y_si, max_y_si = min(y_lims_si), max(y_lims_si)
    min_x_tail, max_x_tail = min(x_lims_tail), max(x_lims_tail)
    min_y_tail, max_y_tail = min(y_lims_tail), max(y_lims_tail)
    
    si = df[min_x_si <= df['x'] <= max_x_si]
    si = si[si['x'] >= min_x_si]
    si = si[si['y'] <= max_y_si]
    si = si[si['y'] >= min_y_si]
    tail = df[~df['id'].isin(si['id'].tolist())]
    tail = df[df['x'] <= max_x_tail]
    tail = tail[tail['x'] >= min_x_tail]
    tail = df[df['y'] <= max_y_tail]
    tail = tail[tail['y'] >= min_y_tail]
    not_classified = df[~df['id'].isin(si['id'].tolist())]
    not_classified = not_classified[~not_classified['id'].isin(tail['id'].tolist())]
    return si, tail, not_classified


def get_impactor_com_particles(output_name):
    output_path = base_path + output_name + "/circularized_{}".format(output_name)
    df = pd.read_csv(output_path + "/{}.csv".format(min_iteration))
    impactor_iron = df[df['tag'] == 3]
    impactor_iron = impactor_iron[impactor_iron['radius'] > 1.5e7]
    return impactor_iron['id']


def get_secondary_and_tail():
    sims, titles = get_all_sims(angle=angle, high=False)
    impact_angles = {}
    target_coms = {}
    impactor_coms = {}
    imp_ids = {}
    times = {}
    r_dot_vs = {}
    r_dot_v_angles = {}
    r_vectors = {}
    v_vectors = {}
    sis = {}
    tails = {}
    not_classifieds = {}
    for output_name, title in zip(sims, titles):
        imp_ids.update({title: get_impactor_com_particles(output_name)})
    dfs = {}
    for output_name, title in zip(sims, titles):
        if title not in impact_angles.keys():
            impact_angles.update({title: []})
            impactor_coms.update({title: []})
            target_coms.update({title: []})
            times.update({title: []})
            r_dot_vs.update({title: []})
            r_dot_v_angles.update({title: []})
            r_vectors.update({title: []})
            v_vectors.update({title: []})
            sis.update({title: []})
            tails.update({title: []})
            not_classifieds.update({title: []})
        if title not in dfs.keys():
            dfs.update({title: []})
        output_path = base_path + output_name + "/circularized_{}".format(output_name)
        # df = pd.read_csv(output_path + "/{}.csv".format(iteration))
        path = base_path + "{}/{}".format(output_name, output_name)
        to_fname = "merged_{}_{}.dat".format(iteration, randint(0, 100000))
        cf = CombineFile(num_processes=number_processes, time=iteration, output_path=path, to_fname=to_fname)
        combined_file = cf.combine()
        formatted_time = round(cf.sim_time * 0.000277778, 2)
        times[title].append(formatted_time)
        f = os.getcwd() + "/{}".format(to_fname)
        headers = ["id", "tag", "mass", "x", "y", "z", "vx", "vy", "vz", "density", "internal energy", "pressure",
                   "potential energy", "entropy", "temperature"]
        df = pd.read_csv(f, skiprows=2, header=None, delimiter="\t", names=headers)
        os.remove(f)
        zero_coords = get_com(df[df['tag'] == 1]['x'], df[df['tag'] == 1]['y'],
                              df[df['tag'] == 1]['z'], df[df['tag'] == 1]['mass'])
        df['x'] -= zero_coords[0]
        df['y'] -= zero_coords[1]
        df['z'] -= zero_coords[2]
        # df['radius'] = [(i ** 2 + j ** 2 + k ** 2) ** (1 / 2) for i, j, k in zip(df['x'], df['y'], df['z'])]

        target_iron = df[df['tag'] == 1]
        impactor_iron = df[df['tag'] == 3]

        # impactor_iron = impactor_iron[impactor_iron['id'] in imp_ids[title]]
        impactor_iron = impactor_iron[impactor_iron['id'].isin(imp_ids[title].tolist())]

        target_com = get_com(target_iron['x'], target_iron['y'], target_iron['z'], target_iron['mass'])
        impactor_com = get_com(impactor_iron['x'], impactor_iron['y'], impactor_iron['z'], impactor_iron['mass'])
        si, tail, not_classified = get_secondary_imp_and_tail_particles(title, df)


        dfs[title].append(df)
        target_coms[title].append(target_com)
        impactor_coms[title].append(impactor_com)
        sis[title].append(si)
        tails[title].append(tail)
        not_classifieds[title].append(not_classified)

    to_save_path = "{}_secondary_impact_structures".format(angle)
    if not os.path.exists(to_save_path):
        os.mkdir(to_save_path)
    # __build_scene(iteration=iteration, dfs=dfs, sims=sims, titles=titles,
    #               target_coms=target_coms, impactor_coms=impactor_coms, to_save_path=to_save_path, times=times)
    __plot_secondary(iteration=iteration, dfs=dfs, sims=sims, titles=titles,
                  target_coms=target_coms, impactor_coms=impactor_coms, to_save_path=to_save_path, times=times,
                     sis=sis, tails=tails, not_classifieds=not_classifieds)

get_secondary_and_tail()
