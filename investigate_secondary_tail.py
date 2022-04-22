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
max_iteration = 1800
angle = 'b073'
cutoff_densities = [5, 500, 1000, 2000]
number_processes = 200
vsquare_scale = 5e7
hsquare_scale = 4e7
base_path = "/home/theia/scotthull/Paper1_SPH/gi/"

secondary_impact_lims = {
    '5b073n': {
        'si_min_x': -0.7e7,
        'si_max_x': -1.4e7,
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
        'tail_min_x': -0.6e7,
        'tail_max_x': -3e7,
        'tail_min_y': -2e7,
        'tail_max_y': -4e7,
    },
    '5b073o': {
        'si_min_x': -0.8e7,
        'si_max_x': -1.55e7,
        'si_min_y': -1.85e7,
        'si_max_y': -2.5e7,
        'tail_min_x': -1e7,
        'tail_max_x': -2.4e7,
        'tail_min_y': -2.5e7,
        'tail_max_y': -4.1e7,
    },
    '500b073o': {
        'si_min_x': -1.1e7,
        'si_max_x': -1.8e7,
        'si_min_y': -2.1e7,
        'si_max_y': -3e7,
        'tail_min_x': -1.1e7,
        'tail_max_x': -3e7,
        'tail_min_y': -3e7,
        'tail_max_y': -4.5e7,
    },
    '1000b073o': {
        'si_min_x': -1.05e7,
        'si_max_x': -1.9e7,
        'si_min_y': -2.2e7,
        'si_max_y': -3.15e7,
        'tail_min_x': -1.1e7,
        'tail_max_x': -2e7,
        'tail_min_y': -2.5e7,
        'tail_max_y': -4.1e7,
    },
    '2000b073o': {
        'si_min_x': -1.05e7,
        'si_max_x': -1.6e7,
        'si_min_y': -2e7,
        'si_max_y': -2.6e7,
        'tail_min_x': -1.1e7,
        'tail_max_x': -2.5e7,
        'tail_min_y': -2.55e7,
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


def __plot_secondary(iteration, times, to_save_path, dfs, sims, titles, target_coms, impactor_coms, sis, tails,
                     not_classifieds):
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

    si = df[df['x'] <= max_x_si]
    si = si[si['x'] >= min_x_si]
    si = si[si['y'] <= max_y_si]
    si = si[si['y'] >= min_y_si]
    tail = df[~df['id'].isin(si['id'].tolist())]
    tail = tail[tail['x'] <= max_x_tail]
    tail = tail[tail['x'] >= min_x_tail]
    tail = tail[tail['y'] <= max_y_tail]
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
    si_and_tail_data = {}
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
            si_and_tail_data.update({title: []})
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

        end_state_file = base_path + "{}/circularized_{}/{}.csv".format(output_name, output_name, max_iteration)
        end_state_df = pd.read_csv(end_state_file, index_col="id")

        si_mass, tail_mass = si['mass'], tail['mass']
        si_position_vector = list(zip(si['x'], si['y'], si['z']))
        si_velocity_vector = list(zip(si['vx'], si['vy'], si['vz']))
        tail_position_vector = list(zip(tail['x'], tail['y'], tail['z']))
        tail_velocity_vector = list(zip(tail['vx'], tail['vy'], tail['vz']))
        si_ang_mom = sum([np.linalg.norm(i * np.cross(j, k)) for i, j, k in
                          list(zip(si_mass, si_position_vector, si_velocity_vector))])
        tail_ang_mom = sum([np.linalg.norm(i * np.cross(j, k)) for i, j, k in
                            list(zip(tail_mass, tail_position_vector, tail_velocity_vector))])
        si_pct_tar_silicate = len(si[si['tag'] == 0]) / len(si)
        si_pct_tar_iron = len(si[si['tag'] == 1]) / len(si)
        si_pct_imp_silicate = len(si[si['tag'] == 2]) / len(si)
        si_pct_imp_iron = len(si[si['tag'] == 3]) / len(si)
        tail_pct_tar_silicate = len(tail[tail['tag'] == 0]) / len(tail)
        tail_pct_tar_iron = len(tail[tail['tag'] == 1]) / len(tail)
        tail_pct_imp_silicate = len(tail[tail['tag'] == 2]) / len(tail)
        tail_pct_imp_iron = len(tail[tail['tag'] == 3]) / len(tail)
        si_pct_silicate = len(si[si['tag'] % 2 == 0]) / len(si)
        si_pct_iron = len(si[si['tag'] % 2 != 0]) / len(si)
        tail_pct_silicate = len(tail[tail['tag'] % 2 == 0]) / len(tail)
        tail_pct_iron = len(tail[tail['tag'] % 2 != 0]) / len(tail)
        si_labels = [end_state_df['label'][i] for i in si['id']]
        tail_labels = [end_state_df['label'][i] for i in tail['id']]
        si_planet = len([i for i in si_labels if i == "PLANET"]) / len(si)
        si_disk = len([i for i in si_labels if i == "DISK"]) / len(si)
        si_escape = len([i for i in si_labels if i == "ESCAPE"]) / len(si)
        tail_planet = len([i for i in tail_labels if i == "PLANET"]) / len(tail)
        tail_disk = len([i for i in tail_labels if i == "DISK"]) / len(tail)
        tail_escape = len([i for i in tail_labels if i == "ESCAPE"]) / len(tail)
        planet_from_si = len(
            [i for i in end_state_df.index if i in si['id'] and end_state_df['label'][i] == "PLANET"]) / len(
            end_state_df[end_state_df['label'] == "PLANET"])
        disk_from_si = len(
            [i for i in end_state_df.index if i in si['id'] and end_state_df['label'][i] == "DISK"]) / len(
            end_state_df[end_state_df['label'] == "DISK"])
        escape_from_si = len(
            [i for i in end_state_df.index if i in si['id'] and end_state_df['label'][i] == "ESCAPE"]) / len(
            end_state_df[end_state_df['label'] == "ESCAPE"])
        planet_from_tail = len(
            [i for i in end_state_df.index if i in tail['id'] and end_state_df['label'][i] == "PLANET"]) / len(
            end_state_df[end_state_df['label'] == "PLANET"])
        disk_from_tail = len(
            [i for i in end_state_df.index if i in tail['id'] and end_state_df['label'][i] == "DISK"]) / len(
            end_state_df[end_state_df['label'] == "DISK"])
        escape_from_tail = len(
            [i for i in end_state_df.index if i in tail['id'] and end_state_df['label'][i] == "ESCAPE"]) / len(
            end_state_df[end_state_df['label'] == "ESCAPE"])
        si_and_tail_data[title] = [sum(si_mass), si_ang_mom, sum(tail_mass), tail_ang_mom, si_pct_tar_silicate,
                                   si_pct_tar_iron,
                                   si_pct_imp_silicate, si_pct_imp_iron, tail_pct_tar_silicate, tail_pct_tar_iron,
                                   tail_pct_imp_silicate, tail_pct_imp_iron, si_pct_silicate, si_pct_iron,
                                   tail_pct_silicate, tail_pct_iron, si_planet, si_disk, si_escape, tail_planet,
                                   tail_disk,
                                   tail_escape, planet_from_si, disk_from_si, escape_from_si, planet_from_tail,
                                   disk_from_tail, escape_from_tail]

        dfs[title].append(df)
        target_coms[title].append(target_com)
        impactor_coms[title].append(impactor_com)
        sis[title].append(si)
        tails[title].append(tail)
        not_classifieds[title].append(not_classified)

    index_headers = ["Secondary Impactor Mass", "Secondary Impactor Angular Momentum",
                     "Tail Mass", "Tail Angular Momentum", "SI Frac. Target Silicate",
                     "SI Frac. Target Iron", "SI Frac. Impactor Silicate",
                     "SI Frac. Impactor Iron", "Tail Frac. Target Silicate",
                     "Tail Frac. Target Iron", "Tail Frac. Impactor Silicate",
                     "Tail Frac. Impactor Iron", "SI Frac. Silicate",
                     "SI Frac. Iron", "Tail Frac. Silicate",
                     "Tail Frac. Iron", "SI Frac. Planet", "SI Frac. Disk",
                     "SI Frac. Escape", "Tail Frac. Planet", "Tail Frac. Disk",
                     "Tail Frac. Escape", "Frac. Planet From SI",
                     "Frac. Disk From SI", "Frac. Escape From SI", "Frac. Planet From Tail",
                     "Frac. Disk From Tail", "Frac. Escape From Tail"]

    df_si_and_tail_data = pd.DataFrame(si_and_tail_data, index=index_headers)
    df_si_and_tail_data.to_csv("{}_secondary_impact_structure_data.csv".format(angle))
    to_save_path = "{}_secondary_impact_structures".format(angle)
    if not os.path.exists(to_save_path):
        os.mkdir(to_save_path)
    # __build_scene(iteration=iteration, dfs=dfs, sims=sims, titles=titles,
    #               target_coms=target_coms, impactor_coms=impactor_coms, to_save_path=to_save_path, times=times)
    __plot_secondary(iteration=iteration, dfs=dfs, sims=sims, titles=titles,
                     target_coms=target_coms, impactor_coms=impactor_coms, to_save_path=to_save_path, times=times,
                     sis=sis, tails=tails, not_classifieds=not_classifieds)


get_secondary_and_tail()
