#!/usr/bin/env python3
import os
import csv
import pandas as pd
import numpy as np
from math import atan, pi
from operator import contains
from random import randint
import multiprocessing as mp
import matplotlib.pyplot as plt

from src.combine import CombineFile

plt.style.use("dark_background")

base_path = "/home/theia/scotthull/Paper1_SPH/gi/"
angle = "b073"
cutoff_densities = [5, 500, 1000, 2000]
number_processes = 200
min_iteration = 200
max_iteration = 300
increment = 5
square_scale = 6e7

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
        n = "S"
        if runs == "old":
            n = "N"
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
    return np.array([x_center, y_center, z_center])


def get_impactor_com_particles(output_name):
    output_path = base_path + output_name + "/circularized_{}".format(output_name)
    df = pd.read_csv(output_path + "/{}.csv".format(min_iteration))
    impactor_iron = df[df['tag'] == 3]
    impactor_iron = impactor_iron[impactor_iron['radius'] > 1.5e7]
    return impactor_iron['id']


def get_angle(target_com, impactor_com):
    x_offset = target_com[0] - impactor_com[0]
    y_offset = target_com[1] - impactor_com[1]
    return atan(x_offset / y_offset) * (180 / pi)  # to degrees


def get_radial_vector_between_target_com_and_secondary_impactor_com(target_com, secondary_com):
    return np.array(secondary_com) - np.array(target_com)


def get_mean_vector(vectors):
    return [float(sum(l)) / len(l) for l in zip(*vectors)]


def __build_scene(iteration, times, dfs, sims, titles, impact_angles, target_coms, impactor_coms, mom_vectors,
                  to_save_path):
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
        df, impact_angle, target_com, impactor_com, time, mom_vector = dfs[t][-1], impact_angles[t][-1], \
                                                                       target_coms[t][-1], impactor_coms[t][-1], \
                                                                       times[t][-1], mom_vectors[t][-1]
        df = df[df['z'] < 0]
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

        r_vector = np.array(impactor_com) - np.array(target_com)
        axs[to_index].quiver(impactor_com[0], impactor_com[1], mom_vector[0], mom_vector[1], color='green', label="Spec. Mom. Vector")
        axs[to_index].quiver(target_com[0], target_com[1], r_vector[0], r_vector[1], color='yellow', label="Radial Vector")

        text = "|r\u20D7| = {:.2e}\n|v\u20D7| = {:.2e}\nr\u20D7 $\cdot$ v\u20D7 = {:.2e}\n".format(
            np.linalg.norm(r_vector), np.linalg.norm(mom_vector), np.dot(r_vector, mom_vector)
        )
        axs[to_index].text(square_scale - (0.5 * square_scale), -square_scale + (0.3 * square_scale),
                                text,
                                fontsize=10)

        axs[to_index].set_title("{} {} hrs. ({} - {} deg.)".format(t, time, iteration, round(impact_angle, 2)))
        if "old" in s:
            index_old += 2
        else:
            index_new += 2
    axs[0].legend(loc='upper left')
    plt.savefig(to_save_path + "/{}.png".format(iteration), format='png', dpi=200)


def get_impact_angle_with_velocity_vector():
    sims, titles = get_all_sims(high=False)
    impact_angles = {}
    target_coms = {}
    impactor_coms = {}
    mom_vectors = {}  # specific angular momentum
    imp_ids = {}
    times = {}
    r_dot_vs = {}
    for output_name, title in zip(sims, titles):
        imp_ids.update({title: get_impactor_com_particles(output_name)})
    for iteration in np.arange(min_iteration, max_iteration + increment, increment):
        dfs = {}
        for output_name, title in zip(sims, titles):
            if title not in impact_angles.keys():
                impact_angles.update({title: []})
                impactor_coms.update({title: []})
                target_coms.update({title: []})
                times.update({title: []})
                r_dot_vs.update({title: []})
                mom_vectors.update({title: []})
            if title not in dfs.keys():
                dfs.update({title: []})
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

            r_vector = impactor_com - target_com
            mean_impactor_velocity_vector = get_mean_vector(
                zip(impactor_iron['vx'], impactor_iron['vy'], impactor_iron['vz']))
            r_dot_v = np.dot(r_vector, mean_impactor_velocity_vector)

            imp_angle = get_angle(target_com, impactor_com)
            dfs[title].append(df)
            impact_angles[title].append(imp_angle)
            target_coms[title].append(target_com)
            impactor_coms[title].append(impactor_com)
            mom_vectors[title].append(mean_impactor_velocity_vector)
            r_dot_vs[title].append(r_dot_v)
        to_save_path = "{}_secondary_impact_angles_vectors".format(angle)
        if not os.path.exists(to_save_path):
            os.mkdir(to_save_path)
        __build_scene(iteration=iteration, dfs=dfs, sims=sims, titles=titles, impact_angles=impact_angles,
                      target_coms=target_coms, impactor_coms=impactor_coms, to_save_path=to_save_path,
                      mom_vectors=mom_vectors, times=times)

get_impact_angle_with_velocity_vector()
