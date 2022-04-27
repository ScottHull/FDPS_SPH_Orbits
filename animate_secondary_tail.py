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
si_iteration = 190
max_iteration = 1800
increment = 20
angle = 'b073'
cutoff_densities = [5, 500, 1000, 2000]
number_processes = 200
square_scale = 6e7
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
        to_fname = "merged_{}_{}.dat".format(si_iteration, randint(0, 100000))
        cf = CombineFile(num_processes=number_processes, time=si_iteration, output_path=path, to_fname=to_fname)
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

        si, tail, not_classified = get_secondary_imp_and_tail_particles(title, df)

        dfs[title].append(df)
        sis[title].append(si)
        tails[title].append(tail)
        not_classifieds[title].append(not_classified)

    return sis, tails


def plot_time(dfs, sis, tails, endstates, to_path, iteration, time):
    num_new = len([i for i in dfs.keys() if "n" in i])
    num_old = len([i for i in dfs.keys() if "o" in i])
    num_cols = max([num_new, num_old])
    plt.style.use("dark_background")
    fig, axs = plt.subplots(2, num_cols, figsize=(20, 10), sharex='all',
                            gridspec_kw={"hspace": 0.10, "wspace": 0.12})
    fig.patch.set_facecolor('xkcd:black')
    axs = axs.flatten()
    for ax in axs:
        ax.set_xlim(-square_scale, square_scale)
        ax.set_ylim(-square_scale, square_scale)
        ax.grid(alpha=0.1)
    index_new, index_old = 0, num_cols
    for t, df in dfs.items():
        to_index = index_new
        if "o" in t:
            to_index = index_old
        endstate = endstates[t]
        endstate_disk = endstate[endstate['label'] == "DISK"]
        curr_disk = df[df['id'].isin(endstate_disk.index.tolist())]

        si = sis[t][-1]
        tail = tails[t][-1]

        curr_si = df[df['id'].isin(si['id'].tolist())]
        curr_tail = df[df['id'].isin(tail['id'].tolist())]

        rest = df[~df['id'].isin(si['id'].tolist())]
        rest = df[~df['id'].isin(tail['id'].tolist())]
                
        # tail_in_disk = curr_disk[curr_disk['id'].isin(tail['id'].tolist())]
        # tail_not_in_disk = tail[~tail['id'].isin(tail_in_disk['id'].tolist())]
        # tail_not_in_disk = df[df['id'].isin(tail_not_in_disk['id'].to_list())]
        # disk_rest = curr_disk[~curr_disk['id'].isin(tail_in_disk['id'].tolist())]
        # disk_rest = disk_rest[~disk_rest['id'].isin(tail_not_in_disk['id'].tolist())]
        # not_disk = df[~df['id'].isin(endstate_disk.index.tolist())]
        # not_disk = not_disk[~not_disk['id'].isin(tail_in_disk['id'].tolist())]
        # axs[to_index].scatter(
        #     not_disk['x'], not_disk['y'], s=2, alpha=0.2, label="Other"
        # )
        # for d, label in zip([disk_rest, tail_in_disk, tail_not_in_disk], ["DISK NOT IN TAIL", "TAIL IN DISK", "TAIL NOT IN DISK"]):
        #     axs[to_index].scatter(
        #         d['x'], d['y'], marker=".", s=2, alpha=1.0, label=label
        #     )
        #     axs[to_index].set_title("{} {} hrs. ({})".format(t, time, iteration))
        for d, label in zip([curr_si, curr_tail, rest], ["si", "tail", "rest"]):
            axs[to_index].scatter(
                    d['x'], d['y'], marker=".", s=2, alpha=1.0, label=label
                )
            axs[to_index].set_title("{} {} hrs. ({})".format(t, time, iteration))
        if "o" in t:
            index_old += 1
        else:
            index_new += 1

    legend = axs[0].legend(loc='upper left', fontsize=8)
    for handle in legend.legendHandles:
        try:
            handle.set_sizes([30.0])
        except:
            pass

    plt.savefig(to_path + "/{}.png".format(iteration), format='png', dpi=200)


def get_end_states(angle, high):
    endstates = {}
    sims, titles = get_all_sims(angle, high)
    for s, t in zip(sims, titles):
        end_state_file = base_path + "{}/circularized_{}/{}.csv".format(s, s, max_iteration)
        end_state_df = pd.read_csv(end_state_file, index_col="id")
        endstates.update({t: end_state_df})
    return endstates

sis, tails = get_secondary_and_tail()
endstates = get_end_states(angle=angle, high=False)

def run_proc(args):
    iteration, to_path = args
    if not os.path.exists(to_path):
        os.mkdir(to_path)
    data = {}
    high = False
    sims, titles = get_all_sims(angle, high=high)
    formatted_time = None

    for s, t in zip(sims, titles):
        print("{} - {}".format(iteration, t))

        path = base_path + "{}/{}".format(s, s)
        to_fname = "merged_{}_{}.dat".format(iteration, randint(0, 100000))
        cf = CombineFile(num_processes=number_processes, time=iteration, output_path=path, to_fname=to_fname)
        combined_file = cf.combine()
        formatted_time = round(cf.sim_time * 0.000277778, 2)
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

        data.update({t: df})

    plot_time(
        dfs=data, sis=sis, tails=tails, endstates=endstates, to_path=to_path, iteration=iteration, time=formatted_time
    )

def run():
    image_path = "{}_tail_in_disk".format(angle)
    if not os.path.exists(image_path):
        os.mkdir(image_path)
    pool = mp.Pool(10)
    pool.map(run_proc, [[iteration, image_path] for iteration in
                        np.arange(min_iteration, max_iteration + increment, increment)])
    pool.close()
    pool.join()


run()
