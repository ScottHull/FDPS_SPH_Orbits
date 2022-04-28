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


def get_end_states(angle, high):
    endstates = {}
    sims, titles = get_all_sims(angle, high)
    for s, t in zip(sims, titles):
        end_state_file = base_path + "{}/circularized_{}/{}.csv".format(s, s, max_iteration)
        end_state_df = pd.read_csv(end_state_file, index_col="id")
        endstates.update({t: end_state_df})
    return endstates


def get_secondary_and_tail():
    sims, titles = get_all_sims(angle=angle, high=False)
    impact_angles = {}
    times = {}
    sis = {}
    tails = {}
    not_classifieds = {}
    dfs = {}
    for output_name, title in zip(sims, titles):
        if title not in impact_angles.keys():
            sis.update({title: []})
            tails.update({title: []})
            not_classifieds.update({title: []})
        if title not in dfs.keys():
            dfs.update({title: []})
            times.update({title: []})
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

    return sis, tails, not_classifieds


def reformat_dict(d: dict):
    index_l = list(d[list(d.keys())[0]].keys())
    data = {}
    for t in d.keys():
        data.update({t: [d[t][h] for h in index_l]})
    return pd.DataFrame(data, index=index_l)


sis, tails, not_classifieds = get_secondary_and_tail()
endstates = get_end_states(angle=angle, high=False)

def profile_time():

    sims, titles = get_all_sims(angle, False)

    data = {}
    for s, t in zip(sims, titles):
        data.update({t: {}})
        endstate = endstates[t]
        si = sis[t][-1]
        tail = tails[t][-1]
        not_classified = not_classifieds[t][-1]

        planet, disk, escape = endstate[endstate['label'] == "PLANET"], endstate[endstate['label'] == "DISK"], \
                               endstate[endstate['label'] == "ESCAPE"]
        not_classified['radius'] = [(i ** 2 + j ** 2 + k ** 2) ** (1 / 2) for i, j, k in
                                    zip(not_classified['x'], not_classified['y'], not_classified['z'])]

        si_in_planet = si[si['id'].isin(planet.index.tolist())]
        si_in_disk = si[si['id'].isin(disk.index.tolist())]
        si_in_escape = si[si['id'].isin(escape.index.tolist())]
        tail_in_planet = tail[tail['id'].isin(planet.index.tolist())]
        tail_in_disk = tail[tail['id'].isin(disk.index.tolist())]
        tail_in_escape = tail[tail['id'].isin(escape.index.tolist())]

        data[t].update({
            "% PLANET FROM SI": len(si_in_planet) / len(planet) * 100.0,
            "% DISK FROM SI": len(si_in_disk) / len(disk) * 100.0,
            "% ESCAPE FROM SI": len(si_in_escape) / len(escape) * 100.0,
            "% PLANET FROM TAIL": len(tail_in_planet) / len(planet) * 100.0,
            "% DISK FROM TAIL": len(tail_in_disk) / len(disk) * 100.0,
            "% ESCAPE FROM TAIL": len(tail_in_escape) / len(escape) * 100.0,

        })

        # for particles not in tail or SI
        nc_planet, nc_disk, nc_escape = not_classified[not_classified['id'].isin(planet.index)], \
                                        not_classified[not_classified['id'].isin(disk.index)], \
                                        not_classified[not_classified['id'].isin(escape.index)]

        # disk from inside planet
        nc_disk_interior = nc_disk[nc_disk['radius'] < 1.5e7]
        # disk from far away planet (outside tail)
        nc_disk_exterior = nc_disk[nc_disk['radius'] >= 1.5e7]

        data[t].update({
            "% DISK FROM INSIDE PLANET": len(nc_disk_interior) / len(disk) * 100.0,
            "% DISK FROM OUTSIDE PLANET (not tail/si)": len(nc_disk_exterior) / len(disk) * 100.0,
        })

        data[t].update({
            "DISK % CHECK": data[t]["% DISK FROM SI"] + data[t]["% DISK FROM TAIL"] +
                            data[t]["% DISK FROM INSIDE PLANET"] + data[t]["% DISK FROM OUTSIDE PLANET (not tail/si)"]
        })

        # check1 = disk[disk.index.isin(si_in_disk['id'].tolist())]
        # check2 = disk[disk.index.isin(tail_in_disk['id'].tolist())]
        # check3 = disk[disk.index.isin(nc_disk_interior['id'].tolist())]
        # check4 = disk[disk.index.isin(nc_disk_exterior['id'].tolist())]
        #
        # data[t].update({
        #     "CHECK1": len(check1) / len(disk) * 100.0,
        #     "CHECK2": len(check2) / len(disk) * 100.0,
        #     "CHECK3": len(check3) / len(disk) * 100.0,
        #     "CHECK4": len(check4) / len(disk) * 100.0,
        # })
        # data[t].update({
        #     "CHECK % CHECK": data[t]["CHECK1"] + data[t]["CHECK2"] + data[t]["CHECK3"] + data[t]["CHECK4"]
        # })

    d = reformat_dict(data)
    d.to_csv("{}_secondary_impact_struc_data.csv".format(angle))


def plot_iteration(iteration, time, dfs, tail_in_disk, planet_in_disk, exterior_in_disk, disk_from_si, other, to_path):
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
        axs[to_index].scatter(
            other[t]['x'], other[t]['y'], s=2, alpha=0.2, label="Other"
        )
        for d, label in zip(
                [tail_in_disk[t], planet_in_disk[t], disk_from_si[t], exterior_in_disk[t]],
                ["Disk From Tail", "Disk from Planet", "Disk from SI", "Disk from Exterior Debris"]):
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



def run_proc(args):
    iteration, to_path = args

    dfs, tail_in_disks, planet_in_disks, exterior_in_disks, other, si_in_disks, times = {}, {}, {}, {}, {}, {}, {}

    sims, titles = get_all_sims(angle, False)

    data = {}
    for s, t in zip(sims, titles):

        print("{} - {}".format(iteration, t))

        path = base_path + "{}/{}".format(s, s)
        to_fname = "merged_{}_{}.dat".format(iteration, randint(0, 100000))
        cf = CombineFile(num_processes=number_processes, time=iteration, output_path=path, to_fname=to_fname)
        combined_file = cf.combine()
        formatted_time = round(cf.sim_time * 0.000277778, 2)
        times.update({t: formatted_time})
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

        endstate = endstates[t]
        si = sis[t][-1]
        tail = tails[t][-1]
        not_classified = not_classifieds[t][-1]

        planet, disk, escape = endstate[endstate['label'] == "PLANET"], endstate[endstate['label'] == "DISK"], \
                               endstate[endstate['label'] == "ESCAPE"]
        not_classified['radius'] = [(i ** 2 + j ** 2 + k ** 2) ** (1 / 2) for i, j, k in
                                    zip(not_classified['x'], not_classified['y'], not_classified['z'])]

        si_in_planet = si[si['id'].isin(planet.index.tolist())]
        si_in_disk = si[si['id'].isin(disk.index.tolist())]
        si_in_escape = si[si['id'].isin(escape.index.tolist())]
        tail_in_planet = tail[tail['id'].isin(planet.index.tolist())]
        tail_in_disk = tail[tail['id'].isin(disk.index.tolist())]
        tail_in_escape = tail[tail['id'].isin(escape.index.tolist())]

        # for particles not in tail or SI
        nc_planet, nc_disk, nc_escape = not_classified[not_classified['id'].isin(planet.index)], \
                                        not_classified[not_classified['id'].isin(disk.index)], \
                                        not_classified[not_classified['id'].isin(escape.index)]

        # disk from inside planet
        nc_disk_interior = nc_disk[nc_disk['radius'] < 1.5e7]
        # disk from far away planet (outside tail)
        nc_disk_exterior = nc_disk[nc_disk['radius'] >= 1.5e7]

        other_particles = df[~df['id'].isin(tail_in_disk['id'].tolist())]
        other_particles = other_particles[~other_particles['id'].isin(nc_disk_interior['id'].tolist())]
        other_particles = other_particles[~other_particles['id'].isin(nc_disk_exterior['id'].tolist())]

        dfs.update({t: df})
        tail_in_disks.update({t: tail_in_disk})
        planet_in_disks.update({t: nc_disk_interior})
        exterior_in_disks.update({t: nc_disk_exterior})
        si_in_disks.update({t: si_in_disk})
        other.update({t: other_particles})

    plot_iteration(iteration=iteration, time=times, dfs=dfs, tail_in_disk=tail_in_disks, planet_in_disk=planet_in_disks,
                   exterior_in_disk=exterior_in_disks, disk_from_si=si_in_disks, other=other, to_path=to_path)


def run():
    image_path = "{}_tail_in_disk_simple".format(angle)
    if not os.path.exists(image_path):
        os.mkdir(image_path)
    pool = mp.Pool(10)
    pool.map(run_proc, [[iteration, image_path] for iteration in
                        np.arange(min_iteration, max_iteration + increment, increment)])
    pool.close()
    pool.join()


run()

