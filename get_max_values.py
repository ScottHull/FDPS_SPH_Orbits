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

min_iteration = 0
max_iteration = 300
end_state_iteration = 1800
increment = 1
angle = 'b073'
cutoff_densities = [5, 500, 1000, 2000]
base_path = "/home/theia/scotthull/Paper1_SPH/gi/"


def get_all_sims(angle, high=True):
    fformat = "{}_{}_{}"
    tformat = "{}{}{}"
    names = []
    titles = []
    for runs in ["new", "old"]:
        n = "s"
        if runs == "old":
            n = "n"
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


def get_endstate(s):
    path = base_path + "{}/circularized_{}".format(s, s)
    df = pd.read_csv(path + "/{}.csv".format(end_state_iteration))
    disk = df[df['label'] == "DISK"]
    return disk['id'].tolist()


sims, titles = get_all_sims(angle=angle, high=True)
max_vals = dict(zip(titles, [{} for i in titles]))


def find_max(df, title, curr_maxes, iteration, time):
    for row in df.index:
        if df['pressure'][row] > curr_maxes[df['id'][row]]['pressure']:
            curr_maxes[df['id'][row]]['pressure'] = df['pressure'][row]
            curr_maxes[df['id'][row]]['time'] = time
            curr_maxes[df['id'][row]]['iteration'] = iteration
            curr_maxes[df['id'][row]]['density'] = df['density'][row]
            curr_maxes[df['id'][row]]['entropy'] = df['entropy'][row]
            curr_maxes[df['id'][row]]['temperature'] = df['temperature'][row]
            curr_maxes[df['id'][row]]['internal energy'] = df['internal energy'][row]
    return curr_maxes


def __run_proc(args):
    s, t = args
    endstate = get_endstate(s)
    number_processes = 200
    if "high" in s:
        number_processes = 500
    curr_max = {i: {'pressure': -1e99, "time": None, "iteration": None, "density": None, "entropy": None,
                    "internal energy": None} for i in endstate}
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
        os.remove(f)
        curr_max = find_max(disk, t, curr_max, iteration, formatted_time)
    max_vals[t] = curr_max
    f = "{}_{}_max_vals.csv".format(t, angle)
    if os.path.exists(f):
        os.remove(f)
    with open(f, 'w') as infile:
        infile.write(str(curr_max))
    infile.close()


def run():
    pool = mp.Pool(10)
    pool.map(__run_proc, [[s, t] for s, t in zip(sims, titles)])
    pool.close()
    pool.join()
    fname = "{}_max_vals.csv".format(angle)
    if os.path.exists(fname):
        os.remove(fname)
    with open(fname, "w") as f:
        f.write(str(max_vals))
    f.close()


run()
