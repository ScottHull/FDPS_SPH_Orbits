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

sims, titles = get_all_sims(angle=angle, high=True)
max_vals = dict(zip(titles, [{} for i in titles]))

def find_max(df, title, curr_maxes, iteration, time):
    max_rho = max(df['density'])
    max_pressure = max(df['density'])
    max_internal_energy = max(df['internal energy'])
    max_temperature = max(df['temperature'])
    max_entropy = max(df['entropy'])

    for i, j in zip([max_rho, max_pressure, max_internal_energy, max_temperature, max_entropy],
                    ['density', 'pressure', 'internal energy', 'temperature', 'entropy']):
        if j not in curr_maxes[title]:
            curr_maxes[title][j] = {"iteration": None, "time": None, "value": -1e99}
        if i > curr_maxes[title][j]['value']:
            curr_maxes[j]['value'] = i
            curr_maxes[j]['iteration'] = iteration
            curr_maxes[j]['time'] = time
    return curr_maxes

def __run_proc(args):
    s, t = args
    number_processes = 200
    if "high" in s:
        number_processes = 500
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
        os.remove(f)
        find_max(df, t, max_vals, iteration, formatted_time)

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
