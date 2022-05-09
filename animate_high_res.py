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
from src.animate import animate

plt.rcParams.update({'font.size': 8, })
plt.style.use("dark_background")

min_iteration = 0
end_iteration = 1800
increment = 20
square_scale = 6e7
base_path = "/home/theia/scotthull/Paper1_SPH/gi/"
to_path = "animate_high_res"
sim, title = "5_b073_new_high", "5b073n-high"

iterations = list(np.arange(min_iteration, end_iteration + increment, increment))
if not os.path.exists(to_path):
    os.mkdir(to_path)

headers = ["id", "tag", "mass", "x", "y", "z", "vx", "vy", "vz", "density", "internal energy", "pressure",
                   "potential energy", "entropy", "temperature"]

endstate = pd.read_csv(base_path + "{}/circularized_{}/{}.csv".format(sim, sim, end_iteration), index_col='id')
end_planet, end_disk, end_escape = endstate[endstate['label'] == "PLANET"], endstate[endstate['label'] == "DISK"], \
                                   endstate[endstate['label'] == "ESCAPE"]

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

def run_proc(args):
    iteration = args[0]
    if "high" in sim:
        number_processes = 500
    else:
        number_processes = 200
    file_format = "results.{}_{}_{}.dat"
    p2 = base_path + "{}/{}".format(sim, sim)
    base_file = file_format.format(
        str(iteration).zfill(5),
        str(number_processes).zfill(5),
        str(0).zfill(5)
    )
    formatted_time = get_time(p2 + "/" + base_file)

    try:
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111, aspect='equal')
        # df = pd.read_csv(base_path + "{}/circularized_{}/{}.csv".format(sim, sim, iteration))
        # planet = df[df['id'].isin(end_planet.index.tolist())]
        # disk = df[df['id'].isin(end_disk.index.tolist())]
        # escape = df[df['id'].isin(end_escape.index.tolist())]
        path = base_path + "{}/{}".format(sim, sim)
        to_fname = "merged_{}_{}.dat".format(iteration, randint(0, 100000))
        cf = CombineFile(num_processes=number_processes, time=iteration, output_path=path, to_fname=to_fname)
        combined_file = cf.combine()
        formatted_time = round(cf.sim_time * 0.000277778, 2)
        f = os.getcwd() + "/{}".format(to_fname)
        df = pd.read_csv(f, skiprows=2, header=None, delimiter="\t", names=headers)
        planet = df[df['id'].isin(end_planet.index.tolist())]
        disk = df[df['id'].isin(end_disk.index.tolist())]
        escape = df[df['id'].isin(end_escape.index.tolist())]
        os.remove(f)
        for i, label in zip([planet, disk, escape], ['Planet', 'Disk', 'Escape']):
            ax.scatter(
                i['x'], i['y'], alpha=1, marker=".", s=1, label=label
            )
        ax.set_xlim(-square_scale, square_scale)
        ax.set_ylim(-square_scale, square_scale)
        ax.set_title("{} - {} hrs ({})".format(title, formatted_time, iteration))

        legend = ax.legend(loc='upper left', fontsize=14)
        for handle in legend.legendHandles:
            try:
                handle.set_sizes([50.0])
            except:
                pass
        plt.savefig(to_path + "/{}.png".format(iteration), format='png', dpi=200)
    except Exception as e:
        print(e)
        pass

def run():
    pool = mp.Pool(10)
    pool.map(run_proc, [[iteration] for iteration in iterations])
    pool.close()
    pool.join()


run()

animate(
    start_time=min_iteration,
    end_time=end_iteration,
    interval=increment,
    path=to_path,
    fps=10,
    filename="animate_5b073n-high.mp4",
)
