#!/usr/bin/env python3
import os
import csv
import numpy as np
import pandas as pd
from random import randint
import matplotlib.pyplot as plt

from src.combine import CombineFile

plt.style.use("dark_background")

from_path = "/home/theia/scotthull/Paper1_SPH/tar-imp/5_new_high/target_5_new_high"
to_path = "animate_target"
start_iteration = 0
end_iteration = 60
inc = 5
num_proc = 400
square_scale = 1e5 * 1000

def calc_center_of_mass(mass, x, y, z):
    total_mass = sum(mass)
    com_x = (mass * x) / total_mass
    com_y = (mass * y) / total_mass
    com_z = (mass * z) / total_mass
    return com_x, com_y, com_z

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

times, coms = [], []
for iteration in np.arange(start_iteration, end_iteration + inc, inc):
    fname = "merged_{}.dat".format(randint(0, 100000))
    c = CombineFile(num_processes=num_proc, time=iteration, output_path=from_path, to_fname=fname).combine()
    times.append(get_time(fname))
    df = pd.read_csv(fname, skiprows=2, header=None, delimiter="\t")
    mass, x, y, z = df[2], df[3], df[4], df[5]
    com_x, com_y, com_z = calc_center_of_mass(mass, x, y, z)
    coms.append(np.linalg.norm(com_x, com_y, com_y))

fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
ax.plot(
    times, coms, linewidth=2.0
)
ax.set_xlabel("Time (hrs)")
ax.set_ylabel("Center of Mass")
ax.set_title("Center of Mass vs Time (Mode 2 Target)")
ax.grid(alpha=0.4)
plt.savefig("mode2_com_tracking.png", format='png', dpi=200)
