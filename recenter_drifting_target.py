#!/usr/bin/env python3
import os
import csv
import numpy as np
import pandas as pd
from random import randint
import matplotlib.pyplot as plt
plt.style.use("dark_background")

from src.combine import CombineFile

from_path = "/home/theia/scotthull/Paper1_SPH/tar-imp/5_new_high/target2_5_new_high"
iteration = 40
num_proc = 500

def calc_center_of_mass(mass, x, y, z):
    total_mass = sum(mass)
    com_x = sum(mass * x) / total_mass
    com_y = sum(mass * y) / total_mass
    com_z = sum(mass * z) / total_mass
    return com_x, com_y, com_z

fname = "merged_{}.dat".format(randint(0, 100000))
c = CombineFile(num_processes=num_proc, time=iteration, output_path=from_path, to_fname=fname).combine()
df = pd.read_csv(fname, skiprows=2, header=None, delimiter="\t")
mass, x, y, z = df[2], df[3], df[4], df[5]
com_x, com_y, com_z = calc_center_of_mass(mass, x, y, z)
x = x - com_x
y = y - com_y
z = z - com_z
distance_re = (((x ** 2 + y ** 2 + z ** 2) ** (1/2)) / (6371 * 1000))

fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
ax.scatter(
    distance_re, df['density'], s=2
)
ax.set_xlabel(r"Radius ($R_{\bigoplus}$)")
ax.set_ylabel("Density")
ax.grid(alpha=0.4)
plt.savefig("recentered_target.png", format='png', dpi=200)

os.remove(fname)
