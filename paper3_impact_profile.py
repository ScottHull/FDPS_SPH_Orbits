#!/usr/bin/env python3
import os
import csv
import shutil
from math import pi, asin, isnan, exp
import numpy as np
import pandas as pd
from random import randint
from statistics import mean
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LogNorm
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import multiprocessing as mp
from matplotlib.ticker import MaxNLocator
import string

from src.vapor import calc_vapor_mass_fraction_from_formatted
from src.geometry import get_impact_geometry_from_formatted, get_velocity_profile_from_formatted
from src.animate import animate
from src.identify import ParticleMap
from src.combine import CombineFile

base_path = "/home/theia/scotthull/Paper2_SPH/gi/"
num_processes = 600
iterations = [50, 100, 200, 500, 1000]
paths = [['500_mars', "Mars " + r"($b=0.73$)"]]

begin_iteration = 0
end_iteration = 100
increment = 5

for s, t in paths:
    cd = int(s.split("_")[0])
    times = []
    target_velocity = []
    impactor_velocity = []
    for iteration in np.arange(begin_iteration, end_iteration + increment, increment):
        path = base_path + "{}/{}".format(s, s)
        to_fname = "merged_{}_{}.dat".format(iteration, randint(0, 100000))
        cf = CombineFile(num_processes=num_processes, time=iteration, output_path=path, to_fname=to_fname)
        df = cf.combine_df()
        formatted_time = round(cf.sim_time * 0.000277778, 2)
        df['velocity'] = np.sqrt(df['vx'] ** 2 + df['vy'] ** 2 + df['z'] ** 2)
        target = df[df['tag'] < 2]
        impactor = df[df['tag'] >= 2]
        target_velocity.append(target['velocity'].mean())
        impactor_velocity.append(impactor['velocity'].mean())
        times.append(formatted_time)
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)
    ax.plot(times, np.array(target_velocity) / 1000, linewidth=2.0, label="Target")
    ax.plot(times, np.array(impactor_velocity) / 1000, linewidth=2.0, label="Impactor")
    ax.set_xlabel("Time (hours)")
    ax.set_ylabel("Velocity (km/s)")
    ax.grid()
    ax.legend()
    ax.set_title(f"{t} Impact Velocity Profile")
    plt.savefig("{}_impact_velocity_profile.png".format(s), dpi=300)
