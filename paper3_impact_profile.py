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
end_iteration = 20
increment = 4

headers = ["id", "tag", "mass", "x", "y", "z", "vx", "vy", "vz", "density", "internal energy", "pressure",
                   "potential energy", "entropy", "temperature"]

for s, t in paths:
    cd = int(s.split("_")[0])
    times = []
    target_velocity = []
    impactor_velocity = []
    imp_vel_from_imp = []
    imp_vel_from_tar = []
    for iteration in np.arange(begin_iteration, end_iteration + increment, increment):
        print(iteration)
        path = base_path + "{}/{}".format(s, s)
        to_fname = "merged_{}_{}.dat".format(iteration, randint(0, 100000))
        cf = CombineFile(num_processes=num_processes, time=iteration, output_path=path, to_fname=to_fname)
        df = cf.combine_df()
        df.columns = headers
        formatted_time = round(cf.sim_time * 0.000277778, 2)
        df['velocity'] = np.sqrt(df['vx'] ** 2 + df['vy'] ** 2 + df['z'] ** 2)
        target = df[df['tag'] < 2]
        impactor = df[df['tag'] >= 2]
        mean_target_velocity = target['velocity'].mean()
        mean_impactor_velocity = impactor['velocity'].mean()
        target_mass = target['mass'].sum()
        impactor_mass = impactor['mass'].sum()
        total_mass = target_mass + impactor_mass
        v_imp_from_tar = (total_mass * mean_target_velocity) / impactor_mass
        v_imp_from_imp = (total_mass * mean_impactor_velocity) / target_mass
        target_velocity.append(mean_target_velocity)
        impactor_velocity.append(mean_impactor_velocity)
        times.append(formatted_time)
        imp_vel_from_imp.append(v_imp_from_imp)
        imp_vel_from_tar.append(v_imp_from_tar)
        print("Mean target velocity: {} km/s".format(mean_target_velocity / 1000))
        print("Mean impactor velocity: {} km/s".format(mean_impactor_velocity / 1000))
        print("Target mass: {} kg".format(target_mass))
        print("Impactor mass: {} kg".format(impactor_mass))
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)
    ax.plot(times, np.array(target_velocity) / 1000, linewidth=2.0, label="Target")
    ax.plot(times, np.array(impactor_velocity) / 1000, linewidth=2.0, label="Impactor")
    ax.plot(times, np.array(imp_vel_from_tar) / 1000, linewidth=2.0, label="Impact Velocity (from target)")
    ax.plot(times, np.array(imp_vel_from_imp) / 1000, linewidth=2.0, label="Impact Velocity (from impactor)")
    ax.set_xlabel("Time (hours)")
    ax.set_ylabel("Velocity (km/s)")
    ax.grid()
    ax.legend()
    ax.set_title(f"{t} Impact Velocity Profile")
    plt.savefig("{}_impact_velocity_profile.png".format(s), dpi=200)
