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

base_path = "/home/theia/scotthull/Paper3_SPH/gi/"
num_processes = 600
paths = [
    ['500_mars_b073_2v_esc', "Mars (v_imp=2v_esc, b=0.73)"],
    ['500_mars_b073_1v_esc', "Mars (v_imp=1v_esc, b=0.73)"],
    ['500_mars_b050_1v_esc', "Mars (v_imp=1 v_esc, b=0.50)"]
]

begin_iteration = 0
end_iteration = 20
increment = 4

headers = ["id", "tag", "mass", "x", "y", "z", "vx", "vy", "vz", "density", "internal energy", "pressure",
                   "potential energy", "entropy", "temperature"]
mass_mars = 6.39e23  # kg

data = {
    i[1]: {
        'times': [],
        'target_velocity': [],
        'impactor_velocity': [],
    } for i in paths
}

def calculate_v_esc(R_tar=3363.53 * 1000, R_imp=336.04 * 1000, M_tot = 1.001 * mass_mars):
    """
    Calculate the mutual escape velocity.
    :return:
    """
    G = 6.67408e-11  # m3 kg-1 s-2
    return np.sqrt((2 * G * M_tot) / (R_tar + R_imp))

v_esc = calculate_v_esc()

for index, (s, t) in enumerate(paths):
    cd = int(s.split("_")[0])
    # times = []
    # target_velocity = []
    # impactor_velocity = []
    # imp_vel_from_imp = []
    # imp_vel_from_tar = []
    for index2, iteration in enumerate(np.arange(begin_iteration, end_iteration + increment, increment)):
        print(iteration)
        path = base_path + "{}/{}".format(s, s)
        to_fname = "merged_{}_{}.dat".format(iteration, randint(0, 100000))
        cf = CombineFile(num_processes=num_processes, time=iteration, output_path=path, to_fname=to_fname)
        df = cf.combine_df()
        print("# of particles: {}".format(len(df)))
        df.columns = headers
        formatted_time = round(cf.sim_time * 0.000277778, 2)
        df['velocity'] = np.sqrt(df['vx'] ** 2 + df['vy'] ** 2 + df['vz'] ** 2)
        target = df[df['tag'] < 2]
        impactor = df[df['tag'] >= 2]
        mean_target_velocity = target['velocity'].mean()
        mean_impactor_velocity = impactor['velocity'].mean()
        target_mass = target['mass'].sum()
        impactor_mass = impactor['mass'].sum()
        total_mass = target_mass + impactor_mass
        v_imp_from_tar = (total_mass * mean_target_velocity) / impactor_mass
        v_imp_from_imp = (total_mass * mean_impactor_velocity) / target_mass
        data[t]['times'].append(formatted_time)
        data[t]['target_velocity'].append(mean_target_velocity)
        data[t]['impactor_velocity'].append(mean_impactor_velocity)
        # target_velocity.append(mean_target_velocity)
        # impactor_velocity.append(mean_impactor_velocity)
        # times.append(formatted_time)
        # imp_vel_from_imp.append(v_imp_from_imp)
        # imp_vel_from_tar.append(v_imp_from_tar)
        print("Mean target velocity: {} km/s".format(mean_target_velocity / 1000))
        print("Mean impactor velocity: {} km/s".format(mean_impactor_velocity / 1000))
        print("Target mass: {} kg".format(target_mass))
        print("Impactor mass: {} kg".format(impactor_mass))

    # fig = plt.figure(figsize=(10, 10))
    # ax = fig.add_subplot(111)
    # ax.plot(times, np.array(target_velocity) / 1000, linewidth=1.0, color='red', label="Target")
    # ax.plot(times, np.array(impactor_velocity) / 1000, linewidth=1.0, color='blue', label="Impactor")
    # # scatter and annotate the points on top of the lines
    # ax.scatter(
    #     times, np.array(target_velocity) / 1000, marker='o', color='red', s=10
    # )
    # ax.scatter(
    #     times, np.array(impactor_velocity) / 1000, marker='o', color='blue', s=10
    # )
    # for i, txt in enumerate(times):
    #     ax.annotate(target_velocity[i] / 1000, (times[i], np.array(target_velocity)[i] / 1000), fontsize=8)
    #     ax.annotate(impactor_velocity[i] / 1000, (times[i], np.array(impactor_velocity)[i] / 1000), fontsize=8)
    # # ax.plot(times, np.array(imp_vel_from_tar) / 1000, linewidth=1.0, label="Impact Velocity (from target)")
    # # ax.plot(times, np.array(imp_vel_from_imp) / 1000, linewidth=1.0, label="Impact Velocity (from impactor)")
    # ax.set_xlabel("Time (hours)")
    # ax.set_ylabel("Velocity (km/s)")
    # ax.grid()
    # ax.legend()
    # ax.set_title(f"{t} Impact Velocity Profile")
    # plt.savefig("{}_impact_velocity_profile.png".format(s), dpi=200)


fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)
# get a series of unique colors of the length of the data keys array
colors = cm.jet(np.linspace(0, 1, len(data.keys())))
for index, p in enumerate(data.keys()):
    d = data[p]
    ax.plot(d['times'], np.array(d['target_velocity']) / 1000, linewidth=2.0, color=colors[index], label=f'{p} (Target)')
    ax.plot(d['times'], np.array(d['impactor_velocity']) / 1000, linewidth=2.0, color=colors[index], linestyle='dashed',
            label=f'{p} (Impactor)')
    # ax.plot(times, np.array(target_velocity) / 1000, linewidth=2.0, color='red', label="Target")
    # ax.plot(times, np.array(impactor_velocity) / 1000, linewidth=2.0, color='blue', label="Impactor")
# set a horizontal line at v_esc
ax.axhline(y=v_esc / 1000, color='black', linestyle='dashed', label=r"$v_{esc}$")
ax.set_xlabel("Time (hours)")
ax.set_ylabel("Velocity (km/s)")
ax.grid()
ax.legend()
ax.set_title(f"Impact Velocity Profiles")
plt.savefig("mars_impact_velocity_profile.png".format(s), dpi=200)
