#!/usr/bin/env python3
import os
import csv
import pandas as pd
import numpy as np
from random import randint
import multiprocessing as mp
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 18, })
plt.style.use("dark_background")

secondary_impact_times = {
    '5b073n': {
        'iteration': 220,
        'time': 6.11,
    },
    '500b073n': {
        'iteration': 255,
        'time': 7.08,
    },
    '1000b073n': {
        'iteration': 235,
        'time': 6.53,
    },
    '2000b073n': {
        'iteration': 205,
        'time': 5.69,
    },
    '5b073o': {
        'iteration': 235,
        'time': 6.53,
    },
    '500b073o': {
        'iteration': 260,
        'time': 7.22,
    },
    '1000b073o': {
        'iteration': 265,
        'time': 7.36,
    },
    '2000b073o': {
        'iteration': 245,
        'time': 6.81,
    },
}


def get_impact_momentum_vector_norm(df):
    mass, vx, vy, vz = df['mass'], df['vx'], df['vy'], df['vz']
    px, py, pz = mass * vx, mass * vy, mass * vz
    return np.array([px, py, pz])


def get_surface_normal_vector(df, secondary_imp_com):
    secondary_imp_com = np.array(secondary_imp_com)
    secondary_imp_com_to_earth_com = secondary_imp_com - np.array([0, 0, 0])
    surface_normal = secondary_imp_com_to_earth_com / np.linalg.norm(secondary_imp_com_to_earth_com)
    return surface_normal


def plot_impact_angles_vs_time():
    df = pd.read_csv("/Users/scotthull/Desktop/b073_secondary_impact_angles.csv", index_col="Iteration")
    df = df[df.index < 300]
    fig, axs = plt.subplots(1, 2, figsize=(16, 9), sharex="all", sharey='all')
    axs = axs.flatten()
    new_index, old_index = 0, 1
    for run in df.keys():
        to_index = new_index
        if "o" in run:
            to_index = old_index
        axs[to_index].plot(
            df.index, df[run], linewidth=2.0, label=run
        )
    for ax in axs:
        ax.grid(alpha=0.4)
        ax.set_xlabel("Iteration")
        ax.legend()
    axs[0].set_ylabel("Secondary Impact Angle (deg.)")
    axs[0].set_title("Stewart M-ANEOS")
    axs[1].set_title("GADGET M-ANEOS")
    plt.show()


def plot_moment_of_impact_vs_angle():
    df = pd.read_csv("/Users/scotthull/Desktop/b073_secondary_impact_angles.csv", index_col="Iteration")
    fig, axs = plt.subplots(1, 2, figsize=(16, 9), sharex="all", sharey='all')
    axs = axs.flatten()
    new_index, old_index = 0, 1
    for run in df.keys():
        to_index = new_index
        if "o" in run:
            to_index = old_index
        iteration, time = secondary_impact_times[run]['iteration'], secondary_impact_times[run]['time']
        angle = df[run][iteration]
        secondary_impact_times[run].update({'angle': angle})
        axs[to_index].scatter(
            time, angle, s=80, label=run
        )
    for ax in axs:
        ax.grid(alpha=0.4)
        ax.set_xlabel("Time (hrs)")
        ax.legend()
    axs[0].set_ylabel("Secondary Impact Angle (deg.)")
    axs[0].set_title("Stewart M-ANEOS")
    axs[1].set_title("GADGET M-ANEOS")
    plt.show()


# plot_impact_angles_vs_time()
plot_moment_of_impact_vs_angle()
