import os
import csv
from math import pi, asin
import numpy as np
import pandas as pd
from statistics import mean
import matplotlib.pyplot as plt

from src.vapor import calc_vapor_mass_fraction_from_formatted
from src.geometry import get_impact_geometry_from_formatted, get_velocity_profile_from_formatted
from src.animate import animate


def get_time(f):
    formatted_time = None
    with open(f, 'r') as infile:
        reader = csv.reader(infile, delimiter="\t")
        formatted_time = float(next(reader)[0])
    infile.close()
    return round(formatted_time * 0.000277778, 2)  # seconds -> hours

def __get_vmf_timeplot_data(path, phase_path, start_iteration, end_iteration, increment):
    max_time = get_time(path + "/{}.csv".format(end_iteration))
    times = []
    vmfs = []
    disk_particle_count = []
    avg_disk_entropy = []
    for time in np.arange(start_iteration, end_iteration + increment, increment):
        f = path + "/{}.csv".format(time)
        times.append(get_time(f))
        df = pd.read_csv(f, skiprows=2)
        vmf = calc_vapor_mass_fraction_from_formatted(df=df, phase_path=phase_path) * 100.0
        vmfs.append(vmf)

        disk_particles = df[df['label'] == "DISK"]
        try:
            avg_disk_entropy_at_time = mean(disk_particles['entropy'])
        except:
            avg_disk_entropy_at_time = 0
        avg_disk_entropy.append(avg_disk_entropy_at_time)
        num_disk_particles = len(disk_particles['entropy'])
        disk_particle_count.append(num_disk_particles)
    return times, vmfs, disk_particle_count, avg_disk_entropy, max_time

def build_vmf_timeplots(meta, start_iteration, end_iteration, increment, label_header='label',
                        output_fname="vmf_profile_output.txt"):
    """
    Builds VMF timeseries plots from formatted outputs.
    :param names:
    :param paths:
    :return:
    """
    plt.style.use("dark_background")
    fig, axs = plt.subplots(2, 3, figsize=(16, 9), sharex='all', sharey='all',
                            gridspec_kw={"hspace": 0.10, "wspace": 0.10})
    axs = axs.flatten()
    axs[0].set_title("New EoS")
    axs[1].set_title("Old EoS")
    for ax in axs:
        ax.grid(alpha=0.4)
    for i in meta.keys():
        try:
            print("working on {}".format(i))
            n = meta[i]['name']
            p = meta[i]['path']
            if "new" in p:
                phase_path = "src/phase_data/forstSTS__vapour_curve.txt"
            else:
                phase_path = "src/phase_data/duniteN__vapour_curve.txt"
            times, vmfs, disk_particle_count, avg_disk_entropy, max_time = __get_vmf_timeplot_data(
                p, phase_path, start_iteration, end_iteration, increment)
            if "new" in n:
                axs[0].plot(times, vmfs, linewidth=2.0, label=n)
                axs[2].plot(times, avg_disk_entropy, linewidth=2.0, label=n)
                axs[4].plot(times, disk_particle_count, linewidth=2.0, label=n)
            else:
                axs[1].plot(times, vmfs, linewidth=2.0, label=n)
                axs[3].plot(times, avg_disk_entropy, linewidth=2.0, label=n)
                axs[5].plot(times, disk_particle_count, linewidth=2.0, label=n)
        except FileNotFoundError:
            print(i)
    for ax in axs:
        ax.legend(loc='upper left')
    plt.savefig("vmf_timeseries.png", format='png', dpi=200)


def build_impact_angle_geometries(meta, start_iteration, end_iteration, specified_imp_angle, increment=1):
    plt.style.use("dark_background")
    """
    :param names:
    :param paths:
    :param start_iteration:
    :param end_iteration:
    :param increment:
    :return:
    """
    imp_ang_fig, imp_ang_axs = plt.subplots(1, 2, figsize=(16, 9), sharex='all', sharey='all',
                            gridspec_kw={"hspace": 0.10, "wspace": 0.10})
    imp_ang_axs = imp_ang_axs.flatten()
    imp_ang_axs[0].set_title("Impact Parameter (New EoS)"), imp_ang_axs[1].set_title("Impact Parameter (Old EoS)")
    imp_ang_axs[0].set_ylabel("Impact Angle (deg)")

    for ax in imp_ang_axs:
        ax.axhline(asin(specified_imp_angle) * (180 / pi), linewidth=2.0, linestyle="--", c='red', label="{} deg.".format(round(asin(specified_imp_angle) * (180 / pi), 2)))
        ax.set_xlabel("Time (hrs)")
        ax.grid(alpha=0.4)

    for i in meta.keys():
        try:
            print("working on {}".format(i))
            n = meta[i]['name']
            p = meta[i]['path']
            times, imp_angles = [], []
            for time in np.arange(start_iteration, end_iteration + increment, increment):
                t = get_time(p + "/{}.csv".format(time))
                times.append(t)
                df = pd.read_csv(p + "/{}.csv".format(time), skiprows=2)
                imp_angles.append(get_impact_geometry_from_formatted(df=df, name=i, time=t, iteration=time))
            if "new" in n.lower():
                imp_ang_axs[0].plot(
                    times, imp_angles, linewidth=2.0, label=n
                )
            else:
                imp_ang_axs[1].plot(
                    times, imp_angles, linewidth=2.0, label=n
                )
        except FileNotFoundError:
            print(i)
        animate(
            start_time=start_iteration,
            end_time=end_iteration,
            interval=increment,
            path=i + "_tmp_geometry",
            filename="{}_imp_geometry.mp4".format(i)
        )
    for ax in imp_ang_axs:
        ax.legend(loc='upper left')
    plt.savefig("impact_angle_profile.png", format='png', dpi=200)

def build_impact_velocity_charts(meta, start_iteration, end_iteration, increment=1):
    imp_vel_fig, imp_vel_axs = plt.subplots(3, int(len(meta.keys()) / 3), figsize=(16, 32), sharex='all', sharey='all',
                                            gridspec_kw={"hspace": 0.10, "wspace": 0.10})
    imp_vel_axs = imp_vel_axs.flatten()

    fig_index = 0

    for ax in imp_vel_axs:
        ax.set_xlabel("Time (hrs)")
        ax.grid(alpha=0.4)

    for i in meta.keys():
        try:
            print("working on {}".format(i))
            n = meta[i]['name']
            p = meta[i]['path']
            times, imp_vels = [], []
            for time in np.arange(start_iteration, end_iteration + increment, increment):
                times.append(get_time(p + "/{}.csv".format(time)))
                df = pd.read_csv(p + "/{}.csv".format(time), skiprows=2)
                imp_vels.append(get_velocity_profile_from_formatted(df) / 1000)
            imp_vel_axs[fig_index].plot(
                times, imp_vels, linewidth=2.0, label=n
            )
            if fig_index % 2 == 0:
                imp_vel_axs[fig_index].set_ylabel("Impact Velocity (km/s)")
            fig_index += 1
        except FileNotFoundError:
            print(i)
    for ax in imp_vel_axs:
        ax.legend(loc='upper left')
    plt.savefig("impact_velocity_profile.png", format='png', dpi=200)

