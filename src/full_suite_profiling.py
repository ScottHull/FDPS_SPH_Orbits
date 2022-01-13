import os
import csv
from math import pi, asin, isnan
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


def get_setup_file_data(path):
    d = {}
    with open(path, 'r') as infile:
        found_initial_setup = False
        for line in infile:
            if not found_initial_setup:
                if "INITIAL SETUP" in line:
                    found_initial_setup = True
            else:
                if len(line) == 1:
                    break
                else:
                    header, data = line.split(":")
                    d.update({header: data.replace("\n", "")})
    return d


def __get_vmf_timeplot_data(path, phase_path, start_iteration, end_iteration, increment, L_EM=3.5 * 10 ** 34):
    max_time = get_time(path + "/{}.csv".format(end_iteration))
    times = []
    vmfs = []
    disk_particle_count = []
    avg_disk_entropy = []
    spec_disk_ams = []
    for time in np.arange(start_iteration, end_iteration + increment, increment):
        f = path + "/{}.csv".format(time)
        try:
            times.append(get_time(f))
            df = pd.read_csv(f, skiprows=2)
            vmf = calc_vapor_mass_fraction_from_formatted(df=df, phase_path=phase_path) * 100.0
            vmfs.append(vmf)

            disk_particles = df[df['label'] == "DISK"]
            positions = list(zip(disk_particles['x'], disk_particles['y'], disk_particles['z']))
            velocities = list(zip(disk_particles['vx'], disk_particles['vy'], disk_particles['vz']))
            spec_disk_ams.append(sum([sum(np.cross(p, velocities[index])) for index, p in enumerate(positions)]) / L_EM)  # specific angular momentum of the disk
            try:
                avg_disk_entropy_at_time = mean(disk_particles['entropy'])
            except:
                avg_disk_entropy_at_time = 0
            avg_disk_entropy.append(avg_disk_entropy_at_time)
            num_disk_particles = len(disk_particles['entropy'])
            disk_particle_count.append(num_disk_particles)
        except:
            times.append(np.nan)
            vmfs.append(np.nan)
            disk_particle_count.append(np.nan)
            avg_disk_entropy.append(np.nan)
    return times, vmfs, disk_particle_count, avg_disk_entropy, max_time, spec_disk_ams


def build_vmf_timeplots(meta, start_iteration, end_iteration, increment, label_header='label',
                        output_fname="vmf_profile_output.txt", L_EM=3.5 * 10 ** 34):
    """
    Builds VMF timeseries plots from formatted outputs.
    :param names:
    :param paths:
    :return:
    """
    plt.style.use("dark_background")
    fig, axs = plt.subplots(4, 2, figsize=(16, 9), sharex='all',
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
            times, vmfs, disk_particle_count, avg_disk_entropy, max_time, spec_ang_mom = __get_vmf_timeplot_data(
                p, phase_path, start_iteration, end_iteration, increment, L_EM=L_EM)
            if "new" in n:
                axs[0].plot(times, vmfs, linewidth=2.0, label=n)
                axs[0].set_ylabel("VMF (%)")
                axs[2].plot(times, avg_disk_entropy, linewidth=2.0, label=n)
                axs[2].set_ylabel("Avg. Disk Entropy")
                axs[4].plot(times, disk_particle_count, linewidth=2.0, label=n)
                axs[4].set_ylabel("# Disk Particles")
                axs[6].plot(times, spec_ang_mom, linewidth=2.0, label=n)
                axs[6].set_ylabel("Specific Disk Ang. Mom. (L_EM)")
            else:
                axs[1].plot(times, vmfs, linewidth=2.0, label=n)
                axs[3].plot(times, avg_disk_entropy, linewidth=2.0, label=n)
                axs[5].plot(times, disk_particle_count, linewidth=2.0, label=n)
                axs[7].plot(times, spec_ang_mom, linewidth=2.0, label=n)
                axs[7].set_ylabel("Specific Disk Ang. Mom. (L_EM)")
        except FileNotFoundError:
            print(i)
    for ax in axs:
        ax.legend(loc='upper left')
    axs[0].set_ylim(0, 45)
    axs[1].set_ylim(0, 45)
    axs[2].set_ylim(0, 7000)
    axs[3].set_ylim(0, 7000)
    axs[4].set_ylim(0, 35000)
    axs[5].set_ylim(0, 35000)
    axs[4].set_xlabel("Time (hrs)")
    axs[5].set_xlabel("Time (hrs)")
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

    d = {}

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
            d.update({n: {"imp_angles": imp_angles, "times": times}})
            animate(
                start_time=start_iteration,
                end_time=end_iteration,
                interval=increment,
                path=i + "_tmp_geometry",
                filename="{}_imp_geometry.mp4".format(i)
            )
        except FileNotFoundError:
            print(i)

    imp_ang_fig, imp_ang_axs = plt.subplots(1, 2, figsize=(16, 9), sharex='all', sharey='all',
                                            gridspec_kw={"hspace": 0.10, "wspace": 0.10})
    imp_ang_axs = imp_ang_axs.flatten()
    imp_ang_axs[0].set_title("Impact Parameter (New EoS)"), imp_ang_axs[1].set_title("Impact Parameter (Old EoS)")
    imp_ang_axs[0].set_ylabel("Impact Angle (deg)")

    for ax in imp_ang_axs:
        ax.axhline(asin(specified_imp_angle) * (180 / pi), linewidth=2.0, linestyle="--", c='red',
                   label="{} deg.".format(round(asin(specified_imp_angle) * (180 / pi), 2)))
        ax.set_xlabel("Time (hrs)")
        ax.grid(alpha=0.4)

    for n in d.keys():
        if "new" in n.lower():
            imp_ang_axs[0].plot(
                d[n]["times"], d[n]["imp_angles"], linewidth=2.0, label=n
            )
        else:
            imp_ang_axs[1].plot(
                d[n]["times"], d[n]["imp_angles"], linewidth=2.0, label=n
            )
    for ax in imp_ang_axs:
        ax.legend(loc='upper left')
    plt.savefig("impact_angle_profile.png", format='png', dpi=200)


def build_impact_velocity_charts(meta, start_iteration, end_iteration, increment=1):
    plt.style.use("dark_background")
    imp_vel_fig, imp_vel_axs = plt.subplots(len(meta.keys()), 1, figsize=(16, 32), sharex='all', sharey='all',
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
            specified_imp_vel = float(meta[i]['setup']["ESCAPE VELOCITY"])
            times, imp_vels = [], []
            for time in np.arange(start_iteration, end_iteration + increment, increment):
                times.append(get_time(p + "/{}.csv".format(time)))
                df = pd.read_csv(p + "/{}.csv".format(time), skiprows=2)
                imp_vels.append(get_velocity_profile_from_formatted(df) / 1000)
            imp_vel_axs[fig_index].plot(
                times, imp_vels, linewidth=2.0,
                label="{} (Max: {})".format(n, round(max([x for x in imp_vels if isnan(x) == False]), 3))
            )
            imp_vel_axs[fig_index].axhline(specified_imp_vel / 1000, color='red', linewidth=2.0, linestyle="--",
                                           label="Specified ({} km/s)".format(round(specified_imp_vel / 1000, 3)))
            imp_vel_axs[fig_index].set_ylabel("Impact Velocity (km/s)")
            fig_index += 1
        except FileNotFoundError:
            print(i)
    for ax in imp_vel_axs:
        ax.legend(loc='upper left')
    plt.savefig("impact_velocity_profile.png", format='png', dpi=200)


def map_disk_to_phase_profile(meta, end_iteration):
    plt.style.use("dark_background")
    new_phase_path = "src/phase_data/forstSTS__vapour_curve.txt"
    old_phase_path = "src/phase_data/duniteN__vapour_curve.txt"
    new_phase_df = pd.read_fwf(new_phase_path, skiprows=1,
                               names=["temperature", "density_sol_liq", "density_vap", "pressure",
                                      "entropy_sol_liq", "entropy_vap"])
    old_phase_df = pd.read_fwf(old_phase_path, skiprows=1,
                               names=["temperature", "density_sol_liq", "density_vap", "pressure",
                                      "entropy_sol_liq", "entropy_vap"])
    fig, axs = plt.subplots(len(meta.keys()), 1, figsize=(16, 32), sharey='all',
                            gridspec_kw={"hspace": 0.10, "wspace": 0.10})
    axs = axs.flatten()
    fig_index = 0
    for ax in axs:
        ax.set_xlim(0, 15000)
        # ax.set_ylim(0, 15000)
        ax.grid(alpha=0.4)

    for i in meta.keys():
        try:
            n = meta[i]['name']
            p = meta[i]['path']
            df = pd.read_csv(p + "/{}.csv".format(end_iteration), skiprows=2)
            disk = df[df['tag'] % 2 == 0]
            disk = disk[disk['label'] == "DISK"]
            if "new" in n.lower():
                axs[fig_index].plot(
                    new_phase_df['entropy_vap'],
                    new_phase_df['temperature'],
                    linewidth=2.0,
                    label="Vapor"
                )
                axs[fig_index].plot(
                    new_phase_df['entropy_sol_liq'],
                    new_phase_df['temperature'],
                    linewidth=2.0,
                    label="Liquid",
                )
            else:
                axs[fig_index].plot(
                    old_phase_df['entropy_vap'],
                    old_phase_df['temperature'],
                    linewidth=2.0,
                    label="Vapor"
                )
                axs[fig_index].plot(
                    old_phase_df['entropy_sol_liq'],
                    old_phase_df['temperature'],
                    linewidth=2.0,
                    label="Liquid",
                )
            axs[fig_index].scatter(
                disk['entropy'],
                disk['temperature'],
                s=2,
                c='#fa8174',
                label="{} disk particles".format(n)
            )
            axs[fig_index].set_ylabel("Temperature")
            axs[fig_index].legend(loc='upper left')
            fig_index += 1
        except Exception as e:
            print("Problem!", e)
    axs[-1].set_xlabel("Entropy")
    axs[0].set_title("Disk Particles on Phase Curve")
    plt.savefig("disk_on_phase_curve.png", format='png', dpi=200)


def map_disk_to_phase_profile_eos_charts(meta, end_iteration):
    plt.style.use("dark_background")
    new_phase_path = "src/phase_data/forstSTS__vapour_curve.txt"
    old_phase_path = "src/phase_data/duniteN__vapour_curve.txt"
    new_phase_df = pd.read_fwf(new_phase_path, skiprows=1,
                               names=["temperature", "density_sol_liq", "density_vap", "pressure",
                                      "entropy_sol_liq", "entropy_vap"])
    old_phase_df = pd.read_fwf(old_phase_path, skiprows=1,
                               names=["temperature", "density_sol_liq", "density_vap", "pressure",
                                      "entropy_sol_liq", "entropy_vap"])
    fig, axs = plt.subplots(1, 2, figsize=(16, 9), sharey='all',
                            gridspec_kw={"hspace": 0.10, "wspace": 0.10})
    axs = axs.flatten()
    axs[0].set_title("New EoS")
    axs[1].set_title("Old EoS")
    axs[0].plot(
        new_phase_df['entropy_vap'],
        new_phase_df['temperature'],
        linewidth=2.0,
        label="Vapor"
    )
    axs[0].plot(
        new_phase_df['entropy_sol_liq'],
        new_phase_df['temperature'],
        linewidth=2.0,
        label="Liquid",
    )
    axs[1].plot(
        old_phase_df['entropy_vap'],
        old_phase_df['temperature'],
        linewidth=2.0,
        label="Vapor"
    )
    axs[1].plot(
        old_phase_df['entropy_sol_liq'],
        old_phase_df['temperature'],
        linewidth=2.0,
        label="Liquid",
    )
    axs[0].set_ylabel("Temperature")
    for ax in axs:
        ax.set_xlim(0, 15000)
        # ax.set_ylim(0, 15000)
        ax.grid(alpha=0.4)
        ax.set_xlabel("Entropy")

    for i in meta.keys():
        try:
            fig_index = None
            n = meta[i]['name']
            p = meta[i]['path']
            df = pd.read_csv(p + "/{}.csv".format(end_iteration), skiprows=2)
            disk = df[df['tag'] % 2 == 0]
            disk = disk[disk['label'] == "DISK"]
            if "new" in n.lower():
                fig_index = 0
            else:
                fig_index = 1
            axs[fig_index].scatter(
                disk['entropy'],
                disk['temperature'],
                s=2,
                label="{} disk particles".format(n)
            )
        except Exception as e:
            print("Problem!", e)
    for ax in axs:
        ax.legend(loc='upper left')
    plt.savefig("disk_on_phase_curve_same_eos_plots.png", format='png', dpi=200)
