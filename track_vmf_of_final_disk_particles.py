#!/usr/bin/env python3
import os
import csv
from statistics import mean
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from src.vapor import calc_vapor_mass_fraction_with_circularization_from_formatted, \
    calc_vapor_mass_fraction_without_circularization_from_formatted

# use seaborn colorblind palette
plt.style.use('seaborn-colorblind')

runs = ['500_b073_new', '500_b073_old', '1000_b073_new', '1000_b073_old', '2000_b073_new', '2000_b073_old']

min_iteration = 0
max_iteration = 1800

base_path = "/home/theia/scotthull/Paper1_SPH/gi"
new_phase_path = "src/phase_data/forstSTS__vapour_curve.txt"
old_phase_path = "src/phase_data/duniteN__vapour_curve.txt"

for run in runs:
    plt.style.use('seaborn-colorblind')
    phase_path = new_phase_path if 'new' in run else old_phase_path
    times, iterations, vmfs_w_circ, vmfs_wo_circ, entropies, entropies_w_circ, temperatures, pressures, velocities = [], [], [], [], [], [], [], [], []
    path = base_path + f"/{run}/circularized_{run}"
    path2 = base_path + f"/{run}/{run}_reports"
    end_time_df = pd.read_csv(
        os.path.join(path, f"{max_iteration}.csv"),
    )
    end_time_disk = end_time_df[end_time_df["label"] == "DISK"]
    end_time_disk = end_time_disk[end_time_disk["circ_entropy_delta"] <= 5000]
    end_time_particle_ids = end_time_disk["id"].values

    # loop through all files in path, where the iteration is the file name minus the extension
    for file in os.listdir(path):
        f = path + "/" + file
        iteration = file.split(".")[0]
        if int(iteration) >= min_iteration and int(iteration) <= max_iteration:
            try:
                # read the report file as a pandas df and get the time
                f2 = path2 + "/" + iteration + ".csv"
                df2 = pd.read_csv(f2, sep=",")
                time = df2["TIME_HRS"][0]
                df = pd.read_csv(f, sep=",")
                df["entropy_w_circ"] = df["entropy"] + df["circ_entropy_delta"]
                # use vx, vy, and vz to calculate the magnitude of the velocity
                df["velocity"] = (df["vx"] ** 2 + df["vy"] ** 2 + df["vz"] ** 2) ** 0.5
                # get the initial conditions for the hydrodynamics
                # limit df to only the particles that were in the disk at the end of the simulation
                disk = df[df["id"].isin(end_time_particle_ids)]
                filtered_disk = disk[disk['circ_entropy_delta'] <= 5000]
                vmf_w_circ = calc_vapor_mass_fraction_with_circularization_from_formatted(filtered_disk,
                                                                                          phase_path=phase_path,
                                                                                          restrict_df=False)
                vmf_wo_circ = calc_vapor_mass_fraction_without_circularization_from_formatted(disk,
                                                                                              phase_path=phase_path,
                                                                                              restrict_df=False)
                times.append(time)
                iterations.append(iteration)
                entropies.append(disk['entropy'].mean())
                entropies_w_circ.append(filtered_disk['entropy_w_circ'].mean())
                temperatures.append(disk['temperature'].mean())
                pressures.append(disk['pressure'].mean())
                vmfs_w_circ.append(vmf_w_circ * 100)
                vmfs_wo_circ.append(vmf_wo_circ * 100)
                velocities.append(disk['velocity'].mean())
            except:
                pass

    # sort the times and vmfs by time
    times, iterations, vmfs_w_circ, vmfs_wo_circ, entropies, entropies_w_circ, temperatures, pressures, velocities = zip(
        *sorted(zip(times, iterations, vmfs_w_circ, vmfs_wo_circ, entropies, entropies_w_circ, temperatures, pressures, velocities)))

    # get each value associated with the maximum pressure (time of impact)
    min_velocity = min(velocities[0:5])
    min_velocity_index = velocities.index(min_velocity)
    min_velocity_pressure = pressures[min_velocity_index]
    min_velocity_time = times[min_velocity_index]
    min_velocity_iteration = iterations[min_velocity_index]
    min_velocity_vmf_w_circ = vmfs_w_circ[min_velocity_index]
    min_velocity_vmf_wo_circ = vmfs_wo_circ[min_velocity_index]
    min_velocity_entropy = entropies[min_velocity_index]
    min_velocity_entropy_w_circ = entropies_w_circ[min_velocity_index]
    min_velocity_temperature = temperatures[min_velocity_index]

    # plot the results
    fig, axs = plt.subplots(1, 5, figsize=(24, 8))
    axs = axs.flatten()
    axs[0].plot(times, vmfs_w_circ, label="with circularization")
    axs[0].plot(times, vmfs_wo_circ, label="without circularization")
    axs[0].set_xlabel("Time (hrs)")
    axs[0].set_ylabel("Vapor Mass Fraction (%)")
    axs[0].set_title(f"Vapor Mass Fraction - {run}")

    axs[1].plot(times, entropies, label="without circularization")
    axs[1].plot(times, entropies_w_circ, label="with circularization")
    axs[1].set_xlabel("Time (hrs)")
    axs[1].set_ylabel("Avg. Disk Entropy (J/K)")
    axs[1].set_title(f"Entropy - {run}")

    axs[2].plot(times, temperatures)
    axs[2].set_xlabel("Time (hrs)")
    axs[2].set_ylabel("Temperature (K)")
    axs[2].set_title(f"Avg. Disk Temperature - {run}")

    axs[3].plot(times, np.array(pressures) / 10 ** 9)
    axs[3].set_xlabel("Time (hrs)")
    axs[3].set_ylabel("Temperature (K)")
    axs[3].set_title(f"Avg. Disk Pressure - {run}")

    axs[4].plot(times, np.array(velocities) / 1000)
    axs[4].set_xlabel("Time (hrs)")
    axs[4].set_ylabel("Velocity (km/s)")
    axs[4].set_title(f"Avg. Disk Velocity - {run}")

    # annotate in upper-right corner of the first plot
    axs[0].annotate(f"Time of impact: {min_velocity_time:.2f} hrs\nP: {min_velocity_pressure / 10 ** 9} GPa\n"
                    f"S_wo_circ: {min_velocity_entropy}\n"
                    f"S_w_circ: {min_velocity_entropy_w_circ}\nVMF_w_circ: {min_velocity_vmf_w_circ}\n"
                    f"VMF_wo_circ: {min_velocity_vmf_wo_circ}\nT: {min_velocity_temperature}\nvel: {min_velocity} km/s",
                    xy=(0.95, 0.75), xycoords='axes fraction', horizontalalignment='right', verticalalignment='top')

    for ax in axs:
        ax.legend()
        ax.grid()
        ax.axvline(x=min_velocity_time, color="red", linestyle="--", linewidth=2, label="impact time")

    plt.savefig(f"vmf_{run}_final_disk_particles.png")

    # use dark background color palette
    plt.style.use('dark_background')
    # get the final iteration file and load into a DataFrame
    df = pd.read_csv(base_path + f"/{run}/circularized_{run}" + f"/{min_velocity_iteration}.csv")
    # scatter the iteration according to when the label is PLANET, DISK, or ESCAPE at the final iteration
    fig = plt.figure(figsize=(12, 12))
    ax = fig.add_subplot(111)
    final_planet = end_time_df[end_time_df['label'] == 'PLANET']
    final_disk = end_time_df[end_time_df['label'] == 'DISK']
    final_escape = end_time_df[end_time_df['label'] == 'ESCAPE']
    ax.scatter(df[df.isin(final_planet)].dropna()['x'], df[df.isin(final_planet)].dropna()['y'], s=2, label='PLANET')
    ax.scatter(df[df.isin(final_disk)].dropna()['x'], df[df.isin(final_disk)].dropna()['y'], s=2, label='DISK')
    ax.scatter(df[df.isin(final_escape)].dropna()['x'], df[df.isin(final_escape)].dropna()['y'], s=2,label='ESCAPE')
    ax.set_title(f"Final Disk Particles - {run} - {min_velocity_iteration} ({min_velocity_time:.2f} hrs)")

    # annotate in upper-right corner of the first plot that this is colored based on the final iteration
    ax.annotate("Color coded based on final iteration",
                    xy=(0.95, 0.75), xycoords='axes fraction', horizontalalignment='right', verticalalignment='top')

    plt.savefig(f"{run}_min_vel_time_plotted.png")
