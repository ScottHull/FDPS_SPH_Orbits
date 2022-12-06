#!/usr/bin/env python3
import os
import csv
from statistics import mean
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
    phase_path = new_phase_path if 'new' in run else old_phase_path
    times, vmfs_w_circ, vmfs_wo_circ, entropies, entropies_w_circ, temperatures, pressures = [], [], [], [], [], [], []
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
                # use vx, vy, and vz to calculate the magnitude of the velocity
                df["velocity"] = (df["vx"] ** 2 + df["vy"] ** 2 + df["vz"] ** 2) ** 0.5
                # get the initial conditions for the hydrodynamics
                # limit df to only the particles that were in the disk at the end of the simulation
                disk = df[df["id"].isin(end_time_particle_ids)]
                filtered_disk = disk[disk['circ_entropy_delta'] <= 5000]
                filtered_disk['entropy_w_circ'] = filtered_disk['entropy'] + filtered_disk['circ_entropy_delta']
                vmf_w_circ = calc_vapor_mass_fraction_with_circularization_from_formatted(filtered_disk,
                                                                                          phase_path=phase_path,
                                                                                          restrict_df=False)
                vmf_wo_circ = calc_vapor_mass_fraction_without_circularization_from_formatted(disk,
                                                                                              phase_path=phase_path,
                                                                                              restrict_df=False)
                times.append(time)
                entropies.append(disk['entropy'].mean())
                entropies_w_circ.append(filtered_disk['entropy_w_circ'].mean())
                temperatures.append(disk['temperature'].mean())
                vmfs_w_circ.append(vmf_w_circ * 100)
                vmfs_wo_circ.append(vmf_wo_circ * 100)
            except:
                pass

    # sort the times and vmfs by time
    times, vmfs_w_circ, vmfs_wo_circ, entropies, entropies_w_circ, temperatures, pressures = zip(
        *sorted(zip(times, vmfs_w_circ, vmfs_wo_circ, entropies, entropies_w_circ, temperatures, pressures)))

    # get the time associated with the max pressure
    max_pressure_time = times[pressures.index(max(pressures))]
    # get the vmfs, entropies, and temperatures at the max pressure time
    vmf_w_circ_at_max_pressure = vmfs_w_circ[pressures.index(max(pressures))]
    vmf_wo_circ_at_max_pressure = vmfs_wo_circ[pressures.index(max(pressures))]
    entropy_at_max_pressure = entropies[pressures.index(max(pressures))]
    entropy_w_circ_at_max_pressure = entropies_w_circ[pressures.index(max(pressures))]
    temperature_at_max_pressure = temperatures[pressures.index(max(pressures))]
    pressure_at_max_pressure = max(pressures)

    # plot the results
    fig, axs = plt.subplots(1, 3, figsize=(16, 9))
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

    # annotate max pressure time at top right
    axs[0].annotate(f"Time of impact: {max_pressure_time:.2f} hrs\nP: {pressure_at_max_pressure}\n"
                    f"S_wo_circ: {entropy_at_max_pressure}\n"
                    f"S_w_circ: {entropy_w_circ_at_max_pressure}\nVMF_w_circ: {vmf_w_circ_at_max_pressure}\n"
                    f"VMF_wo_circ: {vmf_wo_circ_at_max_pressure}\nT: {temperature_at_max_pressure}", xy=(0.95, 0.95),
                    xycoords='axes fraction',
                    horizontalalignment='right', verticalalignment='top')

    for ax in axs:
        ax.legend()
        ax.grid()

    plt.savefig(f"vmf_{run}_final_disk_particles.png")
