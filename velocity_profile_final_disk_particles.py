#!/usr/bin/env python3
import os
import csv
from random import randint
from statistics import mean
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from src.combine import CombineFile
from src.animate import animate

runs = ['500_b073_new', '500_b073_old', '1000_b073_new', '1000_b073_old', '2000_b073_new', '2000_b073_old']

min_iteration = 0
max_iteration = 1800
max_vel_profile_iteration = 100
increment = 1
number_processes = 200
square_scale = 2e7 / 10 ** 6
base_path = "/home/theia/scotthull/Paper1_SPH/gi/"

for run in runs:
    if not os.path.exists(f"{run}_vel_profile"):
        os.mkdir(f"{run}_vel_profile")
    circ_path = base_path + f"/{run}/circularized_{run}"
    end_time_df = pd.read_csv(
        os.path.join(circ_path, f"{max_iteration}.csv"),
    )
    end_time_disk = end_time_df[end_time_df["label"] == "DISK"]
    end_time_particle_ids = end_time_disk["id"].values
    time, iterations, mean_disk_vel, mean_disk_entropy, mean_disk_temperature = [], [], [], [], []
    for iteration in np.arange(min_iteration, max_vel_profile_iteration, increment):
        iterations.append(iteration)
        path = base_path + f"{run}/{run}"
        to_fname = "merged_{}_{}.dat".format(iteration, randint(0, 100000))
        cf = CombineFile(num_processes=number_processes, time=iteration, output_path=path, to_fname=to_fname)
        c = cf.combine()
        formatted_time = round(cf.sim_time * 0.000277778, 2)

        combined_file = pd.read_csv(to_fname, skiprows=2, header=None, delimiter="\t")

        os.remove(to_fname)

        final_disk_particles = combined_file[combined_file[0].isin(end_time_particle_ids)]

        id, x, y, z, vx, vy, vz, entropy, temperature = combined_file[0], combined_file[3], combined_file[4], \
        combined_file[5], combined_file[6], combined_file[7], combined_file[8], combined_file[13], combined_file[14]

        id_disk, x_disk, y_disk, z_disk, vx_disk, vy_disk, vz_disk, entropy_disk, temperature_disk = \
        final_disk_particles[0], final_disk_particles[3], final_disk_particles[4], final_disk_particles[5], \
        final_disk_particles[6], final_disk_particles[7], final_disk_particles[8], final_disk_particles[13], \
        final_disk_particles[14]

        velocity = np.sqrt(vx_disk ** 2 + vy_disk ** 2 + vz_disk ** 2)

        time.append(formatted_time)
        mean_disk_vel.append(velocity.mean())
        mean_disk_entropy.append(entropy_disk.mean())
        mean_disk_temperature.append(temperature_disk.mean())

        fig, ax = plt.subplots(1, 2, figsize=(10, 5))
        axs = ax.flatten()
        axs[0].scatter(
            np.array(x) / 10 ** 6, np.array(y) / 1000, s=0.1, c="k", alpha=1
        )
        axs[0].scatter(
            np.array(x_disk) / 10 ** 6, np.array(y_disk) / 1000, s=0.1, c="r", alpha=1
        )
        axs[0].set_xlabel("x (1000 km)")
        axs[0].set_ylabel("y (1000 km)")
        axs[0].set_title(f"{run} - Time: {formatted_time} hours (iteration {iteration})")
        axs[1].plot(
            iterations, np.array(mean_disk_vel) / 1000, c="k", alpha=1, label="Mean Disk Velocity"
        )
        axs[0].set_xlim(-square_scale, square_scale)
        axs[0].set_ylim(-square_scale, square_scale)
        axs[1].set_xlim(min_iteration, max_vel_profile_iteration)
        # plot the temperature on the right y-axis
        ax1_2 = axs[1].twinx()
        ax1_2.plot(iterations, np.array(mean_disk_temperature), c="r", alpha=1, label="Mean Disk Temperature")
        axs[1].set_ylim(0, 20)
        ax1_2.set_ylim(1000, 6000)
        axs[1].set_xlabel("Time (hrs)")
        axs[1].set_ylabel("Mean Final Disk Particle Velocity (km/s)")
        axs[1].set_title(f"{run} Final Disk Particles\nTime: {formatted_time} hours (iteration {iteration})")
        axs[1].grid(alpha=0.4)

        # in the upper right corner of axs[1], annotate the mean velocity
        axs[1].annotate(
            f"Mean Velocity: {round(mean_disk_vel[-1] / 1000, 2)} km/s",
            xy=(0.95, 0.95),
            xycoords="axes fraction",
            horizontalalignment="right",
            verticalalignment="top",
        )

        plt.savefig(f"{run}_vel_profile/{iteration}.png")

    fig, axs = plt.subplots(1, 3, figsize=(15, 5))
    axs = axs.flatten()
    axs[0].plot(
        time, np.array(mean_disk_vel) / 1000, c="k", alpha=1
    )
    axs[0].set_xlabel("Time (hours)")
    axs[0].set_ylabel("Velocity (km/s)")
    axs[0].set_title(f"{run} Final Disk Particles - Mean Velocity vs Time")
    axs[0].grid(alpha=0.4)
    axs[1].plot(
        time, np.array(mean_disk_entropy), c="k", alpha=1
    )
    axs[1].set_xlabel("Time (hours)")
    axs[1].set_ylabel("Entropy (J/K)")
    axs[1].set_title(f"{run} Final Disk Particles - Mean Entropy vs Time")
    axs[1].grid(alpha=0.4)
    axs[2].plot(
        time, np.array(mean_disk_temperature), c="k", alpha=1
    )
    axs[2].set_xlabel("Time (hours)")
    axs[2].set_ylabel("Temperature (K)")
    axs[2].set_title(f"{run} Final Disk Particles - Mean Temperature vs Time")
    axs[2].grid(alpha=0.4)

    for ax in axs:
        # plot a vertical line at the time of max velocity on each plot
        ax.axvline(
            time[mean_disk_vel.index(max(mean_disk_vel))],
            linewidth=2.0,
            linestyle="--",
        )
        # annotate each plot with the max velocity and the y value at the time of max velocity
    axs[0].annotate(
        f"Max Velocity: {round(max(mean_disk_vel) / 1000, 2)} km/s\nMax Velocity Time: "
        f"{round(time[mean_disk_vel.index(max(mean_disk_vel))], 2)} hours\nEntropy at Max Velocity: "
        f"{round(mean_disk_entropy[mean_disk_vel.index(max(mean_disk_vel))], 2)} J/K\nTemperature at Max Velocity: "
        f"{round(mean_disk_temperature[mean_disk_vel.index(max(mean_disk_vel))], 2)} K\nFinal Entropy: {mean_disk_entropy[-1]} J/K\nFinal Temperature: {mean_disk_temperature[-1]} K",
        xy=(0.95, 0.75),
        xycoords="axes fraction",
        horizontalalignment="right",
        verticalalignment="top",
    )

    plt.savefig(f"{run}_vel_profile_all_iterations.png")

    # output to a CSV file
    with open(f"{run}_vel_profile.csv", "w") as f:
        f.write("Time (hours),Iteration,Mean Disk Velocity (km/s),Mean Disk Entropy (J/K),Mean Disk Temperature (K)\n")
        for i in range(len(time)):
            f.write(f"{time[i]},{iterations[i]},{mean_disk_vel[i] / 1000},{mean_disk_entropy[i]},{mean_disk_temperature[i]}\n")
    f.close()

    animate(min_iteration, max_iteration, increment, f"{run}_vel_profile", filename=f"{run}_disk_vel_profile.mp4",
            fps=10)
