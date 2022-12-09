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
square_scale = 2e7 / 1000
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
    time, iterations, mean_disk_vel = [], [], []
    for iteration in np.arange(min_iteration, max_vel_profile_iteration, increment):
        path = base_path + f"{run}/{run}"
        to_fname = "merged_{}_{}.dat".format(iteration, randint(0, 100000))
        cf = CombineFile(num_processes=number_processes, time=iteration, output_path=path, to_fname=to_fname)
        c = cf.combine()
        formatted_time = round(cf.sim_time * 0.000277778, 2)

        combined_file = pd.read_csv(to_fname, skiprows=2, header=None, delimiter="\t")

        os.remove(to_fname)

        final_disk_particles = combined_file[combined_file[0].isin(end_time_particle_ids)]

        id, x, y, z, vx, vy, vz = combined_file[0], combined_file[3], combined_file[4], combined_file[5], combined_file[6], combined_file[7], combined_file[8]
        id_disk, x_disk, y_disk, z_disk, vx_disk, vy_disk, vz_disk = final_disk_particles[0], final_disk_particles[3], final_disk_particles[4], final_disk_particles[5], final_disk_particles[6], final_disk_particles[7], final_disk_particles[8]
        velocity = np.sqrt(vx_disk**2 + vy_disk**2 + vz_disk**2)

        time.append(formatted_time)
        mean_disk_vel.append(velocity.mean())
        
        fig, ax = plt.subplots(1, 2, figsize=(10, 5))
        axs = ax.flatten()
        axs[0].scatter(
            np.array(x) / 1000, np.array(y) / 1000, s=0.1, c="k", alpha=1
        )
        axs[0].scatter(
            np.array(x_disk) / 1000, np.array(y_disk) / 1000, s=0.1, c="r", alpha=1
        )
        axs[0].set_xlabel("x")
        axs[0].set_ylabel("y")
        axs[0].set_title(f"Time: {formatted_time} hours (iteration {iteration})")
        axs[1].plot(
            iterations, mean_disk_vel, c="k", alpha=1
        )
        axs[0].set_xlim(-square_scale, square_scale)
        axs[0].set_ylim(-square_scale, square_scale)
        axs[1].set_xlim(min_iteration, max_vel_profile_iteration)
        axs[1].set_ylim(0, 10000)
        axs[1].set_xlabel("Iteration")
        axs[1].set_ylabel("Mean Final Disk Particle Velocity (km/s)")
        axs[1].set_title(f"Time: {formatted_time} hours (iteration {iteration})")

        # in the upper right corner of axs[1], annotate the mean velocity
        axs[1].annotate(
            f"Mean Velocity: {round(mean_disk_vel[-1], 2)} km/s",
            xy=(0.95, 0.95),
            xycoords="axes fraction",
            horizontalalignment="right",
            verticalalignment="top",
        )

        plt.savefig(f"{run}_vel_profile/{iteration}.png")

    animate(min_iteration, max_iteration, increment, f"{run}_vel_profile", filename=f"{run}_disk_vel_profile.mp4", fps=10)


