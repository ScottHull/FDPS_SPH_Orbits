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
from src.vapor import calc_vapor_mass_fraction_without_circularization_from_formatted

runs = ['500_b073_new', '500_b073_old', '1000_b073_new', '1000_b073_old', '2000_b073_new', '2000_b073_old']

min_iteration = 0
max_iteration = 1800
max_vel_profile_iteration = 60
increment = 1
number_processes = 200
square_scale = 2e7 / 10 ** 6
base_path = "/home/theia/scotthull/Paper1_SPH/gi/"
new_phase_path = "src/phase_data/forstSTS__vapour_curve.txt"
old_phase_path = "src/phase_data/duniteN__vapour_curve.txt"

for run in runs:
    if not os.path.exists(f"{run}_vel_profile"):
        os.mkdir(f"{run}_vel_profile")
    circ_path = base_path + f"/{run}/circularized_{run}"
    end_time_df = pd.read_csv(
        os.path.join(circ_path, f"{max_iteration}.csv"),
    )
    end_time_disk = end_time_df[end_time_df["label"] == "DISK"]
    end_time_particle_ids = end_time_disk["id"].values
    time, iterations, mean_disk_vel, mean_disk_entropy, mean_disk_temperature, final_disk_vmf, all_silicate_vmf = [], [], [], [], [], [], []
    phase_path = new_phase_path if "new" in run else old_phase_path
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

        # replace final_disk_particles headers with the correct headers
        final_disk_particles.columns = [
            'id', 'tag', 'mass', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'density', 'internal_energy', 'pressure', 'potential_energy', 'entropy', 'temperature'
        ]
        combined_file.columns = final_disk_particles.columns

        vmf_final_disk = calc_vapor_mass_fraction_without_circularization_from_formatted(
            final_disk_particles[final_disk_particles['tag'] % 2 == 0], phase_path, restrict_df=False
        ) * 100
        vmf_all_silicate = calc_vapor_mass_fraction_without_circularization_from_formatted(
            combined_file[combined_file['tag'] % 2 == 0], phase_path, restrict_df=False) * 100
        final_disk_vmf.append(vmf_final_disk)
        all_silicate_vmf.append(vmf_all_silicate)


        fig, ax = plt.subplots(1, 3, figsize=(16, 5))
        axs = ax.flatten()
        axs[0].scatter(
            np.array(x) / 10 ** 6, np.array(y) / 10 ** 6, s=0.1, c="k", alpha=1
        )
        axs[0].scatter(
            np.array(x_disk) / 10 ** 6, np.array(y_disk) / 10 ** 6, s=0.1, c="r", alpha=1
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
        ax1_2.plot(iterations, np.array(mean_disk_temperature) / 1000, c="r", alpha=1, label="Mean Disk Temperature")
        ax1_2.plot(iterations, np.array(mean_disk_entropy) / 1000, c="green", alpha=1, label="Mean Disk Entropy")
        axs[1].plot([], [], c="r", alpha=1, label="Mean Disk Temperature")
        axs[1].plot([], [], c="green", alpha=1, label="Mean Disk Entropy")
        axs[1].set_ylim(0, 20)
        ax1_2.set_ylim(1, 6)
        axs[1].set_xlabel("Iteration")
        axs[1].set_ylabel("Mean Disk-Bound Particle Velocity (km/s)")
        ax1_2.set_ylabel("Mean Disk-Bound Particle Temperature/Entropy (/1000)")
        axs[1].set_title(f"{run} Disk-Bound Particles\nTime: {formatted_time} hours (iteration {iteration})")
        axs[1].grid(alpha=0.4)
        axs[1].legend(loc="lower right")

        axs[2].plot(
            iterations, np.array(final_disk_vmf), c="k", alpha=1, label="Disk-Bound VMF"
        )
        axs[2].plot(
            iterations, np.array(all_silicate_vmf), c="r", alpha=1, label="All Silicate VMF"
        )
        axs[2].set_xlim(min_iteration, max_vel_profile_iteration)
        axs[2].set_ylim(0, 30)
        axs[2].set_xlabel("Iteration")
        axs[2].set_ylabel("Vapor Mass Fraction (%)")
        axs[2].set_title(f"{run} Disk-Bound Particles\nTime: {formatted_time} hours (iteration {iteration})")
        axs[2].grid(alpha=0.4)
        axs[2].legend(loc="lower right")

        # increase subplot spacing
        fig.subplots_adjust(wspace=0.4)

        # in the upper right corner of axs[1], annotate the mean velocity
        axs[1].annotate(
            f"Mean Velocity: {round(mean_disk_vel[-1] / 1000, 2)} km/s\nMean Temperature: "
            f"{round(mean_disk_temperature[-1], 2)} K\nMean Entropy: {round(mean_disk_entropy[-1], 2)} "
            f"(J/K)\nDisk-Bound VMF: {round(final_disk_vmf[-1], 2)} %\n"
            f"All Silicate VMF: {round(all_silicate_vmf[-1], 2)} %",
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
    axs[0].set_title(f"{run} Disk-Bound Particles - Mean Velocity vs Time")
    axs[0].grid(alpha=0.4)
    axs[1].plot(
        time, np.array(mean_disk_entropy), c="k", alpha=1
    )
    axs[1].set_xlabel("Time (hours)")
    axs[1].set_ylabel("Entropy (J/K)")
    axs[1].set_title(f"{run} Disk-Bound Particles - Mean Entropy vs Time")
    axs[1].grid(alpha=0.4)
    axs[2].plot(
        time, np.array(mean_disk_temperature), c="k", alpha=1
    )
    axs[2].set_xlabel("Time (hours)")
    axs[2].set_ylabel("Temperature (K)")
    axs[2].set_title(f"{run} Disk-Bound Particles - Mean Temperature vs Time")
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
        f.write("Time (hours),Iteration,Mean Disk Velocity (km/s),Mean Disk Entropy (J/K),Mean Disk Temperature (K),All Silicates VMF (%),Disk-Bound Particles VMF (%)\n")
        for i in range(len(time)):
            f.write(f"{time[i]},{iterations[i]},{mean_disk_vel[i] / 1000},{mean_disk_entropy[i]},{mean_disk_temperature[i]},{all_silicate_vmf[i]},{final_disk_vmf[i]}\n")
    f.close()

    animate(min_iteration, max_iteration, increment, f"{run}_vel_profile", filename=f"{run}_disk_vel_profile.mp4",
            fps=10)
