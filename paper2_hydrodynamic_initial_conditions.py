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

# use colorblind-friendly colors from seaborn
# plt.style.use('seaborn-colorblind')

runs = [('/home/theia/scotthull/Paper1_SPH/gi/500_b073_new', 'Canonical'), ('', 'Half-Earths'), ('', 'Mars')]

min_iteration = 0
max_iteration = 1800
max_vel_profile_iteration = 60
increment = 1
number_processes = 200
square_scale = 2e7 / 1000
new_phase_path = "src/phase_data/forstSTS__vapour_curve.txt"
old_phase_path = "src/phase_data/duniteN__vapour_curve.txt"

# make a plot with 3 columns and n rows
fig, axs = plt.subplots(len(runs), 3, figsize=(12, 8))
axs = axs.flatten()

for ax in axs:
    ax.grid(True)

for index, (run, verbose_run_name) in enumerate(runs):
    # if len(run) == 0, then skip this part of the loop
    if len(run) == 0:
        continue
    run_name = run.split('/')[-1]
    base_path = run.split(f"/{run_name}")[0] + "/"
    run = run_name
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
            'id', 'tag', 'mass', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'density', 'internal_energy', 'pressure',
            'potential_energy', 'entropy', 'temperature'
        ]
        combined_file.columns = final_disk_particles.columns

        vmf_final_disk = calc_vapor_mass_fraction_without_circularization_from_formatted(
            final_disk_particles[final_disk_particles['tag'] % 2 == 0], phase_path, restrict_df=False
        ) * 100
        vmf_all_silicate = calc_vapor_mass_fraction_without_circularization_from_formatted(
            combined_file[combined_file['tag'] % 2 == 0], phase_path, restrict_df=False) * 100
        final_disk_vmf.append(vmf_final_disk)
        all_silicate_vmf.append(vmf_all_silicate)

    gi_index = 0 + (index * len(runs))
    velocity_index = 1 + (index * len(runs))
    t_s_index = 2 + (index * len(runs))
    ic_iteration, ic_time = iterations[max(enumerate(mean_disk_vel), key=lambda x: x[1])[0]], \
        time[max(enumerate(mean_disk_vel), key=lambda x: x[1])[0]]

    # get the file corresponding to the initial condition
    to_fname = "merged_{}_{}.dat".format(ic_iteration, randint(0, 100000))
    cf = CombineFile(num_processes=number_processes, time=ic_iteration,
                     output_path=base_path + f'{run_name}/{run_name}', to_fname=to_fname)
    c = cf.combine()
    formatted_time = round(cf.sim_time * 0.000277778, 2)
    ic_file = pd.read_csv(to_fname, skiprows=2, header=None, delimiter="\t")
    os.remove(to_fname)
    ic_file.columns = [
        'id', 'tag', 'mass', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'density', 'internal_energy', 'pressure',
        'potential_energy', 'entropy', 'temperature'
    ]

    # get the file with the initial condition
    # scatter the simulation at the time of the initial condition
    # get subset of the data with particle ids that are not in the final disk
    not_disk_bound = ic_file[~ic_file["id"].isin(end_time_particle_ids)]
    disk_bound = ic_file[ic_file["id"].isin(end_time_particle_ids)]
    axs[gi_index].scatter(
        not_disk_bound["x"] / 1000, not_disk_bound["y"] / 1000, s=2, alpha=1, color='black', marker='.', label="Not Disk Bound"
    )
    axs[gi_index].scatter(
        disk_bound["x"] / 1000, disk_bound["y"] / 1000, s=2, alpha=1, color='red', marker='.', label="Disk Bound"
    )
    # set axis labels
    axs[gi_index].set_ylabel("y (km)")
    # annotate the upper left corner with the run name
    axs[gi_index].annotate(
        verbose_run_name, xy=(0.05, 0.95), xycoords="axes fraction", horizontalalignment="left", verticalalignment="top", fontweight="bold", fontsize=18
    )
    # axs[gi_index].set_xlim(-square_scale, square_scale)

    # plot the velocity profile
    axs[velocity_index].plot(
        time, np.array(mean_disk_vel) / 1000, linewidth=3, color='black', label="Mean Disk Velocity"
    )
    # set vertical line at the time of the initial condition
    axs[velocity_index].axvline(x=ic_time, color="black", linestyle="--", linewidth=2, label="Initial Condition")
    # set axis labels
    axs[velocity_index].set_ylabel("Velocity (km/s)")

    # plot the temperature profile
    axs[t_s_index].plot(
        time, mean_disk_temperature, linewidth=3, color='blue', label="Mean Disk Temperature"
    )
    # twin the x-axis and plot entropy
    axs[t_s_index].twinx().plot(
        time, mean_disk_entropy, linewidth=3, color='red', label="Mean Disk Entropy"
    )
    # set vertical line at the time of the initial condition
    axs[t_s_index].axvline(x=ic_time, color="black", linestyle="--", linewidth=2, label="Initial Condition")
    # set axis labels
    axs[t_s_index].set_ylabel("Temperature (K)")
    axs[t_s_index].twinx().set_ylabel("Entropy (J/K)")

    axs[-3].set_xlabel("x (km)")
    axs[-2].set_xlabel("Time (hrs.)")
    axs[-1].set_xlabel("Time (hrs.)")

# make tight layout with no hspace
plt.tight_layout(h_pad=0)
plt.savefig("paper2_hydrodynamic_initial_condition", dpi=300)
