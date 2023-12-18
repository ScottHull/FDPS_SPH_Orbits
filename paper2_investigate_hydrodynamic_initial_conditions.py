#!/usr/bin/env python3
import os
import csv
import string
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

runs = [
    ('/home/theia/scotthull/Paper1_SPH/gi/500_b073_new', 'Canonical', 25),
    ('/home/theia/scotthull/Paper2_SPH/gi/500_half_earths', 'Half-Earths', 17),
    # ('', 'Mars')
]

min_iteration = 0
max_iteration = 1800
max_vel_profile_iteration = 60
increment = 1
number_processes = 200
square_scale = 2e7 / 1000
new_phase_path = "src/phase_data/forstSTS__vapour_curve.txt"
old_phase_path = "src/phase_data/duniteN__vapour_curve.txt"

fig, axs = plt.subplots(4, 2, figsize=(10, 20))
axs = axs.flatten()

for index, (run, verbose_run_name, iteration) in enumerate(runs):
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
    final_disk_particles['velocity'] = np.sqrt(final_disk_particles['vx'] ** 2 + final_disk_particles['vy'] ** 2 + final_disk_particles['vz'] ** 2)
    df2 = final_disk_particles[final_disk_particles['tag'] % 2 == 0]

    vmf_final_disk = calc_vapor_mass_fraction_without_circularization_from_formatted(
        df2, phase_path, restrict_df=False
    ) * 100

    axs[index].scatter(
        df2['velocity'] / 1000, df2['vmf_wo_circ'] * 100, s=5, label=verbose_run_name
    )
    axs[index].axvline(mean(df2['velocity'] / 1000), color='black', linestyle='--', label=f"Mean velocity: {mean(df2['velocity'] / 1000):.2f} km/s")
    axs[index].axhline(df2['vmf_wo_circ'].sum() / len(df2) * 100, color='black', linestyle='--', label=f"Mean VMF: {mean(df2['vmf_wo_circ'] * 100):.2f} %")
    axs[index].set_title(f"Disk-bound particles at initial condition ({run_name})")

    # on the bottom row, plot a CDF of the VMFs
    sorted_vmf = df2['vmf_wo_circ'].sort_values()
    cdf = sorted_vmf.rank(method='average', pct=True)
    axs[index + 2].plot(sorted_vmf * 100, cdf, linewidth=2.0)
    axs[index + 2].axvline(df2['vmf_wo_circ'].sum() / len(df2) * 100, color='black', linestyle='--',
                       label=f"Mean VMF: {mean(df2['vmf_wo_circ'] * 100):.2f} %")

    # plot a PDF of the VMFs
    axs[index + 4].hist(sorted_vmf * 100, bins=100, density=True)
    axs[index + 4].axvline(df2['vmf_wo_circ'].sum() / len(df2) * 100, color='black', linestyle='--',
                       label=f"Mean VMF: {mean(df2['vmf_wo_circ'] * 100):.2f} %")
    axs[index + 4].text(
        0.05, 0.90, f"vmf@0%: {len(df2[df2['vmf_wo_circ'] == 0])}\n"
                    f"0<vmf<1%: {len(df2[(df2['vmf_wo_circ'] > 0) & (df2['vmf_wo_circ'] < .01)])}\n"
                    f"1<vmf<99.99%: {len(df2[(df2['vmf_wo_circ'] > .01) & (df2['vmf_wo_circ'] < 0.9999)])}\n"
                    f"vmf@100%: {len(df2[df2['vmf_wo_circ'] > 0.9999])}",
        transform=axs[index + 4].transAxes, verticalalignment='top'
    )
    axs[index + 6].scatter(
        df2['entropy'], df2['temperature'], s=5, marker="."
    )

for ax in axs:
    ax.grid()
    ax.legend()
for ax in axs[:2]:
    ax.set_xlabel("Velocity (km/s)")
    ax.set_ylabel("Vapor Mass Fraction (%)")
    # ax.set_yscale('log')
for ax in axs[2:4]:
    ax.set_xlabel("Vapor Mass Fraction (%)")
    ax.set_ylabel("CDF")
    # ax.set_xscale('log')
for ax in axs[4:6]:
    ax.set_xlabel("Vapor Mass Fraction (%)")
    ax.set_ylabel("PDF")
    # ax.set_xscale('log')
for ax in axs[6:]:
    ax.set_xlabel("Entropy (J/kg/K)")
    ax.set_ylabel("Temperature (K)")
    # ax.set_xscale('log')
# make tight layout with no hspace
plt.tight_layout()
plt.savefig("paper2_initial_condition_velocity_vs_vmf.png", format='png', dpi=200)






fig, axs = plt.subplots(4, 2, figsize=(10, 20))
axs = axs.flatten()

for index, (run, verbose_run_name, iteration) in enumerate(runs):
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

    axs[index].set_title(f"Disk-bound particles at initial condition ({run_name})")

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
    final_disk_particles['velocity'] = np.sqrt(final_disk_particles['vx'] ** 2 + final_disk_particles['vy'] ** 2 + final_disk_particles['vz'] ** 2)
    df2 = final_disk_particles[final_disk_particles['tag'] % 2 == 0]

    vmf_final_disk = calc_vapor_mass_fraction_without_circularization_from_formatted(
        df2, phase_path, restrict_df=False
    ) * 100

    # draw a pdf of velocity
    axs[index].hist(df2['velocity'] / 1000, bins=100, density=True)
    axs[index].axvline(mean(df2['velocity'] / 1000), color='black', linestyle='--', label=f"Mean velocity: {mean(df2['velocity'] / 1000):.2f} km/s")

    # draw a pdf of entropy
    axs[index + 2].hist(df2['entropy'], bins=100, density=True)
    axs[index + 2].axvline(mean(df2['entropy']), color='black', linestyle='--', label=f"Mean entropy: {mean(df2['entropy']):.2f} J/kg/K")

    # draw a pdf of temperature
    axs[index + 4].hist(df2['temperature'], bins=100, density=True)
    axs[index + 4].axvline(mean(df2['temperature']), color='black', linestyle='--', label=f"Mean temperature: {mean(df2['temperature']):.2f} K")

    # draw a pdf of the vmf
    axs[index + 6].hist(df2['vmf_wo_circ'] * 100, bins=100, density=True)
    axs[index + 6].axvline(df2['vmf_wo_circ'].sum() / len(df2) * 100, color='black', linestyle='--', label=f"Mean VMF: {mean(df2['vmf_wo_circ'] * 100):.2f} %")

for ax in axs:
    ax.grid()
    # ax.legend()
for ax in axs[:2]:
    ax.set_xlabel("Velocity (km/s)")
    ax.set_ylabel("Probability Density")
for ax in axs[2:4]:
    ax.set_xlabel("Entropy (J/kg/K)")
    ax.set_ylabel("Probability Density")
for ax in axs[4:6]:
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel("Probability Density")
for ax in axs[6:]:
    ax.set_xlabel("VMF (%)")
    ax.set_ylabel("Probability Density")

# make tight layout with no hspace
plt.tight_layout()
plt.savefig("paper2_initial_condition_thermodynamics.png", format='png', dpi=200)
