#!/usr/bin/env python3
import os
import csv
import string
from random import randint
from statistics import mean
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LogNorm
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from src.combine import CombineFile
from src.animate import animate
from src.vapor import calc_vapor_mass_fraction_without_circularization_from_formatted

# use colorblind-friendly colors from seaborn
plt.style.use('seaborn-colorblind')
# increase font size
plt.rcParams.update({'font.size': 18})

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

# fig, axs = plt.subplots(4, 2, figsize=(10, 20))
# axs = axs.flatten()
#
# for index, (run, verbose_run_name, iteration) in enumerate(runs):
#     # if len(run) == 0, then skip this part of the loop
#     if len(run) == 0:
#         continue
#     run_name = run.split('/')[-1]
#     base_path = run.split(f"/{run_name}")[0] + "/"
#     run = run_name
#     if not os.path.exists(f"{run}_vel_profile"):
#         os.mkdir(f"{run}_vel_profile")
#     circ_path = base_path + f"/{run}/circularized_{run}"
#     end_time_df = pd.read_csv(
#         os.path.join(circ_path, f"{max_iteration}.csv"),
#     )
#     end_time_disk = end_time_df[end_time_df["label"] == "DISK"]
#     end_time_particle_ids = end_time_disk["id"].values
#     time, iterations, mean_disk_vel, mean_disk_entropy, mean_disk_temperature, final_disk_vmf, all_silicate_vmf = [], [], [], [], [], [], []
#     phase_path = new_phase_path if "new" in run else old_phase_path
#
#
#     path = base_path + f"{run}/{run}"
#     to_fname = "merged_{}_{}.dat".format(iteration, randint(0, 100000))
#     cf = CombineFile(num_processes=number_processes, time=iteration, output_path=path, to_fname=to_fname)
#     c = cf.combine()
#     formatted_time = round(cf.sim_time * 0.000277778, 2)
#
#     combined_file = pd.read_csv(to_fname, skiprows=2, header=None, delimiter="\t")
#
#     os.remove(to_fname)
#
#     final_disk_particles = combined_file[combined_file[0].isin(end_time_particle_ids)]
#
#     id, x, y, z, vx, vy, vz, entropy, temperature = combined_file[0], combined_file[3], combined_file[4], \
#         combined_file[5], combined_file[6], combined_file[7], combined_file[8], combined_file[13], combined_file[14]
#
#     id_disk, x_disk, y_disk, z_disk, vx_disk, vy_disk, vz_disk, entropy_disk, temperature_disk = \
#         final_disk_particles[0], final_disk_particles[3], final_disk_particles[4], final_disk_particles[5], \
#             final_disk_particles[6], final_disk_particles[7], final_disk_particles[8], final_disk_particles[13], \
#             final_disk_particles[14]
#
#     velocity = np.sqrt(vx_disk ** 2 + vy_disk ** 2 + vz_disk ** 2)
#
#     time.append(formatted_time)
#     mean_disk_vel.append(velocity.mean())
#     mean_disk_entropy.append(entropy_disk.mean())
#     mean_disk_temperature.append(temperature_disk.mean())
#
#     # replace final_disk_particles headers with the correct headers
#     final_disk_particles.columns = [
#         'id', 'tag', 'mass', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'density', 'internal_energy', 'pressure',
#         'potential_energy', 'entropy', 'temperature'
#     ]
#     combined_file.columns = final_disk_particles.columns
#     final_disk_particles['velocity'] = np.sqrt(final_disk_particles['vx'] ** 2 + final_disk_particles['vy'] ** 2 + final_disk_particles['vz'] ** 2)
#     df2 = final_disk_particles[final_disk_particles['tag'] % 2 == 0]
#
#     vmf_final_disk = calc_vapor_mass_fraction_without_circularization_from_formatted(
#         df2, phase_path, restrict_df=False
#     ) * 100
#
#     axs[index].scatter(
#         df2['velocity'] / 1000, df2['vmf_wo_circ'] * 100, s=5, label=verbose_run_name
#     )
#     axs[index].axvline(mean(df2['velocity'] / 1000), color='black', linestyle='--', linewidth=2.0, label=f"Mean velocity: {mean(df2['velocity'] / 1000):.2f} km/s")
#     axs[index].axhline(df2['vmf_wo_circ'].sum() / len(df2) * 100, color='black', linestyle='--', linewidth=2.0, label=f"Mean VMF: {mean(df2['vmf_wo_circ'] * 100):.2f} %")
#     axs[index].set_title(f"Disk-bound particles at initial condition ({run_name})")
#
#     # on the bottom row, plot a CDF of the VMFs
#     sorted_vmf = df2['vmf_wo_circ'].sort_values()
#     cdf = sorted_vmf.rank(method='average', pct=True)
#     axs[index + 2].plot(sorted_vmf * 100, cdf, linewidth=2.0)
#     axs[index + 2].axvline(df2['vmf_wo_circ'].sum() / len(df2) * 100, color='black', linestyle='--', linewidth=2.0,
#                        label=f"Mean VMF: {mean(df2['vmf_wo_circ'] * 100):.2f} %")
#
#     # plot a PDF of the VMFs
#     axs[index + 4].hist(sorted_vmf * 100, bins=100, density=False)
#     axs[index + 4].axvline(df2['vmf_wo_circ'].sum() / len(df2) * 100, color='black', linestyle='--', linewidth=2.0,
#                        label=f"Mean VMF: {mean(df2['vmf_wo_circ'] * 100):.2f} %")
#     axs[index + 4].text(
#         0.05, 0.90, f"vmf@0%: {len(df2[df2['vmf_wo_circ'] == 0])}\n"
#                     f"0<vmf<1%: {len(df2[(df2['vmf_wo_circ'] > 0) & (df2['vmf_wo_circ'] < .01)])}\n"
#                     f"1<vmf<99.99%: {len(df2[(df2['vmf_wo_circ'] > .01) & (df2['vmf_wo_circ'] < 0.9999)])}\n"
#                     f"vmf@100%: {len(df2[df2['vmf_wo_circ'] > 0.9999])}",
#         transform=axs[index + 4].transAxes, verticalalignment='top'
#     )
#     axs[index + 6].scatter(
#         df2['entropy'], df2['temperature'], s=5, marker="."
#     )
#
# for ax in axs:
#     ax.grid()
#     ax.legend()
# for ax in axs[:2]:
#     ax.set_xlabel("Velocity (km/s)")
#     ax.set_ylabel("Vapor Mass Fraction (%)")
#     # ax.set_yscale('log')
# for ax in axs[2:4]:
#     ax.set_xlabel("Vapor Mass Fraction (%)")
#     ax.set_ylabel("CDF")
#     # ax.set_xscale('log')
# for ax in axs[4:6]:
#     ax.set_xlabel("Vapor Mass Fraction (%)")
#     ax.set_ylabel("PDF")
#     # ax.set_xscale('log')
# for ax in axs[6:]:
#     ax.set_xlabel("Entropy (J/kg/K)")
#     ax.set_ylabel("Temperature (K)")
#     # ax.set_xscale('log')
# # make tight layout with no hspace
# plt.tight_layout()
# plt.savefig("paper2_initial_condition_velocity_vs_vmf.png", format='png', dpi=200)






fig, axs = plt.subplots(4, 2, figsize=(18, 25))
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

    axs[index].set_title(verbose_run_name, fontsize=20)

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

    # output df2
    df2.to_csv(f"{run}_df2.csv", index=False)

    df2_intermediate_vmf = df2[df2['vmf_wo_circ'] > 0.01]
    df2_intermediate_vmf = df2_intermediate_vmf[df2_intermediate_vmf['vmf_wo_circ'] < 1]


    # draw a pdf of velocity
    axs[index].hist(df2['velocity'] / 1000, bins=100, density=False)
    axs[index].axvline(mean(df2['velocity'] / 1000), color='black', linestyle='--', linewidth=2.0, label=f"Mean velocity: {mean(df2['velocity'] / 1000):.2f} km/s")
    # on the right axis, plot the CDF of velocity
    sorted_vel = df2['velocity'].sort_values()
    cdf = sorted_vel.rank(method='average', pct=True)
    axs2 = axs[index].twinx()
    axs2.plot(sorted_vel / 1000, cdf, linewidth=2.0, color='red')
    axs2.set_ylabel("CDF", fontsize=16)
    axs2.tick_params(axis='both', which='major', labelsize=16)
    axs[index].text(
        0.05, 0.9, "Avg. Velocity: {:.2f} km/s".format(mean(df2['velocity'] / 1000)),
            transform=axs[index].transAxes, verticalalignment='top', horizontalalignment='left', fontsize=18
    )

    # draw a pdf of entropy
    axs[index + 2].hist(df2['entropy'] / 1000, bins=100, density=False)
    axs[index + 2].axvline(mean(df2['entropy']) / 1000, color='black', linestyle='--', linewidth=2.0, label=f"Mean entropy: {mean(df2['entropy']):.2f} J/kg/K")
    # on the right axis, plot the CDF of entropy
    sorted_entropy = df2['entropy'].sort_values()
    cdf = sorted_entropy.rank(method='average', pct=True)
    axs2 = axs[index + 2].twinx()
    axs2.plot(sorted_entropy / 1000, cdf, linewidth=2.0, color='red')
    axs2.set_ylabel("CDF", fontsize=16)
    axs2.tick_params(axis='both', which='major', labelsize=16)
    # get the PDF of the df2_intermediate_vmf and calculate the x location of the largest peak
    hist_values, bin_edges = np.histogram(df2_intermediate_vmf['entropy'], bins=100)
    largest_peak_index = np.argmax(hist_values)
    largest_peak_x = bin_edges[largest_peak_index + 1]  # Adding 1 to get the upper edge of the bin
    axs[index + 2].axvline(largest_peak_x / 1000, color='red', linestyle='--', linewidth=2.0, label=f"Partially vaporized: {largest_peak_x:.2f} J/kg/K")
    axs[index + 2].text(
        0.25, 0.8, "Avg. Entropy (bulk): {:.2f} J/kg/K\nAvg. Entropy (intermediate): {:.2f} J/kg/K".format(mean(df2['entropy'] / 1000), mean(df2_intermediate_vmf['entropy'] / 1000)),
            transform=axs[index + 2].transAxes, verticalalignment='top', horizontalalignment='left', fontsize=18
    )

    # draw a pdf of temperature
    axs[index + 4].hist(df2['temperature'] / 1000, bins=100, density=False)
    axs[index + 4].axvline(mean(df2['temperature']) / 1000, color='black', linestyle='--', linewidth=2.0, label=f"Mean temperature: {mean(df2['temperature']):.2f} K")
    # on the right axis, plot the CDF of temperature
    sorted_temperature = df2['temperature'].sort_values()
    cdf = sorted_temperature.rank(method='average', pct=True)
    axs2 = axs[index + 4].twinx()
    axs2.plot(sorted_temperature / 1000, cdf, linewidth=2.0, color='red')
    axs2.set_ylabel("CDF", fontsize=16)
    axs2.tick_params(axis='both', which='major', labelsize=16)
    # get the PDF of the df2_intermediate_vmf and calculate the x location of the largest peak
    hist_values, bin_edges = np.histogram(df2_intermediate_vmf['temperature'], bins=100)
    largest_peak_index = np.argmax(hist_values)
    largest_peak_x = bin_edges[largest_peak_index + 1]  # Adding 1 to get the upper edge of the bin
    axs[index + 4].axvline(largest_peak_x / 1000, color='red', linestyle='--', linewidth=2.0, label=f"Partially vaporized: {largest_peak_x:.2f} K")
    axs[index + 4].text(
        0.25, 0.8, "Avg. Temperature (bulk): {:.2f} km/s\nAvg. Temperature (intermediate): {:.2f} km/s".format(mean(df2['temperature'] / 1000), mean(df2_intermediate_vmf['temperature'] / 1000)),
            transform=axs[index + 4].transAxes, verticalalignment='top', horizontalalignment='left', fontsize=18
    )

    # draw a pdf of the vmf
    axs[index + 6].hist(df2['vmf_wo_circ'] * 100, bins=100, density=False)
    axs[index + 6].axvline(df2['vmf_wo_circ'].sum() / len(df2) * 100, color='black', linestyle='--', linewidth=2.0, label=f"Mean VMF: {mean(df2['vmf_wo_circ'] * 100):.2f} %")
    # on the right axis, plot the CDF of vmf
    sorted_vmf = df2['vmf_wo_circ'].sort_values()
    cdf = sorted_vmf.rank(method='average', pct=True)
    axs2 = axs[index + 6].twinx()
    axs2.plot(sorted_vmf * 100, cdf, linewidth=2.0, color='red')
    axs2.set_ylabel("CDF", fontsize=16)
    axs2.tick_params(axis='both', which='major', labelsize=16)
    # get the PDF of the df2_intermediate_vmf and calculate the x location of the largest peak
    hist_values, bin_edges = np.histogram(df2_intermediate_vmf['vmf_wo_circ'] * 100, bins=100)
    largest_peak_index = np.argmax(hist_values)
    largest_peak_x = bin_edges[largest_peak_index + 1]  # Adding 1 to get the upper edge of the bin
    axs[index + 6].axvline(largest_peak_x, color='red', linestyle='--', linewidth=2.0, label=f"Partially vaporized: {largest_peak_x:.2f} %")
    axs[index + 6].text(
        0.25, 0.8, "Avg. VMF (bulk): {:.2f} %\nAvg. VMF (intermediate): {:.2f} %".format(mean(df2['vmf_wo_circ'] * 100), mean(df2_intermediate_vmf['vmf_wo_circ'] * 100)),
            transform=axs[index + 6].transAxes, verticalalignment='top', horizontalalignment='left', fontsize=18
    )

    # axs[index + 8].scatter(
    #     df2['temperature'] / 1000, df2['vmf_wo_circ'] * 100, s=10, marker="."
    # )
    # axs[index + 8].axhline(df2['vmf_wo_circ'].sum() / len(df2) * 100, color='black', linestyle='--')
    # axs[index + 8].axvline(mean(df2['temperature'] / 1000), color='black', linestyle='--')
    # # axs[index + 8].text(
    # #     0.50, 0.8, "Avg. Temperature: {:.2f} K\nAvg. VMF: {:.2f} %".format(mean(df2['temperature']), np.sum(df2['vmf_wo_circ']) / len(df2) * 100),
    # #         transform=axs[index + 8].transAxes, verticalalignment='top', horizontalalignment='left', fontsize=16
    # #     # )



for ax in axs:
    ax.grid()
    # increase axis font size
    ax.tick_params(axis='both', which='major', labelsize=16)
    # turn on minor ticks
    ax.minorticks_on()
    ax.set_yscale('log')
    # ax.legend()
for ax in axs[:2]:
    ax.set_xlabel("Velocity (km/s)", fontsize=16)
    ax.set_ylabel("# Particles", fontsize=16)
for ax in axs[2:4]:
    ax.set_xlabel("Entropy (1000 J/kg/K)", fontsize=16)
    ax.set_ylabel("# Particles", fontsize=16)
for ax in axs[4:6]:
    ax.set_xlabel("Temperature (1000 K)", fontsize=16)
    ax.set_ylabel("# Particles", fontsize=16)
for ax in axs[6:8]:
    ax.set_xlabel("VMF (%)", fontsize=16)
    ax.set_ylabel("# Particles", fontsize=16)
# for ax in axs[8:10]:
#     ax.set_xlabel("Temperature (1000 K)", fontsize=16)
#     ax.set_ylabel("VMF (%)", fontsize=16)
#     ax.set_xscale("log")
#     # ax.legend()

# make tight layout with no hspace
plt.tight_layout()
plt.savefig("paper2_initial_condition_thermodynamics.png", format='png', dpi=200)

def center_of_mass(df: pd.DataFrame):
    """
    Calculate the center of mass of the particles dataframe
    """
    com_x = np.sum(df['x'] * df['mass']) / np.sum(df['mass'])
    com_y = np.sum(df['y'] * df['mass']) / np.sum(df['mass'])
    com_z = np.sum(df['z'] * df['mass']) / np.sum(df['mass'])
    return com_x, com_y, com_z


# make a 2 row 3 column figure
fig, axs = plt.subplots(2, 3, figsize=(15, 10), sharex='all', sharey='all', gridspec_kw={'hspace': 0, 'wspace': 0})
axs = axs.flatten()

temperature_normalizer = Normalize(1800 / 1000, 6000 / 1000)
entropy_normalizer = Normalize(2000 / 1000, 6000 / 1000)
vmf_normalizer = Normalize(0, 10)
cmap = cm.get_cmap('cool')

plot_index = 0
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

    com_x, com_y, com_z = center_of_mass(combined_file[combined_file['tag'] == 1])

    vmf_final_disk = calc_vapor_mass_fraction_without_circularization_from_formatted(
        df2, phase_path, restrict_df=False
    ) * 100

    axs[plot_index].scatter(
        (combined_file['x'] - com_x) / 1000 / 10 ** 3, (combined_file['y'] - com_y / 10 ** 3) / 1000 / 10 ** 3, s=5, marker=".", color='black'
    )
    axs[plot_index].scatter(
        (df2['x'] - com_x) / 1000 / 10 ** 3, (df2['y'] - com_y) / 1000 / 10 ** 3, s=5, marker=".", color=cmap(temperature_normalizer(df2['temperature'] / 1000))
    )
    axs[plot_index].text(
        0.10, 0.10, f"{verbose_run_name}", transform=axs[plot_index].transAxes, verticalalignment='top',
        fontsize=16
    )
    plot_index += 1
    axs[plot_index].scatter(
        (combined_file['x'] - com_x) / 1000 / 10 ** 3, (combined_file['y'] - com_y) / 1000 / 10 ** 3, s=5, marker=".", color='black'
    )
    axs[plot_index].scatter(
        (df2['x'] - com_x) / 1000 / 10 ** 3, (df2['y'] - com_y) / 1000 / 10 ** 3, s=5, marker=".", color=cmap(temperature_normalizer(df2['entropy'] / 1000))
    )
    plot_index += 1
    axs[plot_index].scatter(
        (combined_file['x'] - com_x) / 1000 / 10 ** 3, (combined_file['y'] - com_y) / 1000 / 10 ** 3, s=5, marker=".", color='black'
    )
    axs[plot_index].scatter(
        (df2['x'] - com_x) / 1000 / 10 ** 3, (df2['y'] - com_y) / 1000 / 10 ** 3, s=5, marker=".", color=cmap(vmf_normalizer(df2['vmf_wo_circ'] * 100))
    )
    plot_index += 1

for ax in [axs[0], axs[3]]:
    ax.set_ylabel(r"y ($10^3$ km)", fontsize=16)
for ax in axs[-3:]:
    ax.set_xlabel(r"x ($10^3$ km)", fontsize=16)
for ax, label, normalizer in zip(axs[:3], ['Temperature (1000 K)', 'Entropy (1000 J/kg/K)', 'VMF (%)'],
                                 [temperature_normalizer, entropy_normalizer, vmf_normalizer]):
    sm = cm.ScalarMappable(norm=normalizer, cmap=cmap)
    sm.set_array([])
    cbaxes = inset_axes(ax, width="50%", height="8%", loc=1, borderpad=1.8)
    cbar = plt.colorbar(sm, cax=cbaxes, orientation='horizontal')
    cbar.ax.tick_params(labelsize=14)
    cbar.ax.set_title(label, fontsize=14)

for ax in axs:
    # increase font size
    ax.tick_params(axis='both', which='major', labelsize=14)


plt.savefig("paper2_initial_condition_thermodynamics_as_func_space.png", format='png', dpi=200)
