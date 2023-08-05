import sys
from src.hugoniot import Hugoniot
import numpy as np
import pandas as pd
from math import log10
import matplotlib.pyplot as plt
from sys import exit
from ast import literal_eval
import string
from src.interpolation import GenericTrilinearInterpolation

plt.rcParams.update({'font.size': 14, })
plt.style.use('seaborn-colorblind')

angle = 'b073'
cutoff_densities = [5, 500, 1000, 2000]
# base_max_val_folders_loc = "C:/Users/Scott/OneDrive/Desktop/figures-selected/"
# base_max_val_folders_loc = "C:/Users/Scott/Desktop/"
base_max_val_folders_loc = "/Users/scotthull/Desktop/maxvals_8_5_2023/"

stewart_aneos_path = "src/phase_data/forst_STS.table.txt"
stewart_rho_u = "src/phase_data/forst_STS.rho_u.txt"
stewart_aneos_hugoniot_path = "src/phase_data/forstSTS__hugoniot.txt"
n_sph_aneos_path = "src/phase_data/duniteN.table.txt"
n_sph_aneos_hugoniot_path = "src/phase_data/duniteN__hugoniot.txt"
n_sph_rho_u = "src/phase_data/duniteN.rho_u.txt"

stewart_avg_surface_rho = 3729.72  # average surface density of the Stewart impactors acros cutoff densities
stewart_avg_surface_u = 874952.14  # average surface internal energy of the Stewart impactors acros cutoff densities
stewart_avg_surface_p = 548451811.03  # average surface pressure of the Stewart impactors acros cutoff densities

n_sph_avg_surface_rho = 3609.38  # average surface density of the N-SPH impactors across cutoff densities
n_sph_avg_surface_u = 1035922.78  # average surface internal energy of the N-SPH impactors acros cutoff densities
n_sph_avg_surface_p = 609283928.69  # average surface pressure of the N-SPH impactors across cutoff densities

def get_all_sims(angle, runs, high=True):
    fformat = "{}_{}_{}"
    tformat = "{}{}{}"
    names = []
    titles = []
    n = "S"
    if runs == "old":
        n = "N"
    if high and runs == "new" and angle == "b073":
        output_name = fformat.format(5, angle, "new") + "_high"
        names.append(output_name)
        title_name = tformat.format(5, angle, "S") + "-high"
        titles.append(title_name)
    for cd in cutoff_densities:
        output_name = fformat.format(cd, angle, runs)
        title_name = tformat.format(cd, angle, n)
        titles.append(title_name)
        names.append(output_name)
    if high and runs == "old" and angle == "b075":
        output_name = fformat.format(2000, angle, "old") + "_low"
        names.append(output_name)
        title_name = tformat.format(2000, angle, "N") + "-low"
        titles.append(title_name)
    return names, titles


def get_max_vals(angle, title, path):
    path = path + "/{}_{}_max_vals.csv".format(title, angle)
    print(path)
    with open(path, 'r') as infile:
        data = infile.readline().rstrip()
        return literal_eval(data)


colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
# fig, axs = plt.subplots(2, 2, figsize=(12, 12), sharex='all', sharey='all', gridspec_kw={"hspace": 0.14, "wspace": 0.0})
fig, axs = plt.subplots(2, 2, figsize=(12, 12), sharex='all', sharey='all')
axs = list(axs.flatten())
plotting_index = 0
for i in ["Stewart M-ANEOS", "N-SPH M-ANEOS"]:
    high = True
    # high = False
    # if angle == "b073":
    #     high = True
    runs = "new"
    if i == "N-SPH M-ANEOS":
        runs = "old"
    sims, titles = get_all_sims(angle, runs, high)

    initial_impact_path = base_max_val_folders_loc + "{}_maxvals_at_time_of_primary_impact".format(angle)
    all_time_maxval_path = base_max_val_folders_loc + "{}_maxvals_at_endstate".format(angle)
    title_appendix = "n"
    experimental = stewart_aneos_path
    aneos = stewart_aneos_hugoniot_path
    if i == "N-SPH M-ANEOS":
        title_appendix = "o"
        experimental = n_sph_aneos_path
        aneos = n_sph_aneos_hugoniot_path
    h = Hugoniot()
    h.read_ANEOS(experimental)
    rho1, P1, C1, U1, S1 = h.initial_conditions_aneos(aneos)
    rho1 = stewart_avg_surface_rho
    P1 = stewart_avg_surface_p
    U1 = stewart_avg_surface_u
    if "N-SPH" in i:
        rho1 = n_sph_avg_surface_rho
        P1 = n_sph_avg_surface_p
        U1 = n_sph_avg_surface_u
    T_s, Rho_s, Us_s, Up_s, P_s, U_s, S_s = h.rankine_hugoniot_equations(rho1, P1, U1)
    P_s = np.array(P_s) / 10 ** 9
    Us_s = np.array(Us_s) / 10 ** 3
    Up_s = np.array(Up_s) / 10 ** 3
    U_s = np.array(U_s) / 10 ** 6

    for j in [initial_impact_path, all_time_maxval_path]:
        maxvals = {}
        for s, t in zip(sims, titles):
            maxvals.update({t: get_max_vals(angle, t, j)})
        for j, (title, max_vals) in enumerate(maxvals.items()):
            x, y = [max_vals[p]['pressure'] / 10 ** 9 for p in max_vals.keys() if max_vals[p]['tag'] % 2 == 0], [max_vals[p]['entropy'] for p in max_vals.keys() if max_vals[p]['tag'] % 2 == 0]
            cd = int(title.split("b")[0])
            color_index = cutoff_densities.index(cd)
            if "high" in title or "low" in title:
                color_index = len(cutoff_densities)
            axs[plotting_index].scatter(x, y, marker='o', s=4, color=colors[color_index])
            # min_rho_particles = [
            #     i for i in max_vals.keys() if max_vals[i]['density'] == cd
            # ]
            # x, y = [max_vals[p]['pressure'] / 10 ** 9 for p in min_rho_particles], [max_vals[p]['entropy'] for p in min_rho_particles]
            # axs[plotting_index].scatter(x, y, marker='o', s=4, color='red')
        if plotting_index < 3:
            plotting_index += 2

    if "Stewart" in i:
        axs[0].plot(P_s, S_s, linewidth=4.0, color='black', label="Calculated")
        axs[0].plot(h.P_h / 10 ** 9, h.S_h, linewidth=4.0, color='black', linestyle="--", label="M-ANEOS")
        axs[2].plot(P_s, S_s, linewidth=4.0, color='black', label="Calculated")
        axs[2].plot(h.P_h / 10 ** 9, h.S_h, linewidth=4.0, color='black', linestyle="--", label="M-ANEOS")
    else:
        axs[1].plot(P_s, S_s, linewidth=4.0, color='black', label="Calculated")
        axs[1].plot(h.P_h / 10 ** 9, h.S_h, linewidth=4.0, color='black', linestyle="--", label="M-ANEOS")
        axs[3].plot(P_s, S_s, linewidth=4.0, color='black', label="Calculated")
        axs[3].plot(h.P_h / 10 ** 9, h.S_h, linewidth=4.0, color='black', linestyle="--", label="M-ANEOS")
    plotting_index = 1

if angle == "b073":
    r = ["5b073S-high"] + [r"$\rho_c$ = {} kg/m$^3$".format(cd) for cd in cutoff_densities]
elif angle == "b075":
    r = ["2000b075N-low"] + [r"$\rho_c$ = {} kg/m$^3$".format(cd) for cd in cutoff_densities]
for index, i in enumerate(r):
    if i == "5b073S-high" and angle == "b075":
        pass
    else:
        if "high" in i or "low" in i:
            index = len(cutoff_densities)
        else:
            index -= 1
        axs[0].scatter([], [], marker='o', s=4, label=i, color=colors[index])

axs[0].set_title("Stewart M-ANEOS: Primary Impact")
axs[2].set_title("Stewart M-ANEOS: Final Disk State")
axs[1].set_title("N-SPH M-ANEOS: Primary Impact")
axs[3].set_title("N-SPH M-ANEOS: Final Disk State")
axs[2].set_xlabel("Pressure (GPa)")
axs[3].set_xlabel("Pressure (GPa)")
axs[0].set_ylabel("Entropy (J/kg/K)")
axs[2].set_ylabel("Entropy (J/kg/K)")
letters = list(string.ascii_lowercase)
for index, ax in enumerate(axs):
    ax.grid(alpha=0.6)
    ax.set_xscale("log")
    ax.set_xlim(0, 10 ** 4)
    ax.set_ylim(0, 10000)
    # x1, x2, y1, y2 = ax.axis()
    # # x_loc = x2 - (2 * (x2 - x1))
    # y_loc = y2 - (0.08 * (y2 - y1))
    # ax.text(10 ** 3.5, y_loc, letters[index], fontweight="bold")
    # annotate the letter in the upper right-hand corner
    ax.annotate(letters[index], xy=(0.95, 0.95), xycoords='axes fraction', fontsize=18, fontweight='bold')
legend = axs[0].legend(loc='upper left')
for handle in legend.legendHandles:
    try:
        handle.set_sizes([30.0])
    except:
        pass

plt.tight_layout()
# plt.show()
plt.savefig("{}_hugoniot_with_max_pressure_vals.png".format(angle), format='png', dpi=200)




# for i in ["Stewart M-ANEOS", "N-SPH M-ANEOS"]:
#     high = False
#     if angle == "b073":
#         high = True
#     runs = "new"
#     if i == "N-SPH M-ANEOS":
#         runs = "old"
#     sims, titles = get_all_sims(angle, runs, high)
#
#     maxvals = {}
#     for s, t in zip(sims, titles):
#         maxvals.update({t: get_max_vals(angle, t)})
#     title_appendix = "n"
#     experimental = stewart_aneos_path
#     aneos = stewart_aneos_hugoniot_path
#     if i == "N-SPH M-ANEOS":
#         title_appendix = "o"
#         experimental = n_sph_aneos_path
#         aneos = n_sph_aneos_hugoniot_path
#     h = Hugoniot()
#     h.read_ANEOS(experimental)
#     rho1, P1, C1, U1, S1 = h.initial_conditions_aneos(aneos)
#     rho1 = stewart_avg_surface_rho
#     P1 = stewart_avg_surface_p
#     U1 = stewart_avg_surface_u
#     if "N-SPH" in i:
#         rho1 = n_sph_avg_surface_rho
#         P1 = n_sph_avg_surface_p
#         U1 = n_sph_avg_surface_u
#     T_s, Rho_s, Us_s, Up_s, P_s, U_s, S_s = h.rankine_hugoniot_equations(rho1, P1, U1)
#     P_s = np.array(P_s) / 10 ** 9
#     Us_s = np.array(Us_s) / 10 ** 3
#     Up_s = np.array(Up_s) / 10 ** 3
#     U_s = np.array(U_s) / 10 ** 6
#
#     for j, (title, max_vals) in enumerate(maxvals.items()):
#         fig, axs = plt.subplots(1, 1, figsize=(12, 12), gridspec_kw={"hspace": 0.14, "wspace": 0.18})
#         x, y = [max_vals[p]['pressure'] / 10 ** 9 for p in max_vals.keys()], [max_vals[p]['entropy'] for p in max_vals.keys()]
#         axs.scatter(x, y, marker='o', s=4, color='pink', label=title)
#
#         axs.plot(P_s, S_s, linewidth=4.0, label="Calculated")
#         axs.plot(h.P_h / 10 ** 9, h.S_h, linewidth=4.0, label="ANEOS")
#         axs.grid()
#         axs.set_xlabel("Pressure (GPa)")
#         axs.set_ylabel("Entropy (J/K)")
#         axs.set_xscale("log")
#         axs.legend()
#         axs.set_xlim(0, 10 ** 5)
#         axs.set_ylim(1000, 10000)
#
#         fig.suptitle(i + " " + title)
#
#         plt.savefig("{}_{}_hugoniot_with_max_vals.png".format(i, title), format='png', dpi=200)
#
#     # fig, axs = plt.subplots(1, 1, figsize=(12, 12), gridspec_kw={"hspace": 0.14, "wspace": 0.18})
#     # for j, (title, max_vals) in enumerate(maxvals.items()):
#     #     x, y = [max_vals[p]['time'] for p in max_vals.keys()], [max_vals[p]['pressure'] / 10 ** 9 for p in
#     #                                                                           max_vals.keys()]
#     #     axs.scatter(x, y, marker='o', s=10, label=title)
#     #
#     # axs.set_xlabel("Time (hrs)")
#     # axs.set_ylabel("Pressure (GPa)")
#     # axs.grid()
#     # axs.legend()
#     # axs.set_yscale("log")
#     # fig.suptitle(i)
#     # plt.savefig("{}_hugoniot_times_of_max_vals.png".format(i), format='png', dpi=200)

