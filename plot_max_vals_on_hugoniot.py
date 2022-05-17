from src.hugoniot import Hugoniot
import numpy as np
import pandas as pd
from math import log10
import matplotlib.pyplot as plt
from sys import exit
from ast import literal_eval

plt.rcParams.update({'font.size': 12, })
plt.style.use('seaborn-colorblind')

angle = 'b073'
cutoff_densities = [5, 500, 1000, 2000]
base_max_val_folders_loc = "/Users/scotthull/Desktop/"

stewart_aneos_path = "src/phase_data/forst_STS.table.txt"
stewart_aneos_hugoniot_path = "src/phase_data/forstSTS__hugoniot.txt"
gadget_aneos_path = "src/phase_data/duniteN.table.txt"
gadget_aneos_hugoniot_path = "src/phase_data/duniteN__hugoniot.txt"


def get_all_sims(angle, high=True):
    fformat = "{}_{}_{}"
    tformat = "{}{}{}"
    names = []
    titles = []
    for runs in ["new", "old"]:
        n = "n"
        if runs == "old":
            n = "o"
        for cd in cutoff_densities:
            output_name = fformat.format(cd, angle, runs)
            title_name = tformat.format(cd, angle, n)
            titles.append(title_name)
            names.append(output_name)
    if high:
        output_name = fformat.format(5, angle, "new") + "_high"
        names.append(output_name)
        title_name = tformat.format(5, angle, "n") + "-high"
        titles.append(title_name)
    return names, titles


def get_max_vals(angle, title):
    path = base_max_val_folders_loc + "{}_maxvals/{}_{}_max_vals.csv".format(angle, title, angle)
    with open(path, 'r') as infile:
        data = infile.readline().rstrip()
        return literal_eval(data)


high = False
if angle == "b073":
    high = True
sims, titles = get_all_sims(angle, high)

maxvals = {}
for s, t in zip(sims, titles):
    maxvals.update({t: get_max_vals(angle, t)})


for i in ["Stewart M-ANEOS", "GADGET M-ANEOS"]:
    title_appendix = "n"
    experimental = stewart_aneos_path
    aneos = stewart_aneos_hugoniot_path
    if i == "GADGET M-ANEOS":
        title_appendix = "o"
        experimental = gadget_aneos_path
        aneos = gadget_aneos_hugoniot_path
    h = Hugoniot()
    h.read_ANEOS(experimental)
    rho1, P1, C1, U1, S1 = h.initial_conditions_aneos(aneos)
    T_s, Rho_s, Us_s, Up_s, P_s, U_s, S_s = h.rankine_hugoniot_equations(rho1, P1, U1)
    P_s = np.array(P_s) / 10 ** 9
    Us_s = np.array(Us_s) / 10 ** 3
    Up_s = np.array(Up_s) / 10 ** 3
    U_s = np.array(U_s) / 10 ** 6

    fig, axs = plt.subplots(3, 2, figsize=(12, 18), gridspec_kw={"hspace": 0.14, "wspace": 0.18})
    axs = axs.flatten()

    axs[0].plot(
        Up_s, Us_s, linewidth=2.0, label="Calculated"
    )
    axs[0].plot(
        h.Up_h / 1000, h.Us_h / 1000, linewidth=2.0, linestyle="--", label="ANEOS"
    )
    axs[0].set_xlabel("Up (km/s)")
    axs[0].set_ylabel("Us (km/s)")
    axs[1].plot(
        P_s, T_s, linewidth=2.0, label="Calculated"
    )
    axs[1].plot(
        h.P_h / 10 ** 9, h.T_h, linewidth=2.0, linestyle="--", label="ANEOS"
    )
    axs[1].set_xlabel("P (GPa)")
    axs[1].set_ylabel("T (K)")
    # axs[1].set_xlim(10 ** 0, 10 ** 4)
    axs[2].plot(
        P_s, Rho_s, linewidth=2.0, label="Calculated"
    )
    axs[2].plot(
        h.P_h / 10 ** 9, h.rho_h, linewidth=2.0, linestyle="--", label="ANEOS"
    )
    axs[2].set_xlabel("P (GPa)")
    axs[2].set_ylabel("Density (kg/m^3)")
    # axs[2].set_xlim(10 ** 0, 10 ** 4)
    axs[3].plot(
        P_s, U_s, linewidth=2.0, label="Calculated"
    )
    axs[3].plot(
        h.P_h / 10 ** 9, h.U_h / 10 ** 6, linewidth=2.0, linestyle="--", label="ANEOS"
    )
    axs[3].set_xlabel("P (GPa)")
    axs[3].set_ylabel("Internal Energy (J/K/kg)")
    # axs[3].set_xlim(10 ** 0, 10 ** 4)
    axs[4].plot(
        Rho_s, T_s, linewidth=2.0, label="Calculated"
    )
    axs[4].plot(
        h.rho_h, h.T_h, linewidth=2.0, linestyle="--", label="ANEOS"
    )
    axs[4].set_xlabel("Density (kg/m^3)")
    axs[4].set_ylabel("T (K)")

    axs[5].plot(
        Rho_s, S_s, linewidth=2.0, label="Calculated"
    )
    axs[5].plot(
        h.rho_h, h.S_h, linewidth=2.0, linestyle="--", label="ANEOS"
    )
    axs[5].set_xlabel("Density (kg/m^3)")
    axs[5].set_ylabel("Entropy (J/kg/K)")

    runs = [j for j in titles if "{}{}".format(angle, title_appendix) in j]
    for run in runs:
        axs[1].scatter(
            maxvals[run]['pressure']['value'], maxvals[run]['temperature']['value'],
            marker="o", label=run
        )
        axs[2].scatter(
            maxvals[run]['pressure']['value'], maxvals[run]['density']['value'],
            marker="o", label=run
        )
        axs[3].scatter(
            maxvals[run]['pressure']['value'], maxvals[run]['internal energy']['value'] / 10 ** 6,
            marker="o", label=run
        )
        axs[4].scatter(
            maxvals[run]['density']['value'], maxvals[run]['temperature']['value'],
        )
        axs[5].scatter(
            maxvals[run]['density']['value'], maxvals[run]['entropy']['value'],
            marker="o", label=run
        )

    axs[1].legend(loc="upper left")

    for ax in axs:
        ax.grid(alpha=0.4)
    for ax in axs[1:-2]:
        ax.set_xscale("log")
    for ax in axs[1:3]:
        ax.set_yscale("log")

    fig.suptitle(i)

    plt.savefig("{}_hugoniot_with_max_vals.png".format(i), format='png', dpi=200)

