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
    path = base_max_val_folders_loc + "{}_maxvals_all_disk_particles/{}_{}_max_vals.csv".format(angle, title, angle)
    with open(path, 'r') as infile:
        data = infile.readline().rstrip()
        return literal_eval(data)


high = False
# if angle == "b073":
#     high = True
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

    fig, axs = plt.subplots(1, 1, figsize=(12, 12), gridspec_kw={"hspace": 0.14, "wspace": 0.18})

    for j, (title, max_vals) in enumerate(maxvals.items()):
        x, y = [max_vals[p]['pressure'] / 10 ** 9 for p in max_vals.keys()], [max_vals[p]['entropy'] for p in max_vals.keys()]
        axs.scatter(x, y, marker='o', s=4, label=title)

    axs.plot(P_s, S_s, linewidth=2.0, label="Calculated")
    axs.plot(h.P_h / 10 ** 9, h.S_h, linewidth=2.0, label="ANEOS")
    axs.grid()
    axs.set_xlabel("Pressure (GPa)")
    axs.set_ylabel("Entropy (J/K)")
    axs.set_xscale("log")
    axs.legend()

    fig.suptitle(i)

    plt.savefig("{}_hugoniot_with_max_vals.png".format(i), format='png', dpi=200)

