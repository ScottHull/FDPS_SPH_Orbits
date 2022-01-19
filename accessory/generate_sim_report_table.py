import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# runs = ["5_b073_new", "5_b073_old", "500_b073_new", "500_b073_old", "1000_b073_new", "1000_b073_old",
#         "2000_b073_new", "2000_b073_old"]
runs = ["5_b075_new", "5_b075_old", "500_b075_new", "500_b075_old", "1000_b075_new", "1000_b075_old",
        "2000_b075_new", "2000_b075_old"]
table_rows = [
    "mass_protoplanet (M_E)", "mass_disk (M_L)", "mass_escaped (M_L)", "disk_mass_beyond_roche (M_L)",
    "disk_angular_momentum (L_EM)", "total_angular_momentum (L_EM)", "disk vmf", "iron_disk_mass_fraction",
    "iron_disk_mass_fraction_beyond_roche",
    "a", "b", "average_density", "num_particles_planet", "num_particles_disk", "num_particles_escaping",
    "num_particles_error"
]
row_names = [
    "{Planet Mass ($M_\oplus$)}", "{Disk Mass ($M_L$)}", "{Escaping Mass ($M_L$)}", "{$\%$ Disk Mass $\geq R_{Roche}$}",
    "{Disk AM ($L_{EM}$)}", "{Total AM ($L_{EM}$)}", "{Disk VMF  ($\%$)}", "{Disk Iron $\%$}",
    "{Disk Iron Mass $\geq R_{Roche}$ ($\%$)}", "{Planet $a$ (km)}", "{Planet $b$ (km)}", "{Planet Avg. Density}",
    "{$N_{planet}$}", "{$N_{disk}$}", "{$N_{escape}$}", "{$N_{error}$}"
]
percentages = ["disk vmf", "iron_disk_mass_fraction", "iron_disk_mass_fraction_beyond_roche"]
km = ["a", "b"]
whole_nums = ["num_particles_planet", "num_particles_disk", "num_particles_escaping", "num_particles_error"]

to_path = "/Users/scotthull/Documents - Scott’s MacBook Pro/PhD Research/Paper1/sim_endstate_reports"

dat = {}
fname = "table_formatted_endstates.txt"
if fname in os.listdir(os.getcwd()):
    os.remove(fname)
with open(fname, 'w') as outfile:
    header = "{Simulation Name}"
    for i in runs:
        header += " & {" + i.replace("_", "").replace("new", "n").replace("old", "o") + "}"
    header += " \\\ \midrule\n"
    outfile.write(header)
    for index, tr in enumerate(table_rows):
        # line = "{" + tr + "}"
        line = row_names[index]
        for r in runs:
            r_name = r.replace("_", "").replace("new", "n").replace("old", "o")
            try:
                df = pd.read_csv(to_path + "/{}.txt".format(r), header=None, index_col=0, delimiter="\t")
                d = float(df[1][tr])
                if tr in percentages:
                    d *= 100.0
                elif tr in km:
                    d /= 1000.0
                if tr not in whole_nums:
                    d = round(d, 2)
                else:
                    d = int(d)
                if row_names[index] not in dat.keys():
                    dat.update({row_names[index]: {}})
                dat[row_names[index]].update({r_name: d})
                line += " & {" + str(d) + "}"
            except Exception as e:
                print(e)
                line += " & {}"
        outfile.write(line + " \\\ \midrule\n")


def grouped_bc(dat):
    fig, axs = plt.subplots(4, 4, figsize=(20, 20), gridspec_kw={"hspace": 0.16, "wspace": 0.16})
    axs = axs.flatten()
    for ax in axs:
        ax.grid(alpha=0.4)
    index = 0
    for r in row_names:
        labels = ["5", "500", "1000", "2000"]
        new_dat = [dat[r][i] for i in dat[r].keys() if "n" in i]
        old_dat = [dat[r][i] for i in dat[r].keys() if "o" in i]
        x = np.arange(len(labels))  # the label locations
        width = 0.35  # the width of the bars
        rects1 = axs[index].bar(x - width / 2, new_dat, width, label='New EoS')
        rects2 = axs[index].bar(x + width / 2, old_dat, width, label='Old EoS')
        if index >= 16 - 4:
            axs[index].set_xlabel("Cutoff Density (kg/m3)")
        # axs[index].set_ylabel()
        axs[index].set_title(r.replace("{", "").replace("}", ""))
        axs[index].set_xticks(x)
        axs[index].set_xticklabels(labels)
        axs[index].legend()

        # axs[index].bar_label(rects1, padding=3)
        # axs[index].bar_label(rects2, padding=3)
        index += 1
    # plt.show()
    plt.savefig("gi_bars.png", format='png')


grouped_bc(dat)
