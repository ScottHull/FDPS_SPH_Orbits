import pandas as pd
import numpy as np
import string
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 14, })
# plt.style.use("dark_background")
plt.style.use('seaborn-colorblind')

iteration = 1800
angles = ['b073', 'b075']
other = "pressure"
other_normalizer = 1e9
other_label = "Pressure (GPa)"
other_ylim = (None, None)
cutoff_densities = [5, 500, 1000, 2000]
base_path = "/home/theia/scotthull/Paper1_SPH/gi/"
new_phase_path = "src/phase_data/forstSTS__vapour_curve.txt"
old_phase_path = "src/phase_data/duniteN__vapour_curve.txt"

new_phase_df = pd.read_fwf(new_phase_path, skiprows=1,
                           names=["temperature", "density_sol_liq", "density_vap", "pressure",
                                  "entropy_sol_liq", "entropy_vap"])
old_phase_df = pd.read_fwf(old_phase_path, skiprows=1,
                           names=["temperature", "density_sol_liq", "density_vap", "pressure",
                                  "entropy_sol_liq", "entropy_vap"])


def get_all_sims(angle, high=True):
    fformat = "{}_{}_{}"
    tformat = "{}{}{}"
    names = []
    titles = []
    for runs in ["new", "old"]:
        high_res_name = None
        high_res_title = None
        n = "S"
        if runs == "old":
            n = "N"
        for cd in cutoff_densities:
            output_name = fformat.format(cd, angle, runs)
            title_name = tformat.format(cd, angle, n)
            titles.append(title_name)
            names.append(output_name)
            if cd == 5 and high and runs == "new" and angle == "b073":
                high_res_name = fformat.format(cd, angle, runs) + "_high"
                high_res_title = tformat.format(cd, angle, n) + "-high"
            if cd == 2000 and high and runs == "old" and angle == "b075":
                high_res_name = fformat.format(cd, angle, runs) + "_low"
                high_res_title = tformat.format(cd, angle, n) + "-low"
        if high_res_name is not None and high_res_title is not None:
            names.append(high_res_name)
            titles.append(high_res_title)
    return names, titles


fig, axs = plt.subplots(2, 2, figsize=(16, 9), sharex='all', sharey='all')
axs = axs.flatten()
b073n_index, b073o_index, b075n_index, b075o_index = 0, 1, 2, 3
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

critical_point_new = max(new_phase_df['temperature'])
critical_point_old = max(old_phase_df['temperature'])

def plot_phase_diagrams(other_name, other_norm):
    for angle in angles:
        high = True
        # if angle == "b075":
        #     high = False
        sims, titles = get_all_sims(angle=angle, high=high)
        for s, t in zip(sims, titles):
            marker = "."
            to_path = base_path + "{}/circularized_{}/{}.csv".format(s, s, iteration)
            df = pd.read_csv(to_path)
            disk = df[df['label'] == "DISK"]
            disk_filtered = disk[disk['circ_entropy_delta'] < 5000]
            other, density = disk_filtered[other_name], disk_filtered['density']
            # get min other that is not 0
            min_other = min(other[other != 0])
            # get the number of particles at min density
            min_other_count = len(disk_filtered[disk_filtered[other_name] == min_other])
            cd = int(s.split("_")[0])
            color = colors[cutoff_densities.index(cd)]
            if "high" in s:
                color = colors[cutoff_densities.index(cutoff_densities[-1]) + 1]
            if "low" in s:
                color = colors[cutoff_densities.index(cutoff_densities[-1]) + 1]
            if "old" in s:
                phase_curve = old_phase_df
            to_index = 0
            if "b073" in s and "new" in s:
                to_index = 0
            elif "b073" in s and "old" in s:
                to_index = 1
            elif "b075" in s and "new" in s:
                to_index = 2
            else:
                to_index = 3
            label = r"$\rho_c = {}$ kg/m$^3$".format(cd)
            if "high" in s:
                label = t
            axs[to_index].scatter(
                density, other / other_norm, s=1, marker=marker, alpha=1, color=color
            )
            # axs[to_index].scatter(
            #     min_other / other_norm, min_other_count, s=400, marker=marker, alpha=1, color=color
            # )

for cd in cutoff_densities:
    label = r"$\rho_c = {}$ kg/m$^3$".format(cd)
    color = colors[cutoff_densities.index(cd)]
    axs[0].scatter(
        [], [], marker=".", s=120, alpha=1, color=color, label=label
    )
axs[0].scatter(
        [], [], marker=".", s=120, alpha=1, color=colors[cutoff_densities.index(cutoff_densities[-1]) + 1], label="5b073S-high &\n2000b075N-low"
    )

axs[0].set_title("Stewart M-ANEOS ($b=0.73$)")
axs[1].set_title("N-SPH M-ANEOS ($b=0.73$)")
axs[2].set_title("Stewart M-ANEOS ($b=0.75$)")
axs[3].set_title("N-SPH M-ANEOS ($b=0.75$)")
plot_phase_diagrams(other, other_normalizer)
for ax in [axs[2], axs[3]]:
    # ax.set_xlabel("Min. Pressure (GPa)")
    ax.set_xlabel("Density (kg/m$^3$)")
for ax in [axs[0], axs[2]]:
    # ax.set_ylabel("Num. Disk Particles At Min. Pressure")
    ax.set_ylabel(other_label)
for ax in axs:
    ax.grid(alpha=0.4)
    if None not in other_ylim:
        ax.set_ylim(other_ylim[0], other_ylim[1])
    ax.set_yscale("log")
legend = fig.legend(loc=7, fontsize=16)
for line in legend.get_lines():  # increase line widths in legend
    try:
        line.set_linewidth(4.0)
    except:
        pass
for handle in legend.legendHandles:  # increase marker sizes in legend
    try:
        handle.set_sizes([300.0])
    except:
        pass
fig.subplots_adjust(right=0.80)

letters = list(string.ascii_lowercase)
# annotate with letters in the upper left corner
for i, ax in enumerate(axs):
    ax.text(
        0.05,
        0.95,
        letters[i],
        transform=ax.transAxes,
        fontsize=16,
        fontweight="bold",
        va="top",
    )

plt.savefig("paper1_density_vs_other.png", format='png', dpi=200)
