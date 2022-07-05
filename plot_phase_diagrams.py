import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 14, })
# plt.style.use("dark_background")
plt.style.use('seaborn-colorblind')

iteration = 1800
angles = ['b073', 'b075']
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
        n = "S"
        if runs == "old":
            n = "N"
        for cd in cutoff_densities:
            output_name = fformat.format(cd, angle, runs)
            title_name = tformat.format(cd, angle, n)
            titles.append(title_name)
            names.append(output_name)
    if high:
        output_name = fformat.format(5, angle, "new") + "_high"
        names.append(output_name)
        title_name = tformat.format(5, angle, "S") + "-high"
        titles.append(title_name)
    return names, titles


fig, axs = plt.subplots(2, 2, figsize=(16, 9))
axs = axs.flatten()
b073n_index, b073o_index, b075n_index, b075o_index = 0, 1, 2, 3
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

critical_point_new = max(new_phase_df['temperature'])
critical_point_old = max(old_phase_df['temperature'])



def shade_plot(s, ax):
    phase_df = new_phase_df
    cp = critical_point_new
    if "old" in s:
        phase_df = old_phase_df
        cp = critical_point_old
        # ax.text((3000, 5000), "100% Liquid", fontsize=8)
        # ax.text((3000, 5000), "100% Liquid", fontsize=8)
    ax.plot(
        phase_df['entropy_sol_liq'],
        phase_df['temperature'],
        linewidth=2.0,
        color='black'
    )
    ax.plot(
        phase_df['entropy_vap'],
        phase_df['temperature'],
        linewidth=2.0,
        color='black'
    )
    ax.fill_between(
        x=phase_df['entropy_vap'],
        y1=phase_df['temperature'],
        y2=cp,
        color=colors[-1],
        alpha=0.2,
        label="100% Vapor"
    )
    ax.fill_between(
        x=phase_df['entropy_sol_liq'],
        y1=phase_df['temperature'],
        y2=cp,
        color=colors[-2],
        alpha=0.2,
        label="100% Liquid"
    )
    ax.fill_between(
        x=sorted(list(phase_df['entropy_sol_liq']) + list(new_phase_df['entropy_vap'])),
        y1=cp,
        y2=1e10,
        color=colors[-3],
        alpha=0.2,
        label="Supercritical"
    )
    ax.fill_between(
        x=phase_df['entropy_sol_liq'],
        y1=phase_df['temperature'],
        color=colors[-4],
        edgecolor="none",
        alpha=0.2,
        label="Mixed"
    )
    ax.fill_between(
        x=phase_df['entropy_vap'],
        y1=phase_df['temperature'],
        color=colors[-4],
        edgecolor="none",
        alpha=0.2,
    )


for ax in axs:
    ax.grid(alpha=0.4)
    ax.set_xlim(1800, 12000)
    ax.set_ylim(0, 12500)

for ax in [axs[0], axs[2]]:
    ax.set_ylabel("Temperature (K)")
for ax in [axs[2], axs[3]]:
    ax.set_xlabel("Entropy")
for ax in [axs[0], axs[2]]:
    shade_plot(s="new", ax=ax)
for ax in [axs[1], axs[3]]:
    shade_plot(s="old", ax=ax)


def plot_phase_diagrams():
    for angle in angles:
        high = True
        if angle == "b075":
            high = False
        sims, titles = get_all_sims(angle=angle, high=high)
        for s, t in zip(sims, titles):
            marker = "."
            to_path = base_path + "{}/circularized_{}/{}.csv".format(s, s, iteration)
            df = pd.read_csv(to_path)
            disk = df[df['label'] == "DISK"]
            disk_filtered = disk[disk['circ_entropy_delta'] < 5000]
            temp, entropy = disk_filtered['temperature'], disk_filtered['entropy'] + disk_filtered['circ_entropy_delta']
            cd = int(s.split("_")[0])
            color = colors[cutoff_densities.index(cd)]
            if "high" in s:
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
                temp, entropy, s=1, marker=marker, alpha=0.6, color=color, label=label
            )

for cd in cutoff_densities:
    label = r"$\rho_c = {}$ kg/m$^3$".format(cd)
    color = colors[cutoff_densities.index(cd)]
    axs[0].scatter(
        [], [], marker=".", s=1, alpha=0.6, color=color, label=label
    )
for phase, c in [("100% Vapor", colors[-1]), ("100% Liquid", colors[-2]), ("Supercritical", colors[-3]), ("Mixed", colors[-4])]:
    axs[0].fill_between(
        x=[],
        y1=0,
        y2=0,
        color=c,
        alpha=0.2,
        label=phase
    )


axs[0].set_title("Stewart M-ANEOS ($b=0.73$)")
axs[1].set_title("N-SPH M-ANEOS ($b=0.73$)")
axs[2].set_title("Stewart M-ANEOS ($b=0.75$)")
axs[3].set_title("N-SPH M-ANEOS ($b=0.75$)")
plot_phase_diagrams()
plt.tight_layout()
legend = fig.legend(loc=7, fontsize=16)
for line in legend.get_lines():  # increase line widths in legend
    try:
        line.set_linewidth(4.0)
    except:
        pass
for handle in legend.legendHandles:  # increase marker sizes in legend
    try:
        handle.set_sizes([120.0])
    except:
        pass
fig.subplots_adjust(right=0.82)

plt.savefig("phase_curves.png", format='png', dpi=200)
