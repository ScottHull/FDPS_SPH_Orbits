import pandas as pd
import numpy as np
import string
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 14, })
# plt.style.use("dark_background")
plt.style.use('seaborn-colorblind')

iteration = 1800
runs = [
    ('/home/theia/scotthull/Paper1_SPH/gi/500_b073_new/circularized_500_b073_new', 'Canonical'),
    ('/home/theia/scotthull/Paper2_SPH/gi/500_half_earths/circularized_500_half_earths', 'Half-Earths'),
    # ('', 'Mars')
]
base_path = "/home/theia/scotthull/Paper1_SPH/gi/"
new_phase_path = "src/phase_data/forstSTS__vapour_curve.txt"
old_phase_path = "src/phase_data/duniteN__vapour_curve.txt"

new_phase_df = pd.read_fwf(new_phase_path, skiprows=1,
                           names=["temperature", "density_sol_liq", "density_vap", "pressure",
                                  "entropy_sol_liq", "entropy_vap"])
old_phase_df = pd.read_fwf(old_phase_path, skiprows=1,
                           names=["temperature", "density_sol_liq", "density_vap", "pressure",
                                  "entropy_sol_liq", "entropy_vap"])


fig, axs = plt.subplots(1, 2, figsize=(16, 9), sharex='all', sharey='all')
axs = axs.flatten()
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
        # label="100% Vapor"
    )
    ax.fill_between(
        x=phase_df['entropy_sol_liq'],
        y1=phase_df['temperature'],
        y2=cp,
        color=colors[-2],
        alpha=0.2,
        # label="100% Liquid"
    )
    ax.fill_between(
        x=sorted(list(phase_df['entropy_sol_liq']) + list(new_phase_df['entropy_vap'])),
        y1=cp,
        y2=1e10,
        color=colors[-3],
        alpha=0.2,
        # label="Supercritical"
    )
    ax.fill_between(
        x=phase_df['entropy_sol_liq'],
        y1=phase_df['temperature'],
        color=colors[-4],
        edgecolor="none",
        alpha=0.2,
        # label="Mixed"
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

for ax in [axs[0]]:
    ax.set_ylabel("Temperature (K)")
for ax in axs:
    ax.set_xlabel("Entropy (J/kg/K)")
    shade_plot(s="new", ax=ax)


def plot_phase_diagrams():
        high = True
        # if angle == "b075":
        #     high = False
        for index, (path, t) in enumerate(runs):
            marker = "."
            df = pd.read_csv(path + f"/{iteration}.csv")
            disk = df[df['label'] == "DISK"]
            disk_filtered = disk[disk['circ_entropy_delta'] < 5000]
            temp, entropy = disk_filtered['temperature'], disk_filtered['entropy'] + disk_filtered['circ_entropy_delta']
            axs[index].scatter(
                entropy, temp, s=1, marker=marker, alpha=0.6
            )

for phase, c in [("100% Vapor", colors[-1]), ("100% Liquid", colors[-2]), ("Supercritical", colors[-3]), ("Liquid-Vapor\nMixture", colors[-4])]:
    axs[0].fill_between(
        x=[],
        y1=0,
        y2=0,
        color=c,
        alpha=0.2,
        label=phase
    )


axs[0].set_title("Canonical Impact")
axs[1].set_title("Half-Earths Impact")
plot_phase_diagrams()
plt.tight_layout()
legend = axs[-1].legend(fontsize=16)
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

letters = list(string.ascii_lowercase)
for index, ax in enumerate(axs):
    x1, x2, y1, y2 = ax.axis()
    x_loc = x1 + (0.02 * (x2 - x1))
    y_loc = y2 - (0.08 * (y2 - y1))
    ax.grid(alpha=0.4)
    ax.text(x_loc, y_loc, letters[index], fontweight="bold")

plt.tight_layout()
plt.savefig("phase_curves.png", format='png', dpi=200)
