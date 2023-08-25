import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import labellines

isotherms = "/Users/scotthull/Desktop/isotherms2/"

new_phase_path = "src/phase_data/forstSTS__vapour_curve.txt"
new_phase_df = pd.read_fwf(new_phase_path, skiprows=1,
                           names=["temperature", "density_sol_liq", "density_vap", "pressure",
                                  "entropy_sol_liq", "entropy_vap"])
# get 6 random colors
colors = plt.cm.jet(np.linspace(0, 1, 7))
phase_labels = ["", "", "Liq + Sol + Vap", "", "1P Solid", "2P Melt", "1P Vapor"]

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)
# ax.set_xlabel(r"Density (kg/m$^3$)", fontsize=16)
# ax.set_ylabel(r"Pressure (GPa)", fontsize=16)
ax.set_xlabel(r"Entropy", fontsize=16)
# ax.set_ylabel(r"Temperature", fontsize=16)
# ax.set_xlabel("Energy", fontsize=16)
ax.set_ylabel("Pressure (GPa)", fontsize=16)
# increase axis tick label size
ax.tick_params(axis='both', which='major', labelsize=16)
ax.grid()

# ax.plot(new_phase_df['entropy_sol_liq'], new_phase_df['temperature'], linewidth=2, color='blue')
# ax.plot(new_phase_df['entropy_vap'], new_phase_df['temperature'], linewidth=2, color='red')

# # loop through the isotherms directory and plot each isotherm
for i in os.listdir(isotherms):
    if i.endswith(".txt"):
        temp = int(i.replace("for_STS_iso4", "").replace(".txt", ""))
        df = pd.read_fwf(isotherms + i, skiprows=1, sep='\t', header=None)
        # get all unique values of the last (phase) column
        phases = df[5].unique()
        for p in phases:
            if p > 0:
                phase_df = df[df[5] == p]
                color = colors[int(p)]
                # ax.scatter(phase_df[4], [temp for i in phase_df[4]], color=color, label=f"Phase {p} ({temp} K)")
                ax.scatter(phase_df[4], phase_df[1] / 10 ** 9, color=color, label=f"Phase {p} ({temp} K)")

# labellines.labelLines(ax.get_lines(), zorder=2.5, fontsize=16, align=False)

# ax.set_xlim(0, 3000)
# ax.set_ylim(0, 1.7e8)
ax.set_xlim(1000, 12000)
# ax.set_ylim(500, 12000)
# ax.set_xlim(0.5, 2.75e7)
ax.set_ylim(0, 100)
ax.legend(fontsize=16)
plt.tight_layout()
plt.show()

