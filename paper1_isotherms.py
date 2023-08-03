import os
import pandas as pd
import matplotlib.pyplot as plt
import labellines

new_phase_path = "src/phase_data/forstSTS__vapour_curve.txt"
new_phase_df = pd.read_fwf(new_phase_path, skiprows=1,
                           names=["temperature", "density_sol_liq", "density_vap", "pressure",
                                  "entropy_sol_liq", "entropy_vap"])


fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)
ax.set_xlabel(r"Density (kg/m$^3$)", fontsize=16)
ax.set_ylabel(r"Pressure (GPa)", fontsize=16)
# increase axis tick label size
ax.tick_params(axis='both', which='major', labelsize=16)

# get all unique values of temperature
temps = new_phase_df['temperature'].unique()
# loop through each temperature
for index, t in enumerate(temps):
    if index % 2 == 0:
        # get the subset of the dataframe for this temperature
        df = new_phase_df[new_phase_df['temperature'] == t]
        # get density/presure values for the liquid and vapor
        density_sol_liq = df['density_sol_liq'].values[0]
        density_vap = df['density_vap'].values[0]
        pressure = df['pressure'].values[0] / 10 ** 9

        ax.plot(
            density_sol_liq, pressure, color='black', linewidth=2.0, label=f"{t} K"
        )
        ax.plot(
            density_vap, pressure, color='red', linewidth=2.0
        )

labellines.labelLines(ax.get_lines(), zorder=2.5, align=True, fontsize=12)
ax.grid()
plt.tight_layout()
plt.show()

