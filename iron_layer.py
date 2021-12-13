#!/usr/bin/env python3
import csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from src.new_and_old_eos import get_particles

"""
We want to model the iron layer after the GI as a spherical shell.
Volume of a sphere: 4/3 pi r^3
Volume of a spherical shell: 4/3 pi (r_2 - r_1)^3
We need to calculate r_2 - r_1 (can't get it directly from simulation)
Know that density is rho = M / V, so the density of the shell is rho = M / [4/3 pi (r_2 - r_1)^3]
We can get M from the simulation, but what about density?

THERMAL EXPANSIVITY (alpha) : [1 / K]
TEMPERATURE (T) : [K]
DENSITY (rho) : [kg/m3]

DELTA RHO = RHO_0 * ALPHA * DELTA T
RHO_0_iron = 7500

"""

time = 1500
path = "/home/theia/scotthull/1M/formatted_gi_new_eos_b_073"

plt.style.use("dark_background")

labels = {
    0: "Target Silicate",
    1: "Target Iron",
    2: "Impactor Silicate",
    3: "Impactor Iron"
}


def get_time(f):
    formatted_time = None
    with open(f, 'r') as infile:
        reader = csv.reader(infile, delimiter="\t")
        formatted_time = float(next(reader)[0])
    infile.close()
    return round(formatted_time * 0.000277778, 2)  # seconds -> hours


f = path + "/{}.csv".format(time)
time_hrs = get_time(f)
df = pd.read_csv(f, skiprows=2)
impactor_iron = df.loc[(df['tag'] == 3)]
iron_layer = impactor_iron[impactor_iron['radius'] <= 1e7]
iron_layer = iron_layer.sort_values(by=['radius'])
iron_layer_radius = [i / (6371 * 1000) for i in iron_layer['radius']]


fig, axs = plt.subplots(1, 3, figsize=(16, 9), sharex='all', gridspec_kw={"wspace": 0.20})
ax1, ax2, ax3 = axs.flatten()
ax1.scatter(
    iron_layer_radius,
    iron_layer['density'],
    s=2
)
for i in range(1, 6):
    fit = np.polyfit(iron_layer_radius, iron_layer['density'], i)
    p = np.poly1d(fit)
    ax1.plot(
        iron_layer_radius,
        p(iron_layer_radius),
        linewidth=2.0,
        # color='magenta',
        label=i
    )
ax1.legend()
# ax1.annotate(
#     r"y = %f * $x^{%f}$" % (round(fit_A, 2), round(fit_B, 2)),
#     (max(iron_layer_radius) - (.4 * max(iron_layer_radius)), max(fit_y) - (.2 * max(fit_y))),
# )
ax2.scatter(
    iron_layer_radius,
    [i / 1e9 for i in iron_layer['pressure']]
)
ax3.scatter(
    iron_layer_radius,
    iron_layer['temperature']
)
for ax in axs.flatten():
    ax.set_title("{} hrs".format(time_hrs))
    ax.set_xlabel(r'Radius $R_{\bigoplus}$')
    ax.grid(alpha=0.4)
ax1.set_ylabel("Density (kg/m3)"), ax2.set_ylabel("Pressure (GPa)"), ax3.set_ylabel("Temperature (K)")

plt.savefig("iron_layer_profile.png", format='png')
