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


# Function to calculate the exponential with constants a and b
def exponential(x, a, b):
    return a * np.exp(b * x)


f = path + "/{}.csv".format(time)
time_hrs = get_time(f)
df = pd.read_csv(f, skiprows=2)
impactor_iron = df.loc[(df['tag'] == 3)]
iron_layer = impactor_iron[impactor_iron['radius'] <= 1e7]
iron_layer = iron_layer.sort_values(by=['radius'])

parameters, covariance = curve_fit(exponential, [i / (6371 * 1000) for i in iron_layer['radius']], impactor_iron['density'])
fit_A = parameters[0]
fit_B = parameters[1]
fit_y = exponential([i / (6371 * 1000) for i in iron_layer['radius']], fit_A, fit_B)

fig, axs = plt.subplots(1, 3, figsize=(16, 9), sharex='all', gridspec_kw={"wspace": 0.20})
ax1, ax2, ax3 = axs.flatten()
ax1.scatter(
    [i / (6371 * 1000) for i in iron_layer['radius']],
    iron_layer['density']
)
ax1.plot(
    [i / (6371 * 1000) for i in iron_layer['radius']],
    fit_y,
    linewidth=2.0,
    color='magenta'
)
ax2.scatter(
    [i / (6371 * 1000) for i in iron_layer['radius']],
    iron_layer['entropy']
)
ax3.scatter(
    [i / (6371 * 1000) for i in iron_layer['radius']],
    iron_layer['temperature']
)
for ax in axs.flatten():
    ax.set_title("{} hrs".format(time_hrs))
    ax.set_xlabel(r'Radius $R_{\bigoplus}$')
    ax.grid(alpha=0.4)
ax1.set_ylabel("Density"), ax2.set_ylabel("Entropy"), ax3.set_ylabel("Temperature")

plt.savefig("iron_layer_profile.png", format='png')
