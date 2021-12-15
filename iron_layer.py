#!/usr/bin/env python3
import csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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

start_time = 1
end_time = 1000
# path = "/home/theia/scotthull/1M/formatted_gi_new_eos_b_073"
path = "/Users/scotthull/Desktop"

start_df = pd.read_csv(path + "/{}.csv".format(start_time), skiprows=2, index_col="id")
end_df = pd.read_csv(path + "/{}.csv".format(end_time), skiprows=2, index_col="id")

delta_t_impactor_iron = [end_df['temperature'][index] - start_df['temperature'][index] for index in list(end_df.index)
                         if end_df['radius'][index] < (6371 * 1000) and end_df['tag'][index] == 3]
delta_rho_impactor_iron = [end_df['density'][index] - start_df['density'][index] for index in list(end_df.index) if
                           end_df['radius'][index] < (6371 * 1000) and end_df['tag'][index] == 3]
radius_impactor_iron = [end_df['radius'][index] / (6371 * 1000) for index in list(end_df.index) if
                        end_df['radius'][index] < (6371 * 1000) and end_df['tag'][index] == 3]
delta_t_target_iron = [end_df['temperature'][index] - start_df['temperature'][index] for index in list(end_df.index) if
                       end_df['radius'][index] < (6371 * 1000) and end_df['tag'][index] == 1]
delta_rho_target_iron = [end_df['density'][index] - start_df['density'][index] for index in list(end_df.index) if
                         end_df['radius'][index] < (6371 * 1000) and end_df['tag'][index] == 1]
radius_target_iron = [end_df['radius'][index] / (6371 * 1000) for index in list(end_df.index) if
                      end_df['radius'][index] < (6371 * 1000) and end_df['tag'][index] == 1]
target_iron_mass = sum([
    end_df['mass'][index] for index in list(end_df.index)
                         if end_df['radius'][index] < (6371 * 1000) and end_df['tag'][index] == 1
])
impactor_iron_layer_mass = sum([
    end_df['mass'][index] for index in list(end_df.index)
                         if end_df['radius'][index] < (6371 * 1000) and end_df['tag'][index] == 3
])

print(len([i for i in delta_t_impactor_iron if i > 10 ** 4]) / len(delta_t_impactor_iron))
print(len([i for i in delta_t_target_iron if i > 10 ** 4]) / len(delta_t_target_iron))

print("Core mass: {}\nIron layer mass: {}".format(target_iron_mass, impactor_iron_layer_mass))

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


fig, axs = plt.subplots(1, 2, figsize=(16, 9), gridspec_kw={"hspace": 0.0, "wspace": 0.14})
ax1, ax2 = axs.flatten()

ax1.scatter(
    radius_target_iron, delta_t_target_iron, s=2, label="Target Iron", alpha=0.4
)
ax1.scatter(
    radius_impactor_iron, delta_t_impactor_iron, s=2, label="Impactor Iron", alpha=0.4
)
ax2.scatter(
    radius_target_iron, delta_rho_target_iron, s=2, label="Target Iron", alpha=0.4
)
ax2.scatter(
    radius_impactor_iron, delta_rho_impactor_iron, s=2, label="Impactor Iron", alpha=0.4
)
for ax in axs.flatten():
    ax.grid(alpha=0.4)
    ax.legend()
    ax.set_xlabel(r'Radius $R_{\bigoplus}$')
ax1.set_ylabel("Delta T")
ax2.set_ylabel("Delta Rho")
plt.show()
