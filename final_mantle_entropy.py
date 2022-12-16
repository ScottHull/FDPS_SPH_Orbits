#!/usr/bin/env python3
import os
import pandas as pd
import matplotlib.pyplot as plt

# use the colorblind palette
plt.style.use('seaborn-colorblind')

base_path = "/home/theia/scotthull/Paper1_SPH/gi/"
end_iteration = 1800
# read in all the subdirectory names as the run names
runs = [name for name in os.listdir(base_path) if os.path.isdir(os.path.join(base_path, name))]

fig, axs = plt.subplots(2, 2, figsize=(10, 10), sharex="all", sharey="all")
axs = axs.flatten()

for run in runs:
    circ_path = base_path + f"/{run}/circularized_{run}"
    end_time_df = pd.read_csv(circ_path + f"/{end_iteration}.csv")
    planet = end_time_df[end_time_df["label"] == "PLANET"]
    # get all the particles within 1e7 of the planet
    filtered_planet = end_time_df[end_time_df['radius'] < 1e7]
    filtered_planet = filtered_planet[filtered_planet['tag'] % 2 == 0]  # only silicate particles
    plotting_index = 0
    if "new" in run and "b073" in run:
        plotting_index = 0
    elif "old" in run and "b073" in run:
        plotting_index = 1
    elif "new" in run and "b075" in run:
        plotting_index = 2
    elif "old" in run and "b075" in run:
        plotting_index = 3
    else:
        plotting_index = None

    if plotting_index is not None:
        axs[plotting_index].scatter(filtered_planet['radius'] / 1000, filtered_planet['entropy'], s=2, label=run)

axs[0].set_title("Stewart M-ANEOS, b = 0.73")
axs[1].set_title("Stewart M-ANEOS, b = 0.73")
axs[2].set_title("N-SPH M-ANEOS, b = 0.75")
axs[3].set_title("N-SPH M-ANEOS, b = 0.75")

for ax in axs:
    ax.set_xlabel("Radius (km)")
    ax.set_ylabel("Entropy (J/K)")
    ax.legend()
    ax.grid()

plt.tight_layout()

plt.savefig("final_mantle_entropy.png", dpi=300)
