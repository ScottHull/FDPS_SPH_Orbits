#!/usr/bin/env python3
import os
import pandas as pd
import matplotlib.pyplot as plt

# use the colorblind palette
plt.style.use('seaborn-colorblind')

base_path = "/home/theia/scotthull/Paper1_SPH/gi/"
initial_iteration = 0
end_iteration = 1800
# read in all the subdirectory names as the run names
runs = [name for name in os.listdir(base_path) if os.path.isdir(os.path.join(base_path, name))]

fig, axs = plt.subplots(4, 2, figsize=(10, 20), sharex="all", sharey="all")
axs = axs.flatten()

for run in runs:
    circ_path = base_path + f"/{run}/circularized_{run}"
    start_time_df = pd.read_csv(circ_path + f"/{initial_iteration}.csv")
    end_time_df = pd.read_csv(circ_path + f"/{end_iteration}.csv")
    planet = end_time_df[end_time_df["label"] == "PLANET"]
    # get all the particles within 1e7 of the planet
    filtered_planet = end_time_df[end_time_df['radius'] < 1e7]
    filtered_planet = filtered_planet[filtered_planet['tag'] % 2 == 0]  # only silicate particles
    corresponding_start_time = start_time_df[start_time_df["id"].isin(filtered_planet["id"].values)]
    # calculate the fraction of particles where the entropy has increased from the start time to end time by at least 500
    entropy_increase = filtered_planet["entropy"] - corresponding_start_time["entropy"]
    entropy_increase = entropy_increase[entropy_increase > 500]
    entropy_increase_fraction = entropy_increase.count() / filtered_planet["entropy"].count()
    print(f"{run}: {entropy_increase_fraction}")
    plotting_index = 0
    if "new" in run and "b073" in run and "5" in run and 'high' not in run and "500" not in run:
        plotting_index = 0
    elif "old" in run and "b073" in run and "5" in run and 'high' not in run and "500" not in run:
        plotting_index = 1
    elif "new" in run and "b073" in run and "500" in run:
        plotting_index = 2
    elif "old" in run and "b073" in run and "500" in run:
        plotting_index = 3
    elif "new" in run and "b073" in run and "1000" in run:
        plotting_index = 4
    elif "old" in run and "b073" in run and "1000" in run:
        plotting_index = 5
    elif "new" in run and "b073" in run and "2000" in run and "low" not in run:
        plotting_index = 6
    elif "old" in run and "b073" in run and "2000" in run and "low" not in run:
        plotting_index = 7
    else:
        plotting_index = None

    if plotting_index is not None:
        axs[plotting_index].scatter(filtered_planet['radius'] / 1000, filtered_planet['entropy'], s=2, label=run)

axs[0].set_title(r"Stewart M-ANEOS, $\rho_c=5 kg/m^3$")
axs[1].set_title(r"N-SPH M-ANEOS, $\rho_c=5 kg/m^3$")
axs[2].set_title(r"Stewart M-ANEOS, $\rho_c=500 kg/m^3$")
axs[3].set_title(r"N-SPH M-ANEOS, $\rho_c=500 kg/m^3$")
axs[4].set_title(r"Stewart M-ANEOS, $\rho_c=1000 kg/m^3$")
axs[5].set_title(r"N-SPH M-ANEOS, $\rho_c=1000 kg/m^3$")
axs[6].set_title(r"Stewart M-ANEOS, $\rho_c=2000 kg/m^3$")
axs[7].set_title(r"N-SPH M-ANEOS, $\rho_c=2000 kg/m^3$")

for ax in axs:
    ax.set_xlabel("Radius (km)")
    ax.set_ylabel("Entropy (J/K)")
    ax.legend()
    # make the points on the legend larger
    for handle in ax.get_legend().legendHandles:
        handle.set_sizes([30.0])
    ax.grid()

plt.savefig("final_mantle_entropy.png", dpi=300)
