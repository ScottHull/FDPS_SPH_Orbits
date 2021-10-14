#!/usr/bin/env python3
import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
from matplotlib.colors import Normalize
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.font_manager as fm

from src.identify import ParticleMap
from src.combine import CombineFile
from src.time import get_nearest_iteration_to_time, seconds_to_hours, get_all_iterations_and_times
from src.new_and_old_eos import get_particles, scatter, plot, main_plotting_loop
from src.animate import animate

min_iteration = 0
max_iteration = 3000
sample_interval = 20
parameter = "entropy"
min_normalize = 0
max_normalize = 8000
square_scale = 4e7
number_processes = 100
new_path = "/home/theia/scotthull/sph_simulations/gi_new_eos"
old_path = "/home/theia/scotthull/sph_simulations/gi_old_eos"
new_iron_hugoniot = "src/phase_data/ironSTS_hugoniot.txt"
old_iron_hugoniot = "src/phase_data/iron__C_hugoniot.txt"

hug_headers = ["density", "pressure", "temperature", "energy", "sound speed", "entropy",
                                          "shock speed", "part speed", "phase"]

new_iron_hugoniot_df = pd.read_csv(new_iron_hugoniot, skiprows=1, delimiter="\t", names=hug_headers)
old_iron_hugoniot_df = pd.read_csv(old_iron_hugoniot, skiprows=1, delimiter="\t", names=hug_headers)


plt.style.use("dark_background")
fig, axs = plt.subplots(1, 2, figsize=(10, 10),
                        gridspec_kw={"hspace": 0.0, "wspace": 0.0})
fig.patch.set_facecolor('xkcd:black')
new_particles, new_time = get_particles(path=new_path, number_processes=number_processes, time=max_iteration, solve=False)
old_particles, old_time = get_particles(path=old_path, number_processes=number_processes, time=max_iteration, solve=False)
axs.flatten()[0].plot(
    new_iron_hugoniot_df["pressure"],
    new_iron_hugoniot_df["entropy"],
    linewidth=1.0,
    color="aqua",
)
axs.flatten()[1].plot(
    old_iron_hugoniot_df["pressure"],
    old_iron_hugoniot_df["entropy"],
    linewidth=1.0,
    color="magenta",
)
axs.flatten()[0].scatter(
    [p.position[0] for p in new_particles if p.position[2] < 0 and p.tag == 1],
    [p.position[1] for p in new_particles if p.position[2] < 0 and p.tag == 1],
    s=0.02,
    marker="o",
    label="Target Iron (New EoS)"
)
axs.flatten()[0].scatter(
    [p.position[0] for p in new_particles if p.position[2] < 0 and p.tag == 3],
    [p.position[1] for p in new_particles if p.position[2] < 0 and p.tag == 3],
    s=0.02,
    marker="o",
    label="Impactor Iron (New EoS)"
)
axs.flatten()[1].scatter(
    [p.position[0] for p in old_particles if p.position[2] < 0 and p.tag == 1],
    [p.position[1] for p in old_particles if p.position[2] < 0 and p.tag == 1],
    s=0.02,
    marker="o",
    label="Target Iron (Old EoS)"
)
axs.flatten()[1].scatter(
    [p.position[0] for p in old_particles if p.position[2] < 0 and p.tag == 3],
    [p.position[1] for p in old_particles if p.position[2] < 0 and p.tag == 3],
    s=0.02,
    marker="o",
    label="Impactor Iron (Old EoS)"
)
for ax in axs.flatten():
    ax.set_xlabel("Pressure")
    ax.set_ylabel("Entropy")
    ax.grid(alpha=0.4)
legend = axs.flatten()[0].legend(fontsize=6)
for handle in legend.legendHandles:
    handle.set_sizes([3.0])

plt.savefig("hugoniot.png", format='png')
