#!/usr/bin/env python3
import os
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.font_manager as fm

from src.identify import ParticleMap
from src.combine import CombineFile
from src.time import get_nearest_iteration_to_time, seconds_to_hours, get_all_iterations_and_times
from src.new_and_old_eos import get_particles, get_peak
from src.animate import animate

min_iteration = 0
max_iteration = 300
sample_interval = 5
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

new_iron_hugoniot_df = pd.read_fwf(new_iron_hugoniot, skiprows=1, names=hug_headers)
old_iron_hugoniot_df = pd.read_fwf(old_iron_hugoniot, skiprows=1, names=hug_headers)

new_peak, old_peak = get_peak(save=True, parameter="pressure", min_iteration=min_iteration, max_iteration=max_iteration,
                    interval=sample_interval, new_path=new_path, old_path=old_path, number_processes=number_processes)

def pa_to_gpa(pa):
    """
    Converts from Pa to GPa
    :param pa:
    :return:
    """
    return pa * (10 ** -9)

plt.style.use("dark_background")
fig, axs = plt.subplots(2, 2, figsize=(12, 10),
                            gridspec_kw={"hspace": 0.12, "wspace": 0.12})
fig.patch.set_facecolor('xkcd:black')

axs.flatten()[0].scatter(
    [pa_to_gpa(new_peak[p]["pressure"]) for p in new_peak.keys() if new_peak[p]["tag"] == 0],
    [new_peak[p]["entropy"] for p in new_peak.keys() if new_peak[p]["tag"] == 0],
    s=0.02,
    marker="o",
    label="Target Silicate"
)
axs.flatten()[0].scatter(
    [pa_to_gpa(new_peak[p]["pressure"]) for p in new_peak.keys() if new_peak[p]["tag"] == 1],
    [new_peak[p]["entropy"] for p in new_peak.keys() if new_peak[p]["tag"] == 1],
    s=0.02,
    marker="o",
    label="Target Iron"
)
axs.flatten()[0].scatter(
    [pa_to_gpa(new_peak[p]["pressure"]) for p in new_peak.keys() if new_peak[p]["tag"] == 2],
    [new_peak[p]["entropy"] for p in new_peak.keys() if new_peak[p]["tag"] == 2],
    s=0.02,
    marker="o",
    label="Impactor Silicate"
)
axs.flatten()[0].scatter(
    [pa_to_gpa(new_peak[p]["pressure"]) for p in new_peak.keys() if new_peak[p]["tag"] == 3],
    [new_peak[p]["entropy"] for p in new_peak.keys() if new_peak[p]["tag"] == 3],
    s=0.02,
    marker="o",
    label="Impactor Iron"
)
axs.flatten()[1].scatter(
    [pa_to_gpa(old_peak[p]["pressure"]) for p in old_peak.keys() if old_peak[p]["tag"] == 0],
    [old_peak[p]["entropy"] for p in old_peak.keys() if old_peak[p]["tag"] == 0],
    s=0.02,
    marker="o",
    label="Target Silicate"
)
axs.flatten()[1].scatter(
    [pa_to_gpa(old_peak[p]["pressure"]) for p in old_peak.keys() if old_peak[p]["tag"] == 1],
    [old_peak[p]["entropy"] for p in old_peak.keys() if old_peak[p]["tag"] == 1],
    s=0.02,
    marker="o",
    label="Target Iron"
)
axs.flatten()[1].scatter(
    [pa_to_gpa(old_peak[p]["pressure"]) for p in old_peak.keys() if old_peak[p]["tag"] == 2],
    [old_peak[p]["entropy"] for p in old_peak.keys() if old_peak[p]["tag"] == 2],
    s=0.02,
    marker="o",
    label="Impactor Silicate"
)
axs.flatten()[1].scatter(
    [pa_to_gpa(old_peak[p]["pressure"]) for p in old_peak.keys() if old_peak[p]["tag"] == 3],
    [old_peak[p]["entropy"] for p in old_peak.keys() if old_peak[p]["tag"] == 3],
    s=0.02,
    marker="o",
    label="Impactor Iron"
)

axs.flatten()[0].plot(
    new_iron_hugoniot_df["pressure"],
    new_iron_hugoniot_df["entropy"],
    linewidth=1.0,
    color="magenta",
)
axs.flatten()[1].plot(
    old_iron_hugoniot_df["pressure"],
    old_iron_hugoniot_df["entropy"],
    linewidth=1.0,
    color="magenta",
)

axs.flatten()[2].hist(
    [seconds_to_hours(new_peak[p]["time"]) for p in new_peak.keys()],
    bins='auto'
)
axs.flatten()[3].hist(
    [seconds_to_hours(old_peak[p]["time"]) for p in old_peak.keys()],
    bins='auto'
)

for ax in [axs.flatten()[0], axs.flatten()[1]]:
    ax.set_xlabel("Max. Pressure"),
    ax.grid(alpha=0.4)
    ax.set_box_aspect(1)
    ax.set_xlim(0, 0.5e12)
for ax in [axs.flatten()[2], axs.flatten()[3]]:
    ax.set_xlabel("Time of Max Pressure (hrs)")
    ax.grid(alpha=0.4)
    ax.set_box_aspect(1)
legend = axs.flatten()[0].legend(fontsize=6)
for handle in legend.legendHandles:
    handle.set_sizes([3.0])
axs.flatten()[0].set_ylabel("Corresponding Entropy")

plt.savefig("peaks.png", format='png')
