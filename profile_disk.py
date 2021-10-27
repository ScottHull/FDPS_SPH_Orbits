#!/usr/bin/env python3
import os
import shutil
import numpy as np
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
from src.new_and_old_eos import get_particles, scatter, plot, main_plotting_loop, get_parameter_from_particles, \
    get_common_yrange

sample_time = 3000
fps = 10
parameters = {
    "entropy": {
        "range": [0, 10000],
    },
    "internal_energy": {
        "range": [0, 1e7],
    },
    "pressure": {
        "range": [0, 500e9],
    },
    "temperature": {
        "range": [0, 10000],
    },
}
for key in parameters.keys():
    parameters[key].update({"normalizer": Normalize(parameters[key]['range'][0], parameters[key]['range'][1])})

new_path = "/home/theia/scotthull/sph_simulations/gi_new_eos"
old_path = "/home/theia/scotthull/sph_simulations/gi_new_eos"
number_processes = 100
distance_normalizer = 6371 * 1000

plt.style.use("dark_background")
cmap = cm.get_cmap('jet')

new_particles, new_time = get_particles(path=new_path, number_processes=number_processes, time=sample_time,
                                        solve=True, eos_phase_path="src/phase_data/forstSTS__vapour_curve.txt")
old_particles, old_time = get_particles(path=old_path, number_processes=number_processes, time=sample_time,
                                        solve=True, eos_phase_path="src/phase_data/forstSTS__vapour_curve.txt")
fig, axs = plt.subplots(len(parameters.keys()), 2, figsize=(10, 10), sharex="all",
                            gridspec_kw={"hspace": 0.0, "wspace": 0.08})
fig.patch.set_facecolor('xkcd:black')

tracked_index = 0
for index, parameter in enumerate(parameters.keys()):
    ylims = get_common_yrange(new_particles, old_particles, parameter)
    ax = axs.flatten()[tracked_index].scatter(
        [get_parameter_from_particles(p, "distance") / distance_normalizer for p in new_particles if p.label == "DISK"],
        [get_parameter_from_particles(p, parameter) for p in new_particles if p.label == "DISK"],
        s=0.02,
        marker="o",
    )
    ax.set_ylabel(parameter.replace("_", " ").title())
    ax.grid(alpha=0.4)
    if index + 1 == len(parameters.keys()):
        ax.set_ylabel(r"Radius (1 $R_{\bigoplus}$)")
    tracked_index += 1
    ax = axs.flatten()[tracked_index].scatter(
        [get_parameter_from_particles(p, "distance") / distance_normalizer for p in old_particles if p.label == "DISK"],
        [get_parameter_from_particles(p, parameter) for p in old_particles if p.label == "DISK"],
        s=0.02,
        marker="o",
    )
    ax.grid(alpha=0.4)
    if index + 1 == len(parameters.keys()):
        ax.set_ylabel(r"Radius (1 $R_{\bigoplus}$)")
    tracked_index += 1

axs.flatten()[0].set_title("New EoS")
axs.flatten()[1].set_title("Old EoS")
plt.savefig("disk_profile.png", format='png', dpi=200)

