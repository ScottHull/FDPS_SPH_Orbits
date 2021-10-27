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
from src.animate import animate

start_time = 0
end_time = 1000
sample_interval = 10
to_path = "/home/theia/scotthull/FDPS_SPH_Orbits/profile_disk"
fps = 30
parameters = {
    "entropy": {
        "yrange": [1000, 12000],
    },
    "internal_energy": {
        "yrange": [0.1e7, 6.5e7],
    },
    "pressure": {
        "yrange": [0, 3],
    },
    "temperature": {
        "yrange": [2000, 20000],
    },
}
for key in parameters.keys():
    parameters[key].update({"normalizer": Normalize(parameters[key]['yrange'][0], parameters[key]['yrange'][1])})

new_path = "/home/theia/scotthull/sph_simulations/gi_new_eos"
old_path = "/home/theia/scotthull/sph_simulations/gi_new_eos"
number_processes = 100
distance_normalizer = 6371 * 1000
labels = {
        0: "Target Silicate", 1: "Target Iron", 2: "Impactor Silicate", 3: "Impactor Iron",
    }

for i in [to_path]:
    if os.path.exists(i):
        shutil.rmtree(i)
    os.mkdir(i)

def scatter(ax, particles):
    for i in range(0, 4):
        ax = ax.scatter(
            [get_parameter_from_particles(p, "distance") / distance_normalizer for p in particles if p.label == "DISK" and p.tag == i],
            [get_parameter_from_particles(p, parameter) for p in particles if p.label == "DISK" and p.tag == i],
            s=0.4,
            marker="o",
            c="magenta",
            label=labels[i]
        )
        if i == 0:
            legend = ax.legend(loc='upper left', fontsize=6)
            for handle in legend.legendHandles:
                handle.set_sizes([3.0])
    return ax

plt.style.use("dark_background")
cmap = cm.get_cmap('jet')

for sample_time in np.arange(start_time, end_time + sample_interval, sample_interval):
    new_particles, new_time = get_particles(path=new_path, number_processes=number_processes, time=sample_time,
                                            solve=True, eos_phase_path="src/phase_data/forstSTS__vapour_curve.txt")
    # old_particles, old_time = get_particles(path=old_path, number_processes=number_processes, time=sample_time,
    #                                         solve=True, eos_phase_path="src/phase_data/forstSTS__vapour_curve.txt")
    fig, axs = plt.subplots(len(parameters.keys()), 1, figsize=(12, 20), sharex="all",
                                gridspec_kw={"hspace": 0.0, "wspace": 0.14})
    fig.patch.set_facecolor('xkcd:black')

    tracked_index = 0
    for index, parameter in enumerate(parameters.keys()):
        # ylims = get_common_yrange(new_particles, old_particles, parameter)
        ax = axs.flatten()[tracked_index]
        ax.set_ylabel(parameter.replace("_", " ").title())
        ax.set_xlim(0, 60)
        ax.set_ylim(parameters[parameter]['yrange'][0], parameters[parameter]['yrange'][1])
        ax.grid(alpha=0.4)
        scatter(ax=ax, particles=new_particles)
        if index + 1 == len(parameters.keys()):
            ax.set_xlabel(r"Radius (1 $R_{\bigoplus}$)")
        tracked_index += 1

        # ax = axs.flatten()[tracked_index]
        # ax.grid(alpha=0.4)
        # ax.set_xlim(0, 60)
        # if index + 1 == len(parameters.keys()):
        #     ax.set_xlabel(r"Radius ($R_{\bigoplus}$)")
        # ax = axs.flatten()[tracked_index].scatter(
        #     [get_parameter_from_particles(p, "distance") / distance_normalizer for p in old_particles if p.label == "DISK"],
        #     [get_parameter_from_particles(p, parameter) for p in old_particles if p.label == "DISK"],
        #     s=0.4,
        #     marker="o",
        # )
        # tracked_index += 1

    axs.flatten()[0].set_title("Disk Particles - New EoS ({} hrs)".format(seconds_to_hours(new_time)))
    # axs.flatten()[1].set_title("Disk Particles - Old EoS")
    plt.savefig(to_path + "/disk_profile.png", format='png', dpi=200)

animate(
    start_time=start_time,
    end_time=end_time,
    interval=sample_interval,
    path=to_path,
    fps=fps,
    filename="profile_disk.mp4",
)

