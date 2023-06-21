#!/usr/bin/env python3
import os
import shutil
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.font_manager as fm
from random import randint
import multiprocessing as mp

from src.identify import ParticleMap
from src.combine import CombineFile
from src.time import get_nearest_iteration_to_time, seconds_to_hours, get_all_iterations_and_times
from src.new_and_old_eos import get_particles, scatter, plot, main_plotting_loop, get_parameter_from_particles
from src.animate import animate
from src.plots3D import get_cube_verts

start_iteration = 0
end_iteration = 200
increment = 5
run_name = "500_mars"
path = "/home/theia/scotthull/Paper2_SPH/gi/{}/{}".format(run_name, run_name)
number_processes = 600
min_normalize_parameter = 0
max_normalize_parameter = 1000
square_scale = 8e6
path_3d = "3D_{}".format(run_name)
if os.path.exists(path_3d):
    shutil.rmtree(path_3d)
os.mkdir(path_3d)

plt.style.use("dark_background")
normalizer = Normalize(min_normalize_parameter, max_normalize_parameter)
cmap = cm.get_cmap('jet')

iterations = []
times = []
target_velocity = []
impactor_velocity = []
tar_imp_velocity_ratio = []
target_entropy = []
delta_s = []
target_silicate_initial = None

for index, iteration in enumerate(np.arange(start_iteration, end_iteration + increment, increment)):
    fig = plt.figure(figsize=(10, 10))
    fig.patch.set_facecolor('xkcd:black')
    ax = fig.add_subplot(111, projection='3d')

    to_fname = "merged_{}_{}.dat".format(iteration, randint(0, 100000))
    cf = CombineFile(num_processes=number_processes, time=iteration, output_path=path, to_fname=to_fname)
    combined_file = cf.combine()
    formatted_time = round(cf.sim_time * 0.000277778, 2)
    df = pd.read_csv(to_fname, header=None, skiprows=2, delimiter="\t")
    target = df[df[1] <= 1]
    impactor = df[df[1] > 1]
    os.remove(to_fname)

    ax.scatter(
        target[3], target[4], target[5], s=5, c=target[13], cmap=cmap, norm=normalizer, alpha=1
    )

    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax._axis3don = False
    # ax.set_box_aspect(aspect=(1, 1, 1))
    ax.set_xticks([])
    # for minor ticks
    ax.set_xticks([], minor=True)
    ax.set_yticks([])
    # for minor ticks
    ax.set_yticks([], minor=True)
    ax.set_zticks([])
    # for minor ticks
    ax.set_zticks([], minor=True)
    ax.set_title(
        str(formatted_time) + " hrs",
        c="white",
    )
    sm = cm.ScalarMappable(norm=normalizer, cmap=cmap)
    sm.set_array([])
    cbaxes = inset_axes(ax, width="30%", height="3%", loc=2, borderpad=1.8)
    cbar = plt.colorbar(sm, cax=cbaxes, orientation='horizontal')
    cbar.ax.tick_params(labelsize=6)
    cbar.ax.set_title(f"Delta Entropy", fontsize=6)
    ax.set_title(f"{formatted_time} hrs")
    ax.set_xlim(-square_scale, square_scale)
    ax.set_ylim(-square_scale, square_scale)
    ax.set_zlim(-square_scale, square_scale)
    for i in get_cube_verts(square_scale=square_scale):
        ax.plot(i[0], i[1], i[2], c='white', linewidth=0.3)
    plt.savefig(f"{path_3d}/{iteration}.png", format='png', dpi=200)

    # get the mean velocity of the target and impactor
    iterations.append(iteration)
    times.append(formatted_time)
    target['velocity'] = np.sqrt(target[6]**2 + target[7]**2 + target[8]**2)
    impactor['velocity'] = np.sqrt(impactor[6]**2 + impactor[7]**2 + impactor[8]**2)
    target_velocity.append(target['velocity'].mean())
    impactor_velocity.append(impactor['velocity'].mean())
    tar_imp_velocity_ratio.append(target['velocity'].mean() / impactor['velocity'].mean())

    # if the iteration is the first iteration, get the initial entropy of the target silicate
    if index == 0:
        target_silicate_initial = target[target[1] % 2 == 0]
        # sort the target silicate initial by column 0
        target_silicate_initial = target_silicate_initial.sort_values(by=0)

    target_silicate = df[df[1] % 2 == 0]
    # sort the target silicate
    target_silicate = target_silicate.sort_values(by=0)
    target_silicate['delta S'] = target_silicate[13] - target_silicate_initial[13]
    # calculate the number of particles with delta S > 500 divided by the total number of particles
    delta_s.append(len(target_silicate[target_silicate['delta S'] > 500]) / len(target_silicate))

    # if this is the last iteration, scatter all delta S
    if iteration == end_iteration:
        target_silicate['radius'] = np.sqrt(target_silicate[3]**2 + target_silicate[4]**2 + target_silicate[5]**2)
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111)
        ax.scatter(target_silicate['radius'] / 1000, target_silicate['delta S'], s=0.1, c='k')
        ax.set_xlabel("Radius (km)")
        ax.set_ylabel("Delta S")
        ax.set_title("Delta S vs Radius")
        ax.grid()
        plt.savefig("delta_s_vs_radius.png", dpi=300)

        animate(
            start_time=start_iteration,
            end_time=end_iteration,
            interval=increment,
            path=path_3d,
            fps=10,
            filename=f"{run_name}_impact_delta_s.mp4",
        )

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)
ax.plot(times, np.array(target_velocity) / 1000, label="Target Velocity")
ax.plot(times, np.array(impactor_velocity) / 1000, label="Impactor Velocity")
ax.set_xlabel("Iteration")
ax.set_ylabel("Velocity (km/s)")
ax.set_title("Velocity vs Iteration")
ax.legend()
ax.grid()
plt.savefig("velocity_vs_iteration.png", dpi=300)

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)
ax.plot(times, tar_imp_velocity_ratio, label="Target / Impactor Velocity")
ax.set_xlabel("Iteration")
ax.set_ylabel("Velocity Ratio")
ax.set_title("Velocity Ratio vs Iteration")
ax.legend()
ax.grid()
plt.savefig("velocity_ratio_vs_iteration.png", dpi=300)

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)
ax.plot(times, delta_s, label="Delta S > 500")
ax.set_xlabel("Iteration")
ax.set_ylabel("Fraction of Particles w/ Delta S > 500")
ax.set_title("Fraction of Particles w/ Delta S > 500 vs Iteration")
ax.legend()
ax.grid()
plt.savefig("delta_s_frac_vs_iteration.png", dpi=300)
