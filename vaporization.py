#!/usr/bin/env python3
import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
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
from src.vapor import calc_vapor_mass_fraction, get_particle_vapor_fraction

min_iteration = 0
max_iteration = 3000
sample_interval = 5
parameter = "entropy"
min_normalize = 0
max_normalize = 8000
square_scale = 3e7
number_processes = 100
new_path = "/home/theia/scotthull/sph_simulations/gi_new_eos"
old_path = "/home/theia/scotthull/sph_simulations/gi_old_eos"
phase_curve_new = "/home/theia/scotthull/FDPS_SPH_Orbits/src/phase_data/forstSTS__vapour_curve.txt"
phase_curve_old = "/home/theia/scotthull/FDPS_SPH_Orbits/src/phase_data/duniteN_vapour_curve.txt"
to_path = "/home/theia/scotthull/FDPS_SPH_Orbits/vmf_animate"

plt.style.use("dark_background")
cmap = cm.get_cmap('jet')
normalizer = Normalize(min_normalize, max_normalize)

min_time = seconds_to_hours(CombineFile(num_processes=number_processes, time=min_iteration, output_path=new_path).combine().sim_time)
max_time = seconds_to_hours(CombineFile(num_processes=number_processes, time=max_iteration, output_path=new_path).combine().sim_time)

new_vmfs = []
old_vmfs = []
new_times = []
old_times = []
for time in np.arange(min_iteration, max_iteration + sample_interval, sample_interval):
    new_particles, new_time = get_particles(path=new_path, number_processes=number_processes, time=time)
    old_particles, old_time = get_particles(path=old_path, number_processes=number_processes, time=time)
    new_times.append(seconds_to_hours(new_time))
    old_times.append(seconds_to_hours(old_time))
    vmf_new = calc_vapor_mass_fraction(particles=new_particles, phase_path=phase_curve_new) * 100.0
    vmf_old = calc_vapor_mass_fraction(particles=old_particles, phase_path=phase_curve_old) * 100.0
    new_vmfs.append(vmf_new)
    old_vmfs.append(vmf_old)
    # Create 2x2 sub plots
    gs = gridspec.GridSpec(2, 2)
    fig = plt.figure(figsize=(16, 9))
    fig.patch.set_facecolor('xkcd:black')
    ax1 = fig.add_subplot(gs[0, 0])  # row 0, col 0
    ax2 = fig.add_subplot(gs[0, 1])  # row 0, col 1
    ax3 = fig.add_subplot(gs[1, :])  # row 1, span all columns
    axs = [ax1, ax2, ax3]
    ax1.set_title("New EoS")
    ax2.set_title("Old EoS")
    ax3.set_xlabel("Time (hrs)")
    ax3.set_ylabel("Vapor Mass Fraction (%)")
    ax1 = plot(fig=fig, axs=axs, index=0, time=new_time, particles=new_particles, cmap=cmap,
               normalizer=normalizer,
               parameter=parameter, square_scale=square_scale)
    ax2 = plot(fig=fig, axs=axs, index=1, time=new_time, particles=old_particles, cmap=cmap,
               normalizer=normalizer,
               parameter=parameter, square_scale=square_scale)
    ax3.plot(
        new_times,
        new_vmfs,
        color='blue',
        linewidth=1.0,
        label="New EoS"
    )
    ax3.plot(
        old_times,
        old_vmfs,
        color='pink',
        linewidth=1.0,
        label="Old EoS"
    )
    ax3.text(
        max_time - (max_time * 0.1),
        85,
        "New EoS VMF: {}%\nOld EoS VMF: {}%".format(vmf_new, vmf_old),
        c="white",
        fontsize=10
    )
    ax3.grid()
    ax3.set_xlim(min_iteration, max_iteration)
    ax3.set_ylim(0, 100)
    sm = cm.ScalarMappable(norm=normalizer, cmap=cmap)
    sm.set_array([])
    # cbar = fig.colorbar(sm, ax=axs.flatten()[1])
    cbaxes = inset_axes(ax1, width="30%", height="3%", loc=2, borderpad=1.8)
    cbar = plt.colorbar(sm, cax=cbaxes, orientation='horizontal')
    cbar.ax.tick_params(labelsize=6)
    # cbar.ax.xaxis.set_ticks_position('top')
    cbar.ax.set_title(parameter.replace("_", " ").title(), fontsize=6)
    legend = ax3.legend(fontsize=6, loc='upper left')
    plt.savefig(to_path + "/{}.png".format(time), format='png', dpi=200)

animate(
    start_time=min_iteration,
    end_time=max_iteration,
    interval=sample_interval,
    path=to_path,
    fps=10,
    filename="vmf.mp4",
)
