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
sample_interval = 20
parameter = "entropy"
min_normalize = 0
max_normalize = 8000
square_scale = 4e7
number_processes = 100
new_path = "/home/theia/scotthull/sph_simulations/gi_new_eos"
old_path = "/home/theia/scotthull/sph_simulations/gi_old_eos"
phase_curve_new = "/home/theia/scotthull/FDPS_SPH_Orbits/src/phase_data/forstSTS__vapour_curve.txt"
phase_curve_old = "/home/theia/scotthull/FDPS_SPH_Orbits/src/phase_data/duniteN_vapour_curve.txt"
to_path = "/home/theia/scotthull/FDPS_SPH_Orbits/vmf_animate"

for i in [to_path]:
    if os.path.exists(i):
        shutil.rmtree(i)
    os.mkdir(i)

plt.style.use("dark_background")
cmap = cm.get_cmap('jet')
normalizer = Normalize(min_normalize, max_normalize)

cf = CombineFile(num_processes=number_processes, time=min_iteration, output_path=new_path)
cf.combine()
min_time = seconds_to_hours(cf.sim_time)
cf = CombineFile(num_processes=number_processes, time=max_iteration, output_path=new_path)
cf.combine()
max_time = seconds_to_hours(cf.sim_time)

new_vmfs = []
old_vmfs = []
new_times = []
old_times = []
for time in np.arange(min_iteration, max_iteration + sample_interval, sample_interval):
    new_particles, new_time = get_particles(path=new_path, number_processes=number_processes, time=time, solve=False)
    old_particles, old_time = get_particles(path=old_path, number_processes=number_processes, time=time, solve=False)
    new_times.append(seconds_to_hours(new_time))
    old_times.append(seconds_to_hours(old_time))
    vmf_new = calc_vapor_mass_fraction(particles=new_particles, phase_path=phase_curve_new, only_disk=False) * 100.0
    vmf_old = calc_vapor_mass_fraction(particles=old_particles, phase_path=phase_curve_old, only_disk=False) * 100.0
    new_vmfs.append(vmf_new)
    old_vmfs.append(vmf_old)
    # Create 2x2 sub plots
    gs = gridspec.GridSpec(2, 2)
    fig = plt.figure(figsize=(10, 10))
    fig.patch.set_facecolor('xkcd:black')
    ax1 = fig.add_subplot(gs[0, 0])  # row 0, col 0
    ax2 = fig.add_subplot(gs[0, 1])  # row 0, col 1
    ax3 = fig.add_subplot(gs[1, :])  # row 1, span all columns
    axs = [ax1, ax2, ax3]
    ax1.set_title("New EoS")
    ax2.set_title("Old EoS")
    ax3.set_xlabel("Time (hrs)")
    ax3.set_ylabel("Disk Vapor Mass Fraction (%)")
    ax1 = plot(fig=fig, axs=axs, index=0, time=new_time, particles=new_particles, cmap=cmap,
               normalizer=normalizer,
               parameter=parameter, square_scale=square_scale, flatten=False)
    ax2 = plot(fig=fig, axs=axs, index=1, time=new_time, particles=old_particles, cmap=cmap,
               normalizer=normalizer,
               parameter=parameter, square_scale=square_scale, flatten=False)
    ax3.plot(
        new_times,
        new_vmfs,
        color='aqua',
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
        max_time - (max_time * 0.25),
        50,
        "Disk New EoS VMF: {}%\nDisk Old EoS VMF: {}%".format(round(vmf_new, 2), round(vmf_old, 2)),
        c="white",
        fontsize=10
    )
    ax3.grid(alpha=0.4)
    ax3.set_xlim(min_time, max_time)
    ax3.set_ylim(0, 60)
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
