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
import multiprocessing.dummy as mp

from src.identify import ParticleMap
from src.combine import CombineFile
from src.time import get_nearest_iteration_to_time, seconds_to_hours, get_all_iterations_and_times
from src.new_and_old_eos import get_particles, scatter, plot, main_plotting_loop
from src.animate import animate
from src.vapor import calc_vapor_mass_fraction, get_particle_vapor_fraction
from src.report import BuildReports

min_iteration = 0
max_iteration = 3000
sample_interval = 100
parameter = "entropy"
min_normalize = 0
max_normalize = 8000
square_scale = 4e7
number_processes_new_eos = 200
number_processes_old_eos = 200
new_name = "gi_new_eos_b_073"
old_name = "gi_old_eos_b_073"
new_path = "/home/theia/scotthull/1M/{}".format(new_name)
old_path = "/home/theia/scotthull/1M/{}".format(old_name)
phase_curve_new = "/home/theia/scotthull/FDPS_SPH_Orbits/src/phase_data/forstSTS__vapour_curve.txt"
phase_curve_old = "/home/theia/scotthull/FDPS_SPH_Orbits/src/phase_data/duniteN_vapour_curve.txt"
to_path = "/home/theia/scotthull/FDPS_SPH_Orbits/vmf_animate_{}_{}".format(new_name, old_name)
particle_count_path = "/home/theia/scotthull/FDPS_SPH_Orbits/particle_counts_{}_{}".format(new_name, old_name)
new_accessory_path = "/home/theia/scotthull/FDPS_SPH_Orbits/{}_disk_data".format(new_name)
old_accessory_path = "/home/theia/scotthull/FDPS_SPH_Orbits/{}_disk_data".format(old_name)
report_to_dir_new_eos = "/home/theia/scotthull/1M/{}_at_time".format(new_name)
report_to_dir_old_eos = "/home/theia/scotthull/1M/{}_at_time".format(old_name)

for i in [to_path, particle_count_path]:
    if os.path.exists(i):
        shutil.rmtree(i)
    os.mkdir(i)

for i in [report_to_dir_new_eos, report_to_dir_old_eos]:
    os.mkdir(i)

plt.style.use("dark_background")
cmap = cm.get_cmap('jet')
normalizer = Normalize(min_normalize, max_normalize)

cf = CombineFile(num_processes=number_processes_new_eos, time=min_iteration, output_path=new_path)
cf.combine()
min_time = seconds_to_hours(cf.sim_time)
cf = CombineFile(num_processes=number_processes_new_eos, time=max_iteration, output_path=new_path)
cf.combine()
max_time = seconds_to_hours(cf.sim_time)

new_reports = BuildReports(
    accessory_path=new_accessory_path,
    start_time=min_iteration,
    end_time=max_iteration,
    eos_phase_path=phase_curve_new,
    from_dir=new_path,
    number_processes=number_processes_new_eos,
    to_dir=report_to_dir_new_eos
)
old_reports = BuildReports(
    accessory_path=old_accessory_path,
    start_time=min_iteration,
    end_time=max_iteration,
    eos_phase_path=phase_curve_old,
    from_dir=old_path,
    number_processes=number_processes_old_eos,
    to_dir=report_to_dir_old_eos
)

new_vmfs = []
old_vmfs = []
new_times = []
old_times = []

def __loop(time):
    print("VAPORIZATION CODE AT {}".format(time))
    new_particles, new_time, pm_new = get_particles(path=new_path, number_processes=number_processes_new_eos, time=time,
                                            solve=True)
    old_particles, old_time, pm_old = get_particles(path=old_path, number_processes=number_processes_old_eos, time=time,
                                            solve=True)
    new_disk_particles = []
    new_times.append(seconds_to_hours(new_time))
    old_times.append(seconds_to_hours(old_time))
    vmf_new = calc_vapor_mass_fraction(particles=new_particles, phase_path=phase_curve_new, only_disk=True) * 100.0
    vmf_old = calc_vapor_mass_fraction(particles=old_particles, phase_path=phase_curve_old, only_disk=True) * 100.0
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
    ax3.set_ylabel("Silicate Vapor Mass Fraction (%)")
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
        color='magenta',
        linewidth=1.0,
        label="Old EoS"
    )
    ax3.text(
        max_time - (max_time * 0.25),
        50,
        "Silicate New EoS VMF: {}%\nSilicate Old EoS VMF: {}%".format(round(vmf_new, 2), round(vmf_old, 2)),
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

    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(111)
    rects1 = ax.bar(1, [pm_new.num_particles_planet, pm_old.num_particles_planet], width=0.8, label='Planet')
    rects2 = ax.bar(2, [pm_new.num_particles_disk, pm_old.num_particles_disk], width=0.8, label='Disk')
    rects3 = ax.bar(3, [pm_new.num_particles_escaping, pm_old.num_particles_escaping], width=0.8, label='Escaping')

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Particles')
    ax.set_title('Number of Particles in Planet, Disk, Escaping')
    ax.set_xticks([1, 2], ["New EoS", "Old EoS"])
    ax.legend()

    ax.bar_label(rects1, padding=3)
    ax.bar_label(rects2, padding=3)
    ax.bar_label(rects3, padding=3)
    fig.tight_layout()
    plt.savefig(particle_count_path + "/{}.png".format(time), format='png')

pool = mp.Pool(5)
for time in np.arange(min_iteration, max_iteration + sample_interval, sample_interval):
    pool.map(__loop, [time])
pool.close()
pool.join()

animate(
    start_time=min_iteration,
    end_time=max_iteration,
    interval=sample_interval,
    path=to_path,
    fps=10,
    filename="vmf.mp4",
)
