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
from src.new_and_old_eos import get_particles, scatter, plot
# from src.animate import animate

min_iteration = 0
max_iteration = 3000
sample_interval = 5
parameter = "entropy"
min_normalize = 0
max_normalize = 8000
square_scale = 1e7
number_processes = 100
to_path = "/home/theia/scotthull/FDPS_SPH_Orbits/new_and_old_animate"
new_path = "/home/theia/scotthull/sph_simulations/gi_new_eos"
old_path = "/home/theia/scotthull/sph_simulations/gi_old_eos"

if os.path.exists(to_path):
    shutil.rmtree(to_path)
os.mkdir(to_path)

plt.style.use("dark_background")
cmap = cm.get_cmap('jet')
normalizer = Normalize(min_normalize, max_normalize)

for time in np.arange(min_iteration, max_iteration + sample_interval, sample_interval):
    fig, axs = plt.subplots(2, 2, figsize=(10, 10),
                            gridspec_kw={"hspace": 0.0, "wspace": 0.08})
    fig.patch.set_facecolor('xkcd:black')
    new_particles, new_time = get_particles(path=new_path, number_processes=number_processes, time=time)
    old_particles, old_time = get_particles(path=old_path, number_processes=number_processes, time=time)
    ax1 = plot(fig=fig, axs=axs, index=0, time=new_time, particles=new_particles, cmap=cmap, normalizer=normalizer,
               parameter=parameter, square_scale=square_scale)
    ax2 = plot(fig=fig, axs=axs, index=1, time=new_time, particles=old_particles, cmap=cmap, normalizer=normalizer,
               parameter=parameter, square_scale=square_scale)
    ax3 = scatter(fig=fig, axs=axs, index=2, particles=new_particles, parameter=parameter, square_scale=square_scale)
    ax4 = scatter(fig=fig, axs=axs, index=3, particles=old_particles, parameter=parameter, square_scale=square_scale)
    axs.flatten()[2].set_ylabel(parameter.replace("_", " ").title())

    axs.flatten()[0].set_title("New EoS")
    axs.flatten()[1].set_title("Old EoS")
    sm = cm.ScalarMappable(norm=normalizer, cmap=cmap)
    sm.set_array([])
    # cbar = fig.colorbar(sm, ax=axs.flatten()[1])
    cbaxes = inset_axes(ax1, width="30%", height="3%", loc=2, borderpad=1.8)
    cbar = plt.colorbar(sm, cax=cbaxes, orientation='horizontal')
    cbar.ax.tick_params(labelsize=6)
    # cbar.ax.xaxis.set_ticks_position('top')
    cbar.ax.set_title(parameter.replace("_", " ").title(), fontsize=6)
    legend = ax3.legend(fontsize=6)
    for handle in legend.legendHandles:
        handle.set_sizes([3.0])
    ax4.set_yticks([])
    ax4.set_yticks([], minor=True)
    ax3_ymin, ax3_ymax = ax3.get_ylim()
    ax4_ymin, ax4_ymax = ax4.get_ylim()
    lims = [ax3_ymin, ax3_ymax, ax4_ymin, ax4_ymax]
    scatter_range = [min(lims), max(lims)]
    inc = (scatter_range[1] - scatter_range[0]) * 0.1
    ax3.set_ylim(0, scatter_range[1] + inc)
    ax4.set_ylim(0, scatter_range[1] + inc)
    ax3.set_xlabel("Radius (m)")
    ax4.set_xlabel("Radius (m)")
    plt.savefig(to_path + "/{}.png".format(time), format='png', dpi=200)

animate(
    start_time=min_iteration,
    end_time=max_iteration,
    interval=sample_interval,
    path=to_path,
    fps=10,
    filename="impact_geometry.mp4",
)
