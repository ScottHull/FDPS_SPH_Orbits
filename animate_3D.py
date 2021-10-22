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
from src.new_and_old_eos import get_particles, scatter, plot, main_plotting_loop, get_parameter_from_particles
from src.animate import animate

min_iteration = 0
max_iteration = 3000
sample_interval = 5
fps = 10
parameter = "entropy"
min_normalize_parameter = 1000
max_normalize_parameter = 8000
path = "/home/theia/scotthull/sph_simulations/gi_new_eos"
to_path = "/home/theia/scotthull/FDPS_SPH_Orbits/3D_animation"
number_processes = 100
square_scale = 2e7

for i in [to_path]:
    if os.path.exists(i):
        shutil.rmtree(i)
    os.mkdir(i)

plt.style.use("dark_background")
normalizer = Normalize(min_normalize_parameter, max_normalize_parameter)
cmap = cm.get_cmap('jet')


for time in np.arange(min_iteration, max_iteration + sample_interval, sample_interval):
    particles, seconds = get_particles(path=path, number_processes=number_processes, time=time)
    fig = plt.figure(figsize=(10, 10))
    fig.patch.set_facecolor('xkcd:black')
    ax = fig.add_subplot(projection='3d')
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.set_box_aspect(aspect=(1, 1, 1))
    ax.scatter(
        [p.position[0] for p in particles],
        [p.position[1] for p in particles],
        [p.position[2] for p in particles],
        s=0.02,
        marker="o",
        c=[cmap(normalizer(get_parameter_from_particles(particle=p, parameter=parameter))) for p in particles],
    )
    ax.set_xticks([])
    # for minor ticks
    ax.set_xticks([], minor=True)
    ax.set_yticks([])
    # for minor ticks
    ax.set_yticks([], minor=True)
    ax.set_zticks([])
    # for minor ticks
    ax.set_zticks([], minor=True)
    ax.set_xlim(-square_scale, square_scale)
    ax.set_ylim(-square_scale, square_scale)
    ax.set_ylim(-square_scale, square_scale)
    # ax.set_box_aspect(1)

    scalebar = AnchoredSizeBar(ax.transData,
                               6378.1 * 1000,
                               r'1 $R_{\bigoplus}$',
                               loc=8,
                               pad=0.3,
                               color='white',
                               frameon=False,
                               size_vertical=1,
                               fontproperties=fm.FontProperties(size=6))

    ax.add_artist(scalebar)
    ax.set_title(
        str(round(seconds_to_hours(seconds), 2)) + " hrs",
        c="white",
    )
    sm = cm.ScalarMappable(norm=normalizer, cmap=cmap)
    sm.set_array([])
    cbaxes = inset_axes(ax, width="30%", height="3%", loc=2, borderpad=1.8)
    cbar = plt.colorbar(sm, cax=cbaxes, orientation='horizontal')
    cbar.ax.tick_params(labelsize=6)
    cbar.ax.set_title(parameter.title(), fontsize=6)
    plt.savefig(to_path + "/{}.png".format(time), format='png', dpi=200)

animate(
    start_time=min_iteration,
    end_time=max_iteration,
    interval=sample_interval,
    path=to_path,
    fps=fps,
    filename="3D_impact.mp4",
)

