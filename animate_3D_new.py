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
from random import randint
import multiprocessing as mp

from src.identify import ParticleMap
from src.combine import CombineFile
from src.time import get_nearest_iteration_to_time, seconds_to_hours, get_all_iterations_and_times
from src.new_and_old_eos import get_particles, scatter, plot, main_plotting_loop, get_parameter_from_particles
from src.animate import animate
from src.plots3D import get_cube_verts

min_iteration = 0
max_iteration = 3000
increment = 5
fps = 30
parameter = "entropy"
min_normalize_parameter = 1000
max_normalize_parameter = 10000
run_name = "500_b073_new"
path = "/home/theia/scotthull/Paper1_SPH/gi/{}/{}".format(run_name, run_name)
to_path = "3D_{}".format(run_name)
number_processes = 200
square_scale = 2e7

verts = get_cube_verts(square_scale=square_scale)

for i in [to_path]:
    if not os.path.exists(i):
        os.mkdir(i)

plt.style.use("dark_background")
normalizer = Normalize(min_normalize_parameter, max_normalize_parameter)
cmap = cm.get_cmap('jet')


def plot_iteration(iteration):
    to_fname = "merged_{}_{}.dat".format(iteration, randint(0, 100000))
    cf = CombineFile(num_processes=number_processes, time=iteration, output_path=path, to_fname=to_fname)
    combined_file = cf.combine_df()
    combined_file.columns = [
        'id', 'tag', 'mass', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'density', 'internal_energy', 'pressure',
        'potential_energy', 'entropy', 'temperature'
    ]
    formatted_time = round(cf.sim_time * 0.000277778, 2)
    f = os.getcwd() + "/{}".format(to_fname)

    fig = plt.figure(figsize=(10, 10))
    fig.patch.set_facecolor('xkcd:black')
    ax = fig.add_subplot(111, projection='3d')

    ax.scatter(
        combined_file['x'], combined_file['y'], combined_file['z'], s=5, c=cmap(normalizer(combined_file[parameter]))
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
    cbar.ax.set_title(parameter.title(), fontsize=6)
    ax.set_xlim(-square_scale, square_scale)
    ax.set_ylim(-square_scale, square_scale)
    ax.set_zlim(-square_scale, square_scale)
    for i in get_cube_verts(square_scale=square_scale):
        ax.plot(i[0], i[1], i[2], c='white', linewidth=0.3)
    plt.savefig(to_path + "/{}.png".format(iteration), format='png', dpi=200)

ldir = os.listdir(to_path)
pool = mp.Pool(10)
pool.map(plot_iteration, [iteration for iteration in np.arange(min_iteration, max_iteration + increment, increment) if "{}.png".format(iteration) not in ldir])
pool.close()
pool.join()

animate(
    start_time=min_iteration,
    end_time=max_iteration,
    interval=increment,
    path=to_path,
    fps=fps,
    filename="{}_impact.mp4".format(run_name),
)
