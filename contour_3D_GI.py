#!/usr/bin/env python3
import os
import shutil
from random import randint
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.font_manager as fm

from src.identify import ParticleMap
from src.combine import CombineFile
from src.animate import animate
from src.new_and_old_eos import seconds_to_hours, get_particles, get_parameter_from_particles
from src.plots3D import get_cube_verts

start_time = 0
end_time = 1500
interval = 5
number_processes = 200
min_norm = 0
max_norm = 10000
parameter = "entropy"
square_scale = 1e7
path = "/home/theia/scotthull/gi_new_eos"
eos = "src/phase_data/forst_STS.rho_u.txt"
output = "/home/theia/scotthull/FDPS_SPH_Orbits/3D_contour_GI"
output_3D = "/home/theia/scotthull/FDPS_SPH_Orbits/mapped_3D_contour_GI"

normalizer = Normalize(min_norm, max_norm)
cmap = cm.get_cmap('jet')
sm = cm.ScalarMappable(norm=normalizer, cmap=cmap)

for o in [output, output_3D]:
    if os.path.exists(o):
        shutil.rmtree(o)
    os.mkdir(o)

eos_df = pd.read_fwf(eos, skiprows=2, header=None)
eos_density = eos_df[0]
eos_internal_energy = eos_df[1]
eos_entropy = eos_df[5]

cf_end = CombineFile(num_processes=number_processes, time=end_time, output_path=path)
combined_file_end = cf_end.combine()
formatted_time_end = cf_end.sim_time
f = os.getcwd() + "/merged_{}.dat".format(end_time)
pm_end = ParticleMap(path=f, center=False, relative_velocity=False)
select_particles = pm_end.collect_particles()
pm_end.solve(particles=select_particles, phase_path="src/phase_data/forstSTS__vapour_curve.txt")
os.remove(f)

end = {}
high_entropy = {}
for p in select_particles:
    end.update({p.particle_id: p.label})
    if p.entropy > 8000 and p.label == "DISK":
        high_entropy.update({p.particle_id: p.entropy})
high_entropy_ids = list(high_entropy.keys())
rand_select = [high_entropy_ids[randint(0, len(high_entropy_ids) - 1)] for i in range(0, 5)]
colors = ["black", "red", "blue", "green", "brown"]
c_dict = {}
for index, i in enumerate(rand_select):
    c_dict.update({i: colors[index]})

prev_particles = {}

for time in np.arange(start_time, end_time + interval, interval):
    particles, seconds = get_particles(path=path, number_processes=number_processes, time=time)
    select_particles_at_time = [p for p in particles if p.particle_id in rand_select]
    fig = plt.figure(figsize=(16, 9))
    # fig.patch.set_facecolor('xkcd:black')
    ax1, ax2 = fig.add_subplot(121, projection='3d'), fig.add_subplot(122, projection='3d')
    axs = [ax1, ax2]
    ax1.scatter(
        [p.position[0] for p in particles if
         abs(p.position[0]) <= square_scale and abs(p.position[1]) <= square_scale and abs(
             p.position[2]) <= square_scale],
        [p.position[1] for p in particles if
         abs(p.position[0]) <= square_scale and abs(p.position[1]) <= square_scale and abs(
             p.position[2]) <= square_scale],
        [p.position[2] for p in particles if
         abs(p.position[0]) <= square_scale and abs(p.position[1]) <= square_scale and abs(
             p.position[2]) <= square_scale],
        s=0.04,
        marker="o",
        c=[cmap(normalizer(get_parameter_from_particles(particle=p, parameter=parameter))) for p in particles if
           abs(p.position[0]) <= square_scale and abs(p.position[1]) <= square_scale and abs(
               p.position[2]) <= square_scale],
        alpha=0.3
    )
    ax2.scatter(
        [p.position[0] for p in particles if
         abs(p.position[0]) <= 5 * square_scale and abs(p.position[1]) <= 5 * square_scale and abs(
             p.position[2]) <= 5 * square_scale],
        [p.position[1] for p in particles if
         abs(p.position[0]) <= 5 * square_scale and abs(p.position[1]) <= 5 * square_scale and abs(
             p.position[2]) <= 5 * square_scale],
        [p.position[2] for p in particles if
         abs(p.position[0]) <= 5 * square_scale and abs(p.position[1]) <= 5 * square_scale and abs(
             p.position[2]) <= 5 * square_scale],
        s=0.5,
        marker="o",
        c=[cmap(normalizer(get_parameter_from_particles(particle=p, parameter=parameter))) for p in particles if
           abs(p.position[0]) <= 5 * square_scale and abs(p.position[1]) <= 5 * square_scale and abs(
               p.position[2]) <= 5 * square_scale],
        alpha=0.3
    )
    for ax in axs:
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
            str(round(seconds_to_hours(seconds), 2)) + " hrs",
            c="black",
        )
        sm = cm.ScalarMappable(norm=normalizer, cmap=cmap)
        sm.set_array([])
        cbaxes = inset_axes(ax, width="30%", height="3%", loc=2, borderpad=1.8)
        cbar = plt.colorbar(sm, cax=cbaxes, orientation='horizontal')
        cbar.ax.tick_params(labelsize=6)
        cbar.ax.set_title(parameter.title(), fontsize=6)
        for p in select_particles_at_time:
            ax.scatter(
                [p.position[0]], [p.position[1]], [p.position[2]], s=20, c=c_dict[p.particle_id], marker="*"
            )
    ax1.set_xlim(-square_scale, square_scale)
    ax1.set_ylim(-square_scale, square_scale)
    ax1.set_zlim(-square_scale, square_scale)
    for i in get_cube_verts(square_scale=square_scale):
        ax1.plot(i[0], i[1], i[2], c='black', linewidth=0.3)
        ax2.plot([5 * k for k in i[0]], [5 * k for k in i[1]], [5 * k for k in i[2]], c='black', linewidth=0.5)
    ax2.set_xlim(-5 * square_scale, 5 * square_scale)
    ax2.set_ylim(-5 * square_scale, 5 * square_scale)
    ax2.set_zlim(-5 * square_scale, 5 * square_scale)
    plt.savefig(output_3D + "/{}.png".format(time), format='png', dpi=200)

    fig_contour = plt.figure(fig_contour=(16, 9))
    ax2_countour = fig_contour.add_subplot(121)
    ax_countour = fig_contour.add_subplot(122)
    sc = ax_countour.tricontourf(
        eos_density,
        eos_internal_energy,
        eos_entropy,
        cmap=cmap,
        norm=normalizer,
        levels=50
    )
    for p in select_particles_at_time:
        ax_countour.scatter(
            [p.density],
            [p.internal_energy],
            marker="o",
            linewidths=1,
            # facecolor=[cmap(normalizer(p.entropy)) for p in select_particles],
            facecolor=(0, 0, 0, 0),
            edgecolors=c_dict[p.particle_id],
            label="All select_particles"
        )
    prev_particles_tmp = {}
    if time == start_time:
        for p in select_particles_at_time:
            prev_particles.update({p.particle_id: []})
    for p in select_particles_at_time:
        prev_particles[p.particle_id].append(p)
    if time > start_time:
        for p in select_particles_at_time:
            prev = prev_particles[p.particle_id]
            ax_countour.plot(
                [p2.density for p2 in prev],
                [p2.internal_energy for p2 in prev],
                c=c_dict[p.particle_id],
                linewidth=2.0
            )
            ax2_countour.plot(
                [p2.pressure / (10 ** 9) for p2 in prev],
                [p2.entropy for p2 in prev],
                c=c_dict[p.particle_id],
                linewidth=2.0
            )
    ax2_countour.set_xlabel("Pressure (GPa")
    ax2_countour.set_ylabel("Entropy")
    ax2_countour.grid(alpha=0.4)
    ax_countour.grid(alpha=0.4)
    ax_countour.set_xlim(-5, 2000)
    ax_countour.set_ylim(0, 8e7)
    ax2_countour.set_xlim(0, 10)
    ax2_countour.set_ylim(5000, 12000)
    ax_countour.set_xlabel("Density")
    ax_countour.set_ylabel("Internal Energy")
    ax_countour.set_title("Time: {} hrs (iteration: {})".format(round(seconds_to_hours(seconds), 2), time))
    ax2_countour.set_title("Time: {} hrs (iteration: {})".format(round(seconds_to_hours(seconds), 2), time))

    cbar = fig_contour.colorbar(sm)
    cbar.ax_countour.set_title("Entropy")

    plt.savefig(output + "/{}.png".format(time), format='png', dpi=200)
