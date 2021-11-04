#!/usr/bin/env python3
import os
import shutil
from random import randint
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.cm as cm

from src.identify import ParticleMap
from src.combine import CombineFile
from src.animate import animate
from src.new_and_old_eos import seconds_to_hours

start_time = 80
end_time = 3000
interval = 20
number_processes = 200
min_norm = 0
max_norm = 10000
path = "/home/theia/scotthull/gi_new_eos"
eos = "src/phase_data/forst_STS.rho_u.txt"
output = "/home/theia/scotthull/FDPS_SPH_Orbits/3D_contour_GI"

normalizer = Normalize(min_norm, max_norm)
cmap = cm.get_cmap('jet')
sm = cm.ScalarMappable(norm=normalizer, cmap=cmap)

for o in [output]:
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
particles = pm_end.collect_particles()
pm_end.solve(particles=particles, phase_path="src/phase_data/forstSTS__vapour_curve.txt")
os.remove(f)

end = {}
high_entropy = {}
for p in particles:
    end.update({p.particle_id: p.label})
    if p.entropy > 8000 and p.label == "DISK":
        high_entropy.update({p.particle_id: p.entropy})
high_entropy_ids = list(high_entropy.keys())
rand_select = [high_entropy_ids[randint(0, len(high_entropy_ids) - 1)] for i in range(0, 5)]
colors = ["black", "red", "blue", "green", "white"]
c_dict = {}
for index, i in enumerate(rand_select):
    c_dict.update({i: colors[index]})

prev_particles = {}

for time in np.arange(start_time, end_time + interval, interval):
    cf = CombineFile(num_processes=number_processes, time=time, output_path=path)
    combined_file = cf.combine()
    formatted_time = cf.sim_time
    f = os.getcwd() + "/merged_{}.dat".format(time)
    pm = ParticleMap(path=f, center=True, relative_velocity=False)
    particles = pm.collect_particles(find_orbital_elements=False)
    os.remove(f)
    particles = [p for p in particles if p.particle_id in rand_select]

    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(111)
    sc = ax.tricontourf(
        eos_density,
        eos_internal_energy,
        eos_entropy,
        cmap=cmap,
        norm=normalizer,
        levels=20
    )
    for p in particles:
        ax.scatter(
            [p.density],
            [p.internal_energy],
            marker="o",
            linewidths=1,
            # facecolor=[cmap(normalizer(p.entropy)) for p in particles],
            facecolor=(0, 0, 0, 0),
            edgecolors=c_dict[p.particle_id],
            label="All Particles"
        )
    prev_particles_tmp = {}
    if time == start_time:
        for p in particles:
            prev_particles.update({p.particle_id: []})
    for p in particles:
        prev_particles[p.particle_id].append(p)
    if time > start_time:
        for p in particles:
            prev = prev_particles[p.particle_id]
            ax.plot(
                [p2.density for p2 in prev],
                [p2.internal_energy for p2 in prev],
                c=c_dict[p.particle_id],
                linewidth=2.0
            )
    ax.grid(alpha=0.4)
    ax.set_xlim(-5, 2000)
    ax.set_ylim(0, 8e7)
    ax.set_xlabel("Density")
    ax.set_ylabel("Internal Energy")
    ax.set_title("Time: {} hrs (iteration: {})".format(round(seconds_to_hours(formatted_time), 2), time))
    cbar = fig.colorbar(sm)
    cbar.ax.set_title("Entropy")

    plt.savefig(output + "/{}.png".format(time), format='png', dpi=200)
