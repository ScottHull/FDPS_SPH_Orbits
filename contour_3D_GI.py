#!/usr/bin/env python3
import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.cm as cm

from src.identify import ParticleMap
from src.combine import CombineFile
from src.animate import animate
from src.new_and_old_eos import seconds_to_hours

start_time = 0
end_time = 3000
interval = 20
number_processes = 200
min_norm = 0
max_norm = 10000
path = "/home/theia/scotthull/gi_new_eos"
output = "/home/theia/scotthull/FDPS_SPH_Orbits/3D_contour_GI"

normalizer = Normalize(min_norm, max_norm)
cmap = cmap = cm.get_cmap('jet')
sm = cm.ScalarMappable(norm=normalizer, cmap=cmap)

for o in [output]:
    if os.path.exists(o):
        shutil.rmtree(o)
    os.mkdir(o)

for time in np.arange(start_time, end_time + interval, interval):
    cf = CombineFile(num_processes=number_processes, time=time, output_path=path)
    combined_file = cf.combine()
    formatted_time = cf.sim_time
    f = os.getcwd() + "/merged_{}.dat".format(time)
    pm = ParticleMap(path=f, center=True, relative_velocity=False)
    particles = pm.collect_particles(find_orbital_elements=True)
    os.remove(f)

    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(111)
    density = np.array([p.density for p in particles])
    internal_energy = np.array([p.internal_energy for p in particles])
    entropy = [p.entropy for p in particles]
    sc = ax.tricontourf(
        density, internal_energy, entropy, cmap=cmap, norm=normalizer, levels=10
    )
    ax.scatter(
        density, internal_energy,
        marker="o",
        linewidths=0.2,
        c=[cmap(normalizer(p.entropy)) for p in particles],
        edgecolors='black'
    )
    ax.set_xlabel("Density")
    ax.set_ylabel("Internal Energy")
    ax.set_title("Time: {} hrs (iteration: {})".format(round(seconds_to_hours(formatted_time), 2), time))
    cbar = fig.colorbar(sm)
    cbar.ax.set_title("Entropy")

    plt.savefig(output + "/{}.png".format(time), format='png', dpi=200)
