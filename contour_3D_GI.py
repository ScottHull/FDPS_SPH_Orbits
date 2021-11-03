#!/usr/bin/env python3
import os
import shutil
import numpy as np
import matplotlib.pyplot as plt

from src.identify import ParticleMap
from src.combine import CombineFile
from src.animate import animate
from src.new_and_old_eos import seconds_to_hours

start_time = 0
end_time = 3000
interval = 20
number_processes = 200
path = "/home/theia/scotthull/gi_new_eos"
output = "/home/theia/scotthull/FDPS_SPH_Orbits/3D_contour_GI"

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
    particles = pm.collect_particles(find_orbital_elements=False)
    os.remove(f)

    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(111)
    density = np.array([p.density for p in particles])
    internal_energy = np.array([p.internal_energy for p in particles])
    entropy = [p.entropy for p in particles]
    plt.tricontourf(
        density, internal_energy, entropy
    )
    plt.colorbar()

    plt.savefig(output + "/{}.png", format=time)
