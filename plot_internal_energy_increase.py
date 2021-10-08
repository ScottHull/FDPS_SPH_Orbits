import os
import shutil
from math import sqrt
import numpy as np
import matplotlib.pyplot as plt

from src.identify import ParticleMap
from src.combine import CombineFile

start_time = 0
end_time = 3000
interval = 50
number_processes = 100
path = "/home/theia/scotthull/sph_simulations/gi_old_eos"

total_internal_energy = []
total_silicate_internal_energy = []
total_iron_internal_energy = []

for time in np.arange(0, end_time + interval, interval):
    cf = CombineFile(num_processes=number_processes, time=time, output_path=path)
    formatted_time = cf.sim_time
    combined_file = cf.combine()
    f = os.getcwd() + "/merged_{}.dat".format(time)
    pm = ParticleMap(path=f, center=True, relative_velocity=False)
    particles = pm.collect_particles()
    # pm.solve(particles=particles)
    os.remove(f)

    total_u = sum([p.internal_energy for p in particles])
    total_u_silicate = sum([p.internal_energy for p in particles if p.tag % 2 == 0])
    total_u_iron = sum([p.internal_energy for p in particles if p.tag % 2 != 0])
    total_internal_energy.append(total_u)
    total_silicate_internal_energy.append(total_u_silicate)
    total_iron_internal_energy.append(total_u_iron)

fig, ax = plt.subplots(3, 1, figsize=(16, 8))
ax = ax.flatten()
ax[0].plot(
    np.arange(0, end_time + interval, interval),
    total_internal_energy,
    linewidth=2.0,
    color="black"
)
ax[0].set_ylabel("Total Internal Energy")
ax[1].plot(
    np.arange(0, end_time + interval, interval),
    total_silicate_internal_energy,
    linewidth=2.0,
    color="black"
)
ax[1].set_ylabel("Total Silicate Energy")
ax[2].plot(
    np.arange(0, end_time + interval, interval),
    total_iron_internal_energy,
    linewidth=2.0,
    color="black"
)
ax[2].set_ylabel("Total Iron Energy")
ax[2].set_xlabel("Iteration")

plt.rcParams['axes.grid'] = True

plt.suptitle("GI (Old EoS)")

plt.savefig("internal_energy_increase.png", format='png')

