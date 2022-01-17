import os
import csv
from statistics import mean
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from src.new_and_old_eos import get_particles
from src.vapor import calc_vapor_mass_fraction

plt.style.use("dark_background")

path = "/home/theia/scotthull/1M_high_rho_cutoff/formatted_gi_new_eos_b_073_high_rho_cutoff_1M"
to_path = "/home/theia/scotthull/FDPS_SPH_Orbits/track_vmf"
phase_path = "src/phase_data/forstSTS__vapour_curve.txt"
start_time = 100
end_time = 800
increment = 100

def seconds_to_hours(seconds):
    return round(seconds * 0.000277778, 2)

times, vmfs, mean_entropy, num_disk_particles = [], [], [], []
for time in np.arange(start_time, end_time + increment, increment):
    f = path + "/{}.csv".format(time)

    particles, formatted_time = get_particles(path=path, formatted=True, time=time, number_processes=200)
    formatted_time = seconds_to_hours(formatted_time)
    mean_disk_entropy = mean([p.entropy for p in particles if p.label == "DISK"])
    number_disk_particles = len([p.entropy for p in particles if p.label == "DISK"])

    vmf = calc_vapor_mass_fraction(particles=particles, phase_path=phase_path, only_disk=True) * 100.0
    times.append(formatted_time)
    vmfs.append(vmf)
    mean_entropy.append(mean_disk_entropy)
    num_disk_particles.append(number_disk_particles)

fig, axs = plt.subplots(1, 3, figsize=(16, 9), sharex='all', gridspec_kw={"hspace": 0.10, "wspace": 0.14})
ax1, ax2, ax3 = axs.flatten()
ax1.plot(
    times,
    vmfs,
    linewidth=2.0
)
ax1.set_ylabel("VMF (%)")

ax2.plot(
    times,
    mean_entropy,
    linewidth=2.0
)
ax2.set_ylabel("Avg. Disk Entropy")

ax3.plot(
    times,
    num_disk_particles,
    linewidth=2.0
)
ax3.set_ylabel("Number Disk Particles")

for ax in axs.flatten():
    ax.set_xlabel("Time (hrs)")
    ax.grid(alpha=0.4)

plt.savefig("vmf.png", format='png')

