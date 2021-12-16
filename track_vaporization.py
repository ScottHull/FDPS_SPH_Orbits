import os
import csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from src.new_and_old_eos import get_particles
from src.vapor import calc_vapor_mass_fraction

path = "/home/theia/scotthull/1M_high_rho_cutoff/formatted_gi_new_eos_b_073_high_rho_cutoff_1M"
to_path = "/home/theia/scotthull/FDPS_SPH_Orbits/track_vmf"
phase_path = "src/phase_data/forstSTS__vapour_curve.txt"
start_time = 100
end_time = 600
increment = 100

times, vmfs = [], []
for time in np.arange(start_time, end_time + increment, increment):
    f = path + "/{}.csv".format(time)

    particles, formatted_time = get_particles(path=path, formatted=True, time=time, number_processes=200)
    vmf = calc_vapor_mass_fraction(particles=particles, phase_path=phase_path, only_disk=True) * 100.0
    times.append(formatted_time)
    vmfs.append(vmf)

fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
ax.plot(
    times,
    vmfs,
    linewidth=2.0
)
ax.set_xlabel("Time (hrs)")
ax.set_ylabel("VMF (%)")
ax.grid(alpha=0.4)
plt.savefig("vmf.png", format='png')

