import os
import numpy as np
import csv
import matplotlib.pyplot as plt
from src.report import BuildReports

"""
This code produces cumulative outputs, but only "solves" for the disk at certain time intervals.
Useful for getting a coarse understanding of disk mass change without having to solve at each timestep.
"""

plt.style.use("dark_background")

path = "/home/theia/scotthull/1M_high_rho_cutoff/gi_new_eos_b_073_high_rho_cutoff_1M"
num_processes = 200
assess_at_interval = 100
start_iteration = 0
end_iteration = 2900
end_disk_state_iteration = 1000

def get_time(f):
    formatted_time = None
    with open(f, 'r') as infile:
        reader = csv.reader(infile, delimiter="\t")
        formatted_time = float(next(reader)[0])
    infile.close()
    return round(formatted_time * 0.000277778, 2)  # seconds -> hours

times = []
for i in np.arange(0, 3000 + 1, 1):
    fname = "results.{}_{}_{}.dat".format(str(i).zfill(5), str(num_processes).zfill(5), str(1).zfill(5))
    f = path + "/{}.dat".format(fname)
    time = get_time(f)
    times.append(time)

fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
ax.plot(
    np.arange(0, 3000 + 1, 1),
    times,
    linewidth=2.0
)
ax.axhline(30, color='red', linestyle="--")
ax.set_xlabel("Iteration")
ax.set_ylabel("Time (hrs)")

plt.savefig("iteration_vs_time.png", format='png')
