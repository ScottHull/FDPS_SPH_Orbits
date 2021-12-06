#!/usr/bin/env python3
import csv
import pandas as pd
import numpy as np
from statistics import mean
import matplotlib.pyplot as plt

new_path = "/home/theia/scotthull/1M/gi_new_eos_b_073_at_time"
old_path = "/home/theia/scotthull/1M/gi_old_eos_b_073_at_time"
start_time = 0
end_time = 2500
increment = 100


def get_time(f):
    formatted_time = None
    with open(f, 'r') as infile:
        reader = csv.reader(infile, delimiter="\t")
        formatted_time = float(next(reader)[0])
    infile.close()
    return formatted_time * 0.000277778  # seconds -> hours

new_times, old_times, new_avg_entropies, old_avg_entropies = [], [], [], []
num_particles_disk_new, num_particles_disk_old = [], []
for time in np.arange(start_time, end_time + increment, increment):
    print("At time : {}".format(time))
    new_f, old_f = new_path + "/{}.csv".format(time), old_path + "/{}.csv".format(time)
    new_time, old_time = get_time(new_f), get_time(old_f)
    new_times.append(new_time), old_times.append(old_time)
    new_file = pd.read_csv(new_f, skiprows=2).to_dict()
    old_file = pd.read_csv(old_f, skiprows=2).to_dict()
    new_disk_entropies = [s for index, s in enumerate(new_file['entropy']) if new_file['label'][index] == "DISK"]
    old_disk_entropies = [s for index, s in enumerate(old_file['entropy']) if old_file['label'][index] == "DISK"]
    if len(new_disk_entropies) > 0:
        new_avg_entropies.append(mean(new_disk_entropies))
    else:
        new_avg_entropies.append(0)
    if len(old_disk_entropies) > 0:
        old_avg_entropies.append(mean(old_disk_entropies))
    else:
        old_avg_entropies.append(0)

    num_particles_disk_new.append(len(new_disk_entropies))
    num_particles_disk_old.append(len(old_disk_entropies))

fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
ax.plot(
    new_times, new_avg_entropies, linewidth=2.0, label="New EoS"
)
ax.plot(
    old_times, old_avg_entropies, linewidth=2.0, label="Old EoS"
)
ax.set_xlabel("Time (hrs)"), ax.set_ylabel("Avg. Disk Entropy")
ax.grid(alpha=0.4)
ax.legend()
plt.savefig("avg_entropies_from_formatted.png", format='png')

fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
ax.plot(
    new_times, num_particles_disk_new, linewidth=2.0, label="New EoS"
)
ax.plot(
    old_times, num_particles_disk_old, linewidth=2.0, label="Old EoS"
)
ax.set_xlabel("Time (hrs)"), ax.set_ylabel("Num. Disk Particles")
ax.grid(alpha=0.4)
ax.legend()
plt.savefig("num_disk_particles.png", format='png')
    