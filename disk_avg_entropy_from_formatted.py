#!/usr/bin/env python3
import pandas as pd
import numpy as np
from statistics import mean
import matplotlib.pyplot as plt

new_path = "/home/theia/scotthull/1M/gi_new_eos_b_073_at_time"
old_path = "/home/theia/scotthull/1M/gi_old_eos_b_073_at_time"
start_time = 0
end_time = 2500
increment = 100

times, new_avg_entropies, old_avg_entropies = [], [], []
num_particles_disk_new, num_particles_disk_old = [], []
for time in np.arange(start_time, end_time + increment, increment):
    print("At time : {}".format(time))
    times.append(time)
    new_file = pd.read_csv(new_path + "/{}.csv".format(time), skiprows=2).to_dict()
    old_file = pd.read_csv(old_path + "/{}.csv".format(time), skiprows=2).to_dict()
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
    times, new_avg_entropies, linewidth=2.0, label="New EoS"
)
ax.plot(
    times, old_avg_entropies, linewidth=2.0, label="Old EoS"
)
ax.set_xlabel("Iteration"), ax.set_ylabel("Avg. Entropy")
ax.grid(alpha=0.4)
ax.legend()
plt.savefig("avg_entropies_from_formatted.png", format='png')

fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
ax.plot(
    times, num_particles_disk_new, linewidth=2.0, label="New EoS"
)
ax.plot(
    times, num_particles_disk_old, linewidth=2.0, label="Old EoS"
)
ax.set_xlabel("Iteration"), ax.set_ylabel("Num. Disk Particles")
ax.grid(alpha=0.4)
ax.legend()
plt.savefig("num_disk_particles.png", format='png')
    