#!/usr/bin/env python3
import os
import shutil
import csv
import pandas as pd
import numpy as np
from statistics import mean
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from src.animate import animate

new_path = "/home/theia/scotthull/1M/gi_new_eos_b_073_at_time"
old_path = "/home/theia/scotthull/1M/gi_old_eos_b_073_at_time"
path = "/home/theia/scotthull/FDPS_SPH_Orbits/disk_avgs"
start_time = 0
end_time = 3000
increment = 100

if os.path.exists(path):
    shutil.rmtree(path)
os.mkdir(path)

def get_time(f):
    formatted_time = None
    with open(f, 'r') as infile:
        reader = csv.reader(infile, delimiter="\t")
        formatted_time = float(next(reader)[0])
    infile.close()
    return formatted_time * 0.000277778  # seconds -> hours


def get_particles_by_label(particles):
    return {
        "PLANET": [zip(x, particles['y'][index], particles['z'][index]) for index, x in enumerate(particles['x']) if
                   particles['label'] == "PLANET"],
        "DISK": [zip(x, particles['y'][index], particles['z'][index]) for index, x in enumerate(particles['x']) if
                 particles['label'] == "DISK"],
        "ESCAPE": [zip(x, particles['y'][index], particles['z'][index]) for index, x in enumerate(particles['x']) if
                   particles['label'] == "ESCAPE"]
    }


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

    labeled_particles_new = get_particles_by_label(new_file)
    labeled_particles_old =  get_particles_by_label(old_file)

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

    gs = gridspec.GridSpec(2, 2)
    fig = plt.figure(figsize=(10, 10))
    fig.patch.set_facecolor('xkcd:black')
    ax1 = fig.add_subplot(gs[0, 0])  # row 0, col 0
    ax2 = fig.add_subplot(gs[0, 1])  # row 0, col 1
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])
    axs = [ax1, ax2, ax3, ax4]
    for ax in axs:
        ax.grid(alpha=0.4)
    for label in labeled_particles_new.keys():
        ax1.scatter(
            [i[0] for i in labeled_particles_new[label] if i[2] < 0],
            [i[1] for i in labeled_particles_new[label] if i[2] < 0],
            s=1,
            label=label
        )
    ax1.legend(loc='upper right')
    for label in labeled_particles_old.keys():
        ax2.scatter(
            [i[0] for i in labeled_particles_old[label] if i[2] < 0],
            [i[1] for i in labeled_particles_old[label] if i[2] < 0],
            s=1,
            label=label
        )
    ax2.legend(loc='upper right')
    ax3.plot(
        new_times, new_avg_entropies, linewidth=2.0, label="New EoS"
    )
    ax3.plot(
        old_times, old_avg_entropies, linewidth=2.0, label="Old EoS"
    )
    ax3.set_xlabel("Time (hrs)")
    ax3.set_ylabel("Avg. Disk Entropy")
    ax3.legend(loc='upper right')
    ax4.plot(
        new_times, num_particles_disk_new, linewidth=2.0, label="New EoS"
    )
    ax4.plot(
        old_times, num_particles_disk_old, linewidth=2.0, label="Old EoS"
    )
    ax4.set_xlabel("Time (hrs)")
    ax4.set_ylabel("Num. Disk Particles")
    ax4.legend(loc='upper right')

    plt.savefig(path + "/{}.png".format(time), format='png')

animate(
    start_time=start_time,
    end_time=end_time,
    interval=increment,
    path=path,
    fps=10,
    filename="disk_avg.mp4",
)
    
