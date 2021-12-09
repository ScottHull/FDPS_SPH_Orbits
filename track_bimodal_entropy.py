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
path = "/home/theia/scotthull/FDPS_SPH_Orbits/disk_s_tracking"
start_time = 0
end_time = 3000
increment = 100

for p in [path]:
    if os.path.exists(p):
        shutil.rmtree(p)
    os.mkdir(p)

plt.style.use("dark_background")

def get_time(f):
    formatted_time = None
    with open(f, 'r') as infile:
        reader = csv.reader(infile, delimiter="\t")
        formatted_time = float(next(reader)[0])
    infile.close()
    return formatted_time * 0.000277778  # seconds -> hours

max_time = get_time(new_path + "/{}.csv".format(end_time))

def get_particles(particles):
    return [
        (s, tag, label, x, y, z) for s, tag, label, x, y, z in
        zip(particles['entropy'], particles['tag'], particles['label'], particles['x'], particles['y'], particles['z'])
        if label == "DISK"
    ]


num_high_s_particles = []
num_low_s_particles = []
times = []
for time in np.arange(start_time, end_time + increment, increment):
    new_f = new_path + "/{}.csv".format(time)
    new_time = get_time(new_f)
    times.append(new_time)
    new_file = pd.read_csv(new_f, skiprows=2).to_dict('list')
    disk_particles = get_particles(new_file)
    high_s = [i for i in disk_particles if i[0] >= 8000]
    low_s = [i for i in disk_particles if i[0] < 8000]
    num_low_s_particles.append(len(low_s))
    num_high_s_particles.append(len(high_s))
    fig, axs = plt.subplots(1, 2, figsize=(16, 9), gridspec_kw={"hspace": 0.0, "wspace": 0.14})
    fig.patch.set_facecolor('xkcd:black')
    ax1, ax2 = axs.flatten()

    ax1.scatter(
        [i[3] for i in low_s],
        [i[4] for i in low_s],
        color='aqua',
        s=3,
        label="S < 8000"
    )

    ax1.scatter(
        [i[3] for i in high_s],
        [i[4] for i in high_s],
        color='pink',
        s=3,
        label="S >= 8000"
    )

    ax1.set_xticks([])
    # for minor ticks
    ax1.set_xticks([], minor=True)
    ax1.set_yticks([])
    # for minor ticks
    ax1.set_yticks([], minor=True)

    ax2.plot(
        times,
        num_high_s_particles,
        linewidth=2.0,
        label="S >= 8000"
    )
    ax2.plot(
        times,
        num_low_s_particles,
        linewidth=2.0,
        label="S < 8000"
    )

    ax2.set_xlim(0, max_time)
    ax2.set_xlabel("Time (hrs)")
    ax2.set_ylabel("Num. Disk Particles")
    for ax in axs.flatten():
        ax.grid(alpha=0.4)
        ax.legend(loc='upper right')

    plt.savefig(path + "/{}.png".format(time), format='png')

animate(
    start_time=start_time,
    end_time=end_time,
    interval=increment,
    path=path,
    fps=10,
    filename="disk_s_bimodal.mp4",
)
