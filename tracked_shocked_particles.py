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

start_time = 0
end_time = 3000
increment = 100
start_shock_sample = 400
end_shock_sample = 2500
f_path = "/home/theia/scotthull/1M/gi_new_eos_b_073_at_time"
output = "/home/theia/scotthull/FDPS_SPH_Orbits/track_high_s_shocks"

for p in [output]:
    if os.path.exists(p):
        shutil.rmtree(p)
    os.mkdir(p)


def get_time(f):
    formatted_time = None
    with open(f, 'r') as infile:
        reader = csv.reader(infile, delimiter="\t")
        formatted_time = float(next(reader)[0])
    infile.close()
    return formatted_time * 0.000277778  # seconds -> hours


max_time = get_time(f_path + "/{}.csv".format(end_time))


def identify_shocked_particles(start_sample_time, end_sample_time, num=5, s_cutoff=8000):
    end_f = f_path + "/{}.csv".format(end_sample_time)
    start_f = f_path + "/{}.csv".format(start_sample_time)
    end_df = pd.read_csv(end_f, skiprows=2).to_dict('list')
    start_df = pd.read_csv(start_f, skiprows=2).to_dict('list')
    ids = []
    high_s_ids = [i for index, i in enumerate(end_df['id']) if
                  end_df['label'][index] == "DISK" and end_df['entropy'][index] >= s_cutoff]
    for index, i in enumerate(start_df['id']):
        if len(ids) >= num:
            break
        if start_df['entropy'][index] < s_cutoff and start_df['label'][index] == "DISK" and start_df['id'][index] in high_s_ids:
            ids.append(i)
    return ids


ids = identify_shocked_particles(start_sample_time=start_shock_sample, end_sample_time=end_shock_sample)

plt.style.use("dark_background")
times = []
for time in np.arange(start_time, end_time + increment, increment):
    new_f = f_path + "/{}.csv".format(time)
    new_time = get_time(new_f)
    times.append(new_time)
    new_file = pd.read_csv(new_f, skiprows=2).to_dict('list')
    d = [(i, new_file['entropy'][index], new_file['density'][index], new_file['internal_energy'][index]) for index, i in
         enumerate(new_file['id']) if i in ids]
    fig, axs = plt.subplots(1, 3, figsize=(16, 9), gridspec_kw={"hspace": 0.0, "wspace": 0.14})
    fig.patch.set_facecolor('xkcd:black')
    ax1, ax2, ax3 = axs.flatten()
    for ax in axs.flatten():
        ax.set_xlabel("Time (hrs)")
        ax.grid(alpha=0.4)
    for i in d:
        ax1.plot(
            times,
            i[1],
            linewidth=2.0,
            label=i[0]
        )
        ax2.plot(
            times,
            i[2],
            linewidth=2.0,
            label=i[0]
        )
        ax3.plot(
            times,
            i[3],
            linewidth=2.0,
            label=i[0]
        )
    ax1.set_ylabel("Entropy")
    ax2.set_ylabel("Density")
    ax3.set_ylabel("Internal Energy")

    plt.savefig(output + "/{}.png".format(time), format='png')

animate(
    start_time=start_time,
    end_time=end_time,
    interval=increment,
    path=output,
    fps=10,
    filename="disk_s_shocks.mp4",
)

