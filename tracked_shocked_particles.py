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
increment = 2
start_shock_sample = 200
end_shock_sample = 3000
f_path = "/home/theia/scotthull/Paper1_SPH/gi/5_b073_new/formatted_5_b073_new"
f_fine_path = "/home/theia/scotthull/Paper1_SPH/gi/5_b073_new/formatted_5_b073_new"
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
colors = ['r', 'g', 'b', 'magenta', 'w']
times = []
color_pairs = dict(zip(ids, colors))
data = dict(zip(ids, [{"s": [], "rho": [], "u": []} for i in ids]))
for time in np.arange(start_time, end_time + increment, increment):
    print('At time {}'.format(time))
    new_f = f_fine_path + "/{}.csv".format(time)
    new_time = get_time(new_f)
    times.append(new_time)
    new_file = pd.read_csv(new_f, skiprows=2, index_col='id')
    d = [[index, new_file['entropy'][index], new_file['density'][index], new_file['internal_energy'][index]] for index in ids]
    for i in d:
        data[i[0]]['s'].append(i[1])
        data[i[0]]['rho'].append(i[2])
        data[i[0]]['u'].append(i[3])
    fig, axs = plt.subplots(3, 1, figsize=(9, 16), gridspec_kw={"hspace": 0.12, "wspace": 0})
    fig.patch.set_facecolor('xkcd:black')
    ax1, ax2, ax3 = axs.flatten()
    for ax in axs.flatten():
        ax.set_xlim(0, max_time)
        # ax.set_xlabel("Time (hrs)")
        ax.grid(alpha=0.4)
    for index, i in enumerate(d):
        ax1.plot(
            times,
            data[i[0]]['s'],
            linewidth=2.0,
            label=i[0],
            c=color_pairs[i[0]],
            alpha=0.6
        )
        ax2.plot(
            times,
            data[i[0]]['rho'],
            linewidth=2.0,
            label=i[0],
            c=color_pairs[i[0]],
            alpha=0.6
        )
        ax3.plot(
            times,
            data[i[0]]['u'],
            linewidth=2.0,
            label=i[0],
            c=color_pairs[i[0]],
            alpha=0.6
        )
    ax1.set_ylim(0, 10000)
    ax2.set_ylim(0, 20),
    ax3.set_ylim(0, 4e7)
    ax3.set_xlabel("Time (hrs)")
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

