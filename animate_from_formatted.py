#!/usr/bin/env python3
import os
import shutil
import csv
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.pyplot as plt

from src.animate import animate

plt.style.use("dark_background")

start_time = 1
end_time = 180
increment = 1
path = "/home/theia/scotthull/1M_high_rho_cutoff/formatted_gi_new_eos_b_073_high_rho_cutoff_1M"
to_path = "/home/theia/scotthull/FDPS_SPH_Orbits/animated_shots"
square_scale = 4e7
min_normalize_parameter = 1000
max_normalize_parameter = 6000

for p in [to_path]:
    if os.path.exists(p):
        shutil.rmtree(p)
    os.mkdir(p)

normalizer = Normalize(min_normalize_parameter, max_normalize_parameter)
cmap = cm.get_cmap('jet')

def get_time(f):
    formatted_time = None
    with open(f, 'r') as infile:
        reader = csv.reader(infile, delimiter="\t")
        formatted_time = float(next(reader)[0])
    infile.close()
    return round(formatted_time * 0.000277778, 2)  # seconds -> hours


for time in np.arange(start_time, end_time + increment, increment):
    print("At time {}".format(time))
    f = path + "/{}.csv".format(time)
    formatted_time = get_time(f)
    df = pd.read_csv(f, skiprows=2).to_dict('list')
    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(111)
    ax.scatter(
        [df['x'][index] for index, z in enumerate(df['z']) if z < 0],
        [df['y'][index] for index, z in enumerate(df['z']) if z < 0],
        c=[cmap(normalizer(df['entropy'][index])) for index, z in enumerate(df['z']) if z < 0],
        s=1
    )
    ax.set_xticks([])
    # for minor ticks
    ax.set_xticks([], minor=True)
    ax.set_yticks([])
    # for minor ticks
    ax.set_yticks([], minor=True)
    ax.set_xlim(-square_scale, square_scale)
    ax.set_ylim(-square_scale, square_scale)
    ax.set_title(
        str(formatted_time) + " hrs",
        c="white",
    )
    sm = cm.ScalarMappable(norm=normalizer, cmap=cmap)
    sm.set_array([])
    cbaxes = inset_axes(ax, width="30%", height="3%", loc=2, borderpad=1.8)
    cbar = plt.colorbar(sm, cax=cbaxes, orientation='horizontal')
    cbar.ax.tick_params(labelsize=6)
    cbar.ax.set_title("Entropy", fontsize=6)
    plt.savefig(to_path + "/{}.png".format(time), format='png')

animate(
    start_time=start_time,
    end_time=end_time,
    interval=increment,
    path=to_path,
    fps=20,
    filename="animated_from_formatted.mp4",
)

