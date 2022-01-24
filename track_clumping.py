#!/usr/bin/env python3
import csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
plt.style.use("dark_background")

base_path = "/home/theia/scotthull/Paper1_SPH/gi/"
base_path_setups = "/home/theia/scotthull/Paper1_SPH/setups/"

gi_b073_runs = {
    "5_b073_new": {
        "name": "5b073n",
        "path": base_path + "5_b073_new/formatted_5_b073_new",
    },
    "5_b073_old": {
        "name": "5b073o",
        "path": base_path + "5_b073_old/formatted_5_b073_old",
    },
    "500_b073_new": {
        "name": "500b073n",
        "path": base_path + "500_b073_new/formatted_500_b073_new",
    },
    "500_b073_old": {
        "name": "500b073o",
        "path": base_path + "500_b073_old/formatted_500_b073_old",
    },
    "1000_b073_new": {
        "name": "1000b073n",
        "path": base_path + "1000_b073_new/formatted_1000_b073_new",
    },
    "1000_b073_old": {
        "name": "1000b073o",
        "path": base_path + "1000_b073_old/formatted_1000_b073_old",
    },
    "2000_b073_new": {
        "name": "2000b073n",
        "path": base_path + "2000_b073_new/formatted_2000_b073_new",
    },
    "2000_b073_old": {
        "name": "2000b073o",
        "path": base_path + "2000_b073_old/formatted_2000_b073_old",
    },
}

gi_b075_runs = {
    "5_b075_new": {
        "name": "5b075n",
        "path": base_path + "5_b075_new/formatted_5_b075_new",
    },
    "5_b075_old": {
        "name": "5b075o",
        "path": base_path + "5_b075_old/formatted_5_b075_old",
    },
    "500_b075_new": {
        "name": "500b075n",
        "path": base_path + "500_b075_new/formatted_500_b075_new",
    },
    "500_b075_old": {
        "name": "500b075o",
        "path": base_path + "500_b075_old/formatted_500_b075_old",
    },
    "1000_b075_new": {
        "name": "1000b075n",
        "path": base_path + "1000_b075_new/formatted_1000_b075_new",
    },
    "1000_b075_old": {
        "name": "1000b075o",
        "path": base_path + "1000_b075_old/formatted_1000_b075_old",
    },
    "2000_b075_new": {
        "name": "2000b075n",
        "path": base_path + "2000_b075_new/formatted_2000_b075_new",
    },
    "2000_b075_old": {
        "name": "2000b075o",
        "path": base_path + "2000_b075_old/formatted_2000_b075_old",
    },
}

def get_time(f):
    formatted_time = None
    with open(f, 'r') as infile:
        reader = csv.reader(infile, delimiter="\t")
        formatted_time = float(next(reader)[0])
    infile.close()
    return round(formatted_time * 0.000277778, 2)  # seconds -> hours

seleted = gi_b073_runs
start_iteration = 20
end_iteration = 200
want = "new"
num_cols = 4
square_scale = 4e7
increment = (end_iteration - start_iteration) / num_cols
L_EM = 3.5 * 10 ** 34
normalizer = Normalize(0, 2)
cmap = cm.get_cmap('jet')

fig, axs = plt.subplots(4, num_cols, figsize=(32, 32), sharex='all',
                            gridspec_kw={"hspace": 0.0, "wspace": 0.00})
axs = axs.flatten()
for ax in axs:
    ax.set_xlim(-square_scale, square_scale)
    ax.set_ylim(-square_scale, square_scale)
    ax.set_xticks([], minor=False)
    ax.set_yticks([], minor=False)
for i in seleted.keys():
    ax_index = 0
    if ax_index == 0:
        sm = cm.ScalarMappable(norm=normalizer, cmap=cmap)
        sm.set_array([])
        cbaxes = inset_axes(axs[0], width="30%", height="3%", loc=2, borderpad=1.8)
        cbar = plt.colorbar(sm, cax=cbaxes, orientation='horizontal')
        cbar.ax.tick_params(labelsize=6)
        cbar.ax.set_title("Ang. Momentum ($L_{EM}$", fontsize=6)
    if want in i:
        n, p = seleted[i]['name'], seleted[i]['path']
        times = []
        for time in np.arange(start_iteration, end_iteration + 1, 1):
            if time % increment == 0:
                print("at iteration {} (ax index {})".format(time, ax_index))
                t = get_time(p + "/{}.csv".format(time))
                times.append(t)
                df = pd.read_csv(p + "/{}.csv".format(time), skiprows=2)
                df = df[df['z'] < 0]
                masses = list(df['mass'])
                positions = zip(df['x'], df['y'], df['z'])
                velocities = list(zip(df['vx'], df['vy'], df['vz']))
                angular_momenta = [masses[index] * np.linalg.norm(np.cross(p, velocities[index])) / L_EM for index, p in
                                          enumerate(positions)]
                axs[ax_index].scatter(
                    df['x'], df['y'], s=2,
                    color=[cmap(normalizer(masses[index] * np.linalg.norm(np.cross(p, velocities[index])) / L_EM)) for index, p in
                                          enumerate(positions)]
                )
                axs[ax_index].annotate(n + "\n{} hrs".format(t), (square_scale - (square_scale * 0.15), square_scale - (square_scale * 0.15)), fontsize=8)
                ax_index += 1

plt.savefig("{}_track_clumping.png".format(want), format='png', dpi=200)

            




