#!/usr/bin/env python3
import os
import shutil
from random import randint
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.font_manager as fm
from scipy.interpolate import UnivariateSpline

from src.identify import ParticleMap
from src.combine import CombineFile

time = 3000
number_processes = 200
path = "/home/theia/scotthull/200k/formatted_gi_new_eos/3000.csv"
radius_earth = 6371 * 1000

df = pd.read_csv(path, delimiter=",", skiprows=2)
disk_df = df.loc[df['label'] == "DISK"]
disk_df_target_silicate = disk_df.loc[df['tag'] == 0]
disk_df_target_iron = disk_df.loc[df['tag'] == 1]
disk_df_impactor_silicate = disk_df.loc[df['tag'] == 2]
disk_df_impactor_iron = disk_df.loc[df['tag'] == 3]
disk_dfs = [disk_df_target_silicate, disk_df_target_iron, disk_df_impactor_silicate, disk_df_impactor_iron]
tags = [0, 1, 2, 3]
labels = {
    0: "Target Silicate",
    1: "Target Iron",
    2: "Impactor Silicate",
    3: "Impactor Iron"
}

plt.style.use("dark_background")
fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
fig.patch.set_facecolor('xkcd:black')
# s = UnivariateSpline(
#     radius,
#     temperature,
#     s=5
# )
# avg_temp = s(radius)
for index, d in enumerate(disk_dfs):
    ax.scatter(
        [i / radius_earth for i in d['radius']],
        d['temperature'],
        s=1,
        label=labels[tags[index]]
    )
ax.grid(alpha=0.4)
ax.set_xlabel(r"Radius ($R_{\bigoplus}$)")
ax.set_ylabel("Disk Temperature")
ax.set_title("iteration {}".format(time))
ax.legend()
plt.savefig("disk_temp_profile.png", format='png', dpi=400)

fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
fig.patch.set_facecolor('xkcd:black')
for index, d in enumerate(disk_dfs):
    ax.scatter(
        [i / radius_earth for i in d['radius']],
        d['pressure'],
        s=1,
        label=labels[tags[index]]
    )
ax.grid(alpha=0.4)
ax.set_xlabel(r"Radius ($R_{\bigoplus}$)")
ax.set_ylabel("Disk Pressure")
ax.set_title("iteration {}".format(time))
ax.legend()
plt.savefig("disk_pres_profile.png", format='png', dpi=400)
