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
radius, temperature, pressure = [i / radius_earth for i in disk_df['distance']], disk_df['temperature'], disk_df['pressure']

plt.style.use("dark_background")
fig = plt.figure()
ax = fig.add_subplot(111)
fig.patch.set_facecolor('xkcd:black')


s = UnivariateSpline(
    radius,
    temperature,
    s=5
)
avg_temp = s(radius)

ax.scatter(
    [p.distance / radius_earth for p in disk if p == 0],
    [p.temperature for p in disk if p == 0],
    s=0.2,
    label="Target Silicate"
)
ax.scatter(
    [p.distance / radius_earth for p in disk if p == 1],
    [p.temperature for p in disk if p == 1],
    s=0.2,
    label="Target Iron"
)
ax.scatter(
    [p.distance / radius_earth for p in disk if p == 2],
    [p.temperature for p in disk if p == 2],
    s=0.2,
    label="Impactor Silicate"
)
ax.scatter(
    [p.distance / radius_earth for p in disk if p == 3],
    [p.temperature for p in disk if p == 3],
    s=0.2,
    label="Impactor Iron"
)
ax.plot(
    radius,
    avg_temp,
    linewidth=1.0,
    label="Smoothed Temperature"
)
ax.grid(alpha=0.4)
ax.set_xlabel(r"Radius ($R_{\bigoplus}$)")
ax.set_ylabel("Disk Temperature")
ax.set_title("{} hrs (iteration {})".format(formatted_time, time))
ax.legend()
plt.savefig("disk_temp_profile.png", format='png')
