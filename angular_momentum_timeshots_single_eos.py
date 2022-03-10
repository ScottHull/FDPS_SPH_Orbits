#!/usr/bin/env python3
import os
import csv
import shutil
from math import pi, asin, isnan
import numpy as np
import pandas as pd
from random import randint
from statistics import mean
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import multiprocessing as mp

from src.vapor import calc_vapor_mass_fraction_from_formatted
from src.geometry import get_impact_geometry_from_formatted, get_velocity_profile_from_formatted
from src.animate import animate
from src.identify import ParticleMap
from src.combine import CombineFile

base_path = "/home/theia/scotthull/Paper1_SPH/gi/"
runs = "new"
angle = "b073"
iterations = [50, 100, 200, 300]

new_phase_path = "src/phase_data/forstSTS__vapour_curve.txt"
old_phase_path = "src/phase_data/duniteN__vapour_curve.txt"

phase_path = new_phase_path
if runs == "old":
    phase_path = old_phase_path

run_set = ["{}_{}_{}".format(i, angle, runs) for i in [5, 500, 1000, 2000]]
paths = ["formatted_{}/{}".format(i, i) for i in run_set]

fig, axs = plt.subplots(len(paths), len(iterations), figsize=(16, 42), sharey='all',
                            gridspec_kw={"hspace": 0.10, "wspace": 0.10})

current_index = 0
for iteration in iterations:
    for p in paths:
        df = pd.read_csv(p + "/{}.dat".format(iteration), skiprows=2)
        positions = list(zip(df['x'], df['y'], df['z']))
        velocities = list(zip(df['vx'], df['vy'], df['vz']))
        spec_am = [np.linalg.norm(np.cross(i, j)) for i, j in zip(positions, velocities)]
        

        current_index += 1
