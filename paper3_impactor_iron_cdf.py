#!/usr/bin/env python3
import os
import csv
import shutil
from math import pi, asin, isnan, exp
import numpy as np
import pandas as pd
from random import randint
from statistics import mean
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LogNorm
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import multiprocessing as mp
from matplotlib.ticker import MaxNLocator
import string

from src.vapor import calc_vapor_mass_fraction_from_formatted
from src.geometry import get_impact_geometry_from_formatted, get_velocity_profile_from_formatted
from src.animate import animate
from src.identify import ParticleMap
from src.combine import CombineFile

base_path = "/home/theia/scotthull/Paper3_SPH/gi/500_mars_b073_2v_esc/500_mars_b073_2v_esc"
num_processes = 600
iteration = 100
mars_radius = 3389.5

headers = ["id", "tag", "mass", "x", "y", "z", "vx", "vy", "vz", "density", "internal energy", "pressure",
                   "potential energy", "entropy", "temperature"]

def center_of_mass(x: np.array, y: np.array, z: np.array, mass: np.array):
    """
    Calculate the center of mass of a system of particles
    :param x:
    :param y:
    :param z:
    :param mass:
    :return:
    """
    # calculate the center of mass
    x_com = np.sum(x * mass) / np.sum(mass)
    y_com = np.sum(y * mass) / np.sum(mass)
    z_com = np.sum(z * mass) / np.sum(mass)
    return x_com, y_com, z_com

to_fname = "merged_{}_{}.dat".format(iteration, randint(0, 100000))
cf = CombineFile(num_processes=num_processes, time=iteration, output_path=base_path, to_fname=to_fname)
df = cf.combine_df()
df.columns = headers
# calculate the martian COM

df['radius'] = np.sqrt(df['x'] ** 2 + df['y'] ** 2 + df['z'] ** 2) / 1000
particles_within_mars = df[df['radius'] <= mars_radius]
impactor_iron_within_mars = particles_within_mars[particles_within_mars['tag'] == 3]
# sort particles by radius
impactor_iron_within_mars = impactor_iron_within_mars.sort_values(by=['radius'])
# create a column of the cumulative sum of the mass of the particles within each radius
impactor_iron_within_mars['cumulative_mass'] = impactor_iron_within_mars['mass'].cumsum() / impactor_iron_within_mars['mass'].sum()

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(impactor_iron_within_mars['radius'], impactor_iron_within_mars['cumulative_mass'], color='black', linewidth=2)
ax.set_xlabel("Radius (km)")
ax.set_ylabel("Cumulative Mass Fraction of Impactor Iron Particles")
ax.set_title(f"{round(cf.sim_time * 0.000277778, 2)} hrs.")
ax.grid()
plt.tight_layout()
plt.savefig("iron_cdf.png", dpi=200)
