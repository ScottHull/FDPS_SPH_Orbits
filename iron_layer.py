#!/usr/bin/env python3
from scipy.interpolate import UnivariateSpline
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from src.new_and_old_eos import get_particles

"""
We want to model the iron layer after the GI as a spherical shell.
Volume of a sphere: 4/3 pi r^3
Volume of a spherical shell: 4/3 pi (r_2 - r_1)^3
We need to calculate r_2 - r_1 (can't get it directly from simulation)
Know that density is rho = M / V, so the density of the shell is rho = M / [4/3 pi (r_2 - r_1)^3]
We can get M from the simulation, but what about density?

THERMAL EXPANSIVITY (alpha) : [1 / K]
TEMPERATURE (T) : [K]
DENSITY (rho) : [kg/m3]

DELTA RHO = RHO_0 * ALPHA * DELTA T
RHO_0_iron = 7500

"""

time = 1500
path = "/home/theia/scotthull/1M/gi_new_eos_b_073_formatted"

plt.style.use("dark_background")

labels = {
    0: "Target Silicate",
    1: "Target Iron",
    2: "Impactor Silicate",
    3: "Impactor Iron"
}

f = path + "/{}.csv".format(time)
df = pd.read_csv(f, skiprows=2).tolist('index')
impactor_iron = df.loc[(df['tag'] == 3)
print(impactor_iron)
