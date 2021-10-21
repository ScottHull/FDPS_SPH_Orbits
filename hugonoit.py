#!/usr/bin/env python3
import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
from matplotlib.colors import Normalize
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.font_manager as fm

from src.identify import ParticleMap
from src.combine import CombineFile
from src.time import get_nearest_iteration_to_time, seconds_to_hours, get_all_iterations_and_times
from src.new_and_old_eos import get_particles, scatter, plot, main_plotting_loop
from src.animate import animate
from src.hugonoit import entropy_increase_analytical_hugonoit

plt.style.use("dark_background")
fig, axs = plt.subplots(1, 1, figsize=(10, 10), sharey='all', sharex='all',
                        gridspec_kw={"hspace": 0.0, "wspace": 0.12})
fig.patch.set_facecolor('xkcd:black')

dS, P = entropy_increase_analytical_hugonoit(
    v_p_list=np.arange(0, 20, 0.1),
    s=1.56,
    T=2000,
    P_i=50,
    C_0=9,
    rho_i=4500
)

axs.flatten()[0].plot(
    P,
    dS,
    color='magenta',
    linwidth=1.0
)
axs.flatten()[0].set_xlabel("Pressure")
axs.flatten()[0].set_ylabel("dS")
plt.savefig("hugonoit.png", format='png')
