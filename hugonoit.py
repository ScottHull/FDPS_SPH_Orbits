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
from src.hugonoit import entropy_increase_analytical_hugonoit, entropy_increase_table_hugonoit

plt.style.use("dark_background")
fig, axs = plt.subplots(1, 1, figsize=(10, 10), sharey='all', sharex='all',
                        gridspec_kw={"hspace": 0.0, "wspace": 0.12})
fig.patch.set_facecolor('xkcd:black')

dS, P = entropy_increase_analytical_hugonoit(
    v_p_list=np.arange(0, 1000, 1),
    s=1.56,
    T=2000,
    P_i=50 * (10 ** 9),
    C_0=9 / 1000,
    rho_i=4500
)
axs.plot(
    P,
    dS,
    color='magenta',
    linewidth=1.0
)
dS, P = entropy_increase_table_hugonoit(hugonoit_path="src/phase_data/forstSTS__hugoniot.txt")
# axs.plot(
#     P,
#     dS,
#     color='aqua',
#     linewidth=1.0
# )

axs.set_xlabel("Pressure (GPa)")
axs.set_ylabel("dS (10^3)")
# axs.set_xlim(0, 300)
# plt.savefig("hugonoit.png", format='png')
plt.show()
