#!/usr/bin/env python3
import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.font_manager as fm

from src.identify import ParticleMap
from src.combine import CombineFile
from src.time import get_nearest_iteration_to_time, seconds_to_hours, get_all_iterations_and_times
from src.new_and_old_eos import get_particles, scatter, plot, main_plotting_loop
from src.animate import animate

min_iteration = 0
max_iteration = 3000
sample_interval = 5
min_normalize_entropy = 0
max_normalize_entropy = 8000
min_normalize_internal_energy = 0
max_normalize_internal_energy = 1e7
min_normalize_temperature = 0
max_normalize_temperature = 10000
square_scale = 1e7
number_processes = 100
path = "/home/theia/scotthull/sph_simulations/gi_new_eos"
to_path_entropy = "/home/theia/scotthull/FDPS_SPH_Orbits/animate_entropy"
to_path_internal_energy = "/home/theia/scotthull/FDPS_SPH_Orbits/animate_internal_energy"
to_path_temperature = "/home/theia/scotthull/FDPS_SPH_Orbits/animate_temperature"
entropy_parameter = "entropy"
internal_energy_parameter = "internal_energy"
temperature_parameter = "temperature"

for i in [to_path_entropy, to_path_internal_energy, to_path_temperature]:
    if os.path.exists(i):
        shutil.rmtree(i)
    os.mkdir(i)

plt.style.use("dark_background")
cmap = cm.get_cmap('jet')
normalizer_entropy = Normalize(min_normalize_entropy, max_normalize_entropy)
normalizer_internal_energy = Normalize(min_normalize_internal_energy, max_normalize_internal_energy)
normalizer_temperature = Normalize(min_normalize_temperature, max_normalize_temperature)

for time in np.arange(min_iteration, max_iteration + sample_interval, sample_interval):
    particles, sample_time = get_particles(path=path, number_processes=number_processes, time=time)
    fig, axs = plt.subplots(2, 2, figsize=(10, 10),
                                gridspec_kw={"hspace": 0.0, "wspace": 0.08})

