import pandas as pd
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.pyplot as plt

from src.vapor import get_all_particle_vapor_fractions_from_formatted

plt.style.use("dark_background")


phase_path = "src/phase_data/forstSTS__vapour_curve.txt"
eos_path = "src/phase_data/forst_STS.rho_u.txt"
output = "/home/theia/scotthull/1M_high_rho_cutoff/formatted_gi_new_eos_b_073_high_rho_cutoff_1M"
start_time = 0
end_time = 2900
increment = 1000


