import pandas as pd
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.pyplot as plt

from src.vapor import get_all_particle_vapor_fractions_from_formatted

plt.style.use("dark_background")


phase_path = "src/phase_data/forstSTS__vapour_curve.txt"
eos_path = "src/phase_data/forst_STS.rho_u.txt"
output = "/Users/scotthull/Desktop/2900.csv"
square_scale = 6e7

min_normalize_parameter = 0
max_normalize_parameter = 1
normalizer = Normalize(min_normalize_parameter, max_normalize_parameter)
cmap = cm.get_cmap('jet')

df = pd.read_csv(output, skiprows=2)
disk_particles = df[df['end_label'] == "DISK"]
sil_disk_particles = disk_particles[disk_particles['tag'] % 2 == 0]
vmfs = get_all_particle_vapor_fractions_from_formatted(sil_disk_particles, phase_path, end_label=True)

fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
ax.scatter(
    sil_disk_particles['x'],
    sil_disk_particles['y'],
    color=[cmap(normalizer(f)) for s, t, f, rho, x, y, z in vmfs]
)

ax.scatter(
    [x for s, t, f, rho, x, y, z in vmfs if f == 0],
    [y for s, t, f, rho, x, y, z in vmfs if f == 0],
    color='blue',
    label="100% Liquid"
)
ax.scatter(
    [x for s, t, f, rho, x, y, z in vmfs if 0 < f < 1],
    [y for s, t, f, rho, x, y, z in vmfs if 0 < f < 1],
    color='yellow',
    label="Mixed"
)
ax.scatter(
    [x for s, t, f, rho, x, y, z in vmfs if f == 1],
    [y for s, t, f, rho, x, y, z in vmfs if f == 1],
    color='red',
    label="100% Vapor"
)


ax.set_xlim(-square_scale, square_scale)
ax.set_ylim(-square_scale, square_scale)
ax.legend(loc='upper right')

# sm = cm.ScalarMappable(norm=normalizer, cmap=cmap)
# sm.set_array([])
# cbaxes = inset_axes(ax, width="30%", height="3%", loc=2, borderpad=1.8)
# cbar = plt.colorbar(sm, cax=cbaxes, orientation='horizontal')
# cbar.ax.tick_params(labelsize=6)
# cbar.ax.set_title("VMF", fontsize=6)

plt.show()
