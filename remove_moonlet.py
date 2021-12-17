import pandas as pd
from statistics import mean
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.pyplot as plt

plt.style.use("dark_background")


phase_path = "src/phase_data/forstSTS__vapour_curve.txt"
eos_path = "src/phase_data/forst_STS.rho_u.txt"
output = "/Users/scotthull/Desktop/1000.csv"
square_scale = 4e7

min_normalize_parameter = 1000
max_normalize_parameter = 6000
normalizer = Normalize(min_normalize_parameter, max_normalize_parameter)
cmap = cm.get_cmap('jet')

df = pd.read_csv(output, skiprows=2)
disk = df[df['end_label'] == "DISK"]
disk_low_rho = disk[disk['density'] == 2000]
print(len(disk.index), len(disk_low_rho.index), mean(disk_low_rho['temperature']))

fig, axs = plt.subplots(1, 2, figsize=(16, 9), gridspec_kw={"hspace": 0.10, "wspace": 0.14})
ax1, ax2 = axs.flatten()
ax1.scatter(
    disk['x'],
    disk['y'],
    s=4,
    color=[cmap(normalizer(i)) for i in disk['temperature']]
)
ax1.set_title("All Disk Particles")
ax2.scatter(
    disk_low_rho['x'],
    disk_low_rho['y'],
    s=4,
    color=[cmap(normalizer(i)) for i in disk_low_rho['temperature']]
)
ax2.set_title("Disk Particles w/ Density == 2000")

for ax in axs.flatten():
    ax1.set_xlim(-square_scale, square_scale)
    ax1.set_ylim(-square_scale, square_scale)
sm = cm.ScalarMappable(norm=normalizer, cmap=cmap)
sm.set_array([])
cbaxes = inset_axes(ax1, width="30%", height="3%", loc=2, borderpad=1.8)
cbar = plt.colorbar(sm, cax=cbaxes, orientation='horizontal')
cbar.ax.tick_params(labelsize=6)
cbar.ax.set_title("temperature", fontsize=6)

plt.show()
