import pandas as pd
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.pyplot as plt

plt.style.use("dark_background")


phase_path = "src/phase_data/forstSTS__vapour_curve.txt"
eos_path = "src/phase_data/forst_STS.rho_u.txt"
output = "/Users/scotthull/Desktop/1000.csv"

min_normalize_parameter = 1
max_normalize_parameter = 10
normalizer = Normalize(min_normalize_parameter, max_normalize_parameter)
cmap = cm.get_cmap('jet')

phase_df = pd.read_fwf(phase_path, skiprows=1,
                           names=["temperature", "density_sol_liq", "density_vap", "pressure",
                                  "entropy_sol_liq", "entropy_vap"])
eos_df = pd.read_fwf(eos_path, skiprows=2, header=None)
df = pd.read_csv(output, skiprows=2)
disk = df[df['end_label'] == "DISK"]

fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)

ax.plot(
    phase_df['entropy_sol_liq'],
    phase_df['temperature'],
    linewidth=2.0,
    label="Liquid",
)
ax.plot(
    phase_df['entropy_vap'],
    phase_df['temperature'],
    linewidth=2.0,
    label="Vapor"
)
ax.scatter(
    eos_df[5],
    eos_df[2],
    s=2,
    color='white'
)
ax.scatter(
    disk['entropy'],
    disk['temperature'],
    s=2,
    color=[cmap(normalizer(i / (6371 * 1000))) for i in disk['radius']]
)
ax.set_xlim(0, 10000)
ax.set_ylim(0, 10000)
ax.set_xlabel("Entropy")
ax.set_ylabel("Temperature")
ax.grid(alpha=0.4)
ax.legend(loc='upper right')
ax.set_title("ForstSTS (New EoS)")

sm = cm.ScalarMappable(norm=normalizer, cmap=cmap)
sm.set_array([])
cbaxes = inset_axes(ax, width="30%", height="3%", loc=2, borderpad=1.8)
cbar = plt.colorbar(sm, cax=cbaxes, orientation='horizontal')
cbar.ax.tick_params(labelsize=6)
cbar.ax.set_title("Radius (R_E)", fontsize=6)

plt.show()
