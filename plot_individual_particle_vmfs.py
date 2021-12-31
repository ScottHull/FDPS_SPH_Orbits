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

min_normalize_parameter = 0
max_normalize_parameter = 1
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
    phase_df['density_sol_liq'],
    linewidth=2.0,
    label="Liquid",
)
ax.plot(
    phase_df['entropy_vap'],
    phase_df['density_vap'],
    linewidth=2.0,
    label="Vapor"
)

particles_w_vmf = get_all_particle_vapor_fractions_from_formatted(df=df, phase_path=phase_path, end_label=True)
all_liq = [i for i in particles_w_vmf if i[2] == 0]
all_vap = [i for i in particles_w_vmf if i[2] == 1]
mixed = [i for i in particles_w_vmf if 0 < i[2] < 1]
# ax.scatter(
#     [s for s, t, f, rho in particles_w_vmf],
#     [t for s, t, f, rho in particles_w_vmf],
#     s=10,
#     color=[cmap(normalizer(f)) for s, t, f, rho in particles_w_vmf]
# )
ax.scatter(
    [s for s, t, f, rho in all_liq],
    [rho for s, t, f, rho in all_liq],
    s=10,
    color='blue',
    label="100% Liquid"
)
ax.scatter(
    [s for s, t, f, rho in all_vap],
    [rho for s, t, f, rho in all_vap],
    s=10,
    color='red',
    label="100% Vapor"
)
ax.scatter(
    [s for s, t, f, rho in mixed],
    [rho for s, t, f, rho in mixed],
    s=10,
    color='yellow',
    label="Mixed"
)
ax.set_xlim(0, 12000)
ax.set_ylim(0, 10000)
ax.set_xlabel("Entropy")
ax.set_ylabel("Density")
ax.grid(alpha=0.4)
ax.legend(loc='upper right')
ax.set_title("ForstSTS (New EoS)")

sm = cm.ScalarMappable(norm=normalizer, cmap=cmap)
sm.set_array([])
cbaxes = inset_axes(ax, width="30%", height="3%", loc=2, borderpad=1.8)
cbar = plt.colorbar(sm, cax=cbaxes, orientation='horizontal')
cbar.ax.tick_params(labelsize=6)
cbar.ax.set_title("VMF", fontsize=6)

plt.show()
