import pandas as pd
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

fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
ax.scatter(
    df['x'],
    df['y'],
    s=2,
    color=[cmap(normalizer(i)) for i in df['entropy']]
)
ax.set_xlim(-square_scale, square_scale)
ax.set_ylim(-square_scale, square_scale)
sm = cm.ScalarMappable(norm=normalizer, cmap=cmap)
sm.set_array([])
cbaxes = inset_axes(ax, width="30%", height="3%", loc=2, borderpad=1.8)
cbar = plt.colorbar(sm, cax=cbaxes, orientation='horizontal')
cbar.ax.tick_params(labelsize=6)
cbar.ax.set_title("Entropy", fontsize=6)

plt.show()
