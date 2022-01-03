import pandas as pd
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.pyplot as plt

plt.style.use("dark_background")

base_path = "/home/theia/scotthull/Paper1_SPH/tar-imp/"
output_headers = [
    "id", "tag", "x", "y", "z", "vx", "vy", "vz", "density", "internal_energy", "pressure", "potential_energy",
    "entropy", "temperature"
]

densities = [
    5, 500, 1000, 2000
]

run_names = [
    ["{}_new".format(i), "{}_old".format(i)] for i in densities
]

types = ['tar', 'imp']

fig, axs = plt.subplots(len(run_names), len(run_names[0]), figsize=(16, 32), sharex='all', sharey='all',
                            gridspec_kw={"hspace": 0.10, "wspace": 0.10})
axs = axs.flatten()

for t in types:
    at_index = 0
    for new, old in run_names:
        try:
            new_path, old_path = base_path + "{}".format(new), base_path + "{}".format(old)
            new_df, old_df = pd.read_csv(new_path + "/{}.dat".format(t), skiprows=2, names=output_headers), \
                             pd.read_csv(old_path + "/{}.dat".format(t), skiprows=2, names=output_headers)
            new_radius = [(i ** 2 + j ** 2 + k ** 2) ** (1 / 2) for i, j, k in zip(new_df['x'], new_df['y'], new_df['z'])]
            old_radius = [(i ** 2 + j ** 2 + k ** 2) ** (1 / 2) for i, j, k in zip(old_df['x'], old_df['y'], old_df['z'])]
            new_df['radius'] = new_radius
            old_df['radius'] = old_radius

            axs[at_index].scatter(
                [i / 1000.0 for i in new_df[new_df['tag'] % 2 != 0]['radius']],
                new_df[new_df['tag'] % 2 == 0]['density'],
                s=2,
                label="Iron"
            )
            axs[at_index].scatter(
                [i / 1000.0 for i in new_df[new_df['tag'] % 2 == 0]['radius']],
                new_df[new_df['tag'] % 2 == 0]['density'],
                s=2,
                label="Silicate"
            )
            at_index += 1

            axs[at_index].scatter(
                [i / 1000.0 for i in old_df[old_df['tag'] % 2 != 0]['radius']],
                old_df[old_df['tag'] % 2 == 0]['density'],
                s=2,
                label="Iron"
            )
            axs[at_index].scatter(
                [i / 1000.0 for i in old_df[old_df['tag'] % 2 == 0]['radius']],
                old_df[old_df['tag'] % 2 == 0]['density'],
                s=2,
                label="Silicate"
            )
            at_index += 1
        except Exception as e:
            print(e)


    for ax in axs:
        ax.grid(alpha=0.4)
    axs[0].legend()

    plt.savefig("{}_verify.png".format(t), format='png')

