import pandas as pd
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.pyplot as plt

plt.style.use("dark_background")

base_path = "/home/theia/scotthull/Paper1_SPH/tar-imp/"
output_headers = [
    "id", "tag", "mass", "x", "y", "z", "vx", "vy", "vz", "density", "internal_energy", "pressure", "potential_energy",
    "entropy", "temperature"
]

densities = [
    5, 500, 1000, 2000
]

run_names = [
    ["{}_new".format(i), "{}_old".format(i)] for i in densities
]

types = ['tar', 'imp']
types_formal = ["Target", "Impactor"]

def add_annotation(ax, text, min_rho, max_rho):
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    coord = (xmax - (0.48
                     * xmax), ymax - (0.25 * ymax))
    density, new_or_old = text.split("_")
    t = "{} kg/m3\nMin. Density: {}\nMax. Density: {}".format(density, new_or_old.capitalize(),
                                                                   int(min_rho), int(max_rho))
    ax.annotate(t, coord, fontsize=18)


for index, t in enumerate(types):
    format_label = types_formal[index]
    fig, axs = plt.subplots(len(run_names), len(run_names[0]), figsize=(16, 32), sharex='all', sharey='all',
                            gridspec_kw={"hspace": 0.10, "wspace": 0.10})
    axs = axs.flatten()
    at_index = 0
    has_legend = False
    for new, old in run_names:
        try:
            new_path, old_path = base_path + "{}".format(new), base_path + "{}".format(old)
            new_df, old_df = pd.read_csv(new_path + "/{}.dat".format(t), skiprows=2, names=output_headers, delimiter="\t"), \
                             pd.read_csv(old_path + "/{}.dat".format(t), skiprows=2, names=output_headers, delimiter="\t")
            new_radius = [(i ** 2 + j ** 2 + k ** 2) ** (1 / 2) for i, j, k in zip(new_df['x'], new_df['y'], new_df['z'])]
            old_radius = [(i ** 2 + j ** 2 + k ** 2) ** (1 / 2) for i, j, k in zip(old_df['x'], old_df['y'], old_df['z'])]
            new_df['radius'] = new_radius
            old_df['radius'] = old_radius

            axs[at_index].scatter(
                [i / 1000.0 for i in new_df[new_df['tag'] % 2 != 0]['radius']],
                new_df[new_df['tag'] % 2 != 0]['density'],
                s=2,
                label="Iron",
                color="#8dd3c7"
            )
            axs[at_index].scatter(
                [i / 1000.0 for i in new_df[new_df['tag'] % 2 == 0]['radius']],
                new_df[new_df['tag'] % 2 == 0]['density'],
                s=2,
                label="Silicate",
                color="#feffb3"
            )
            if not has_legend:
                legend = axs[at_index].legend()
                has_legend = True
                for handle in legend.legendHandles:
                    handle.set_sizes([3.0])
            add_annotation(ax=axs[at_index], text=new, min_rho=min(new_df['density']), max_rho=max(new_df['density']))
            if at_index % 2 == 0:
                axs[at_index].set_ylabel("Density (kg/m3)")
            at_index += 1

            axs[at_index].scatter(
                [i / 1000.0 for i in old_df[old_df['tag'] % 2 != 0]['radius']],
                old_df[old_df['tag'] % 2 != 0]['density'],
                s=2,
                # label="Iron",
                # color="#8dd3c7"
            )
            axs[at_index].scatter(
                [i / 1000.0 for i in old_df[old_df['tag'] % 2 == 0]['radius']],
                old_df[old_df['tag'] % 2 == 0]['density'],
                s=2,
                # label="Silicate",
                # color="#feffb3"
            )
            add_annotation(ax=axs[at_index], text=old, min_rho=min(old_df['density']), max_rho=max(old_df['density']))
            at_index += 1
        except Exception as e:
            at_index += 2
            print(e)


    for ax in axs:
        ax.grid(alpha=0.4)

    axs[-2].set_xlabel(r'Radius $R_{\bigoplus}$')
    axs[-1].set_xlabel(r'Radius $R_{\bigoplus}$')
    axs[0].set_title("{} (New EoS)".format(format_label))
    axs[1].set_title("{} (Old EoS)".format(format_label))

    plt.savefig("{}_verify.png".format(t), format='png')

