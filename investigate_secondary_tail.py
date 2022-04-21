import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 8, })
# plt.style.use("dark_background")
plt.style.use('seaborn-colorblind')

iteration = 190
angles = ['b073', 'b075']
cutoff_densities = [5, 500, 1000, 2000]
vsquare_scale = 6e7
hsquare_scale = 1e7
base_path = "/home/theia/scotthull/Paper1_SPH/gi/"
new_phase_path = "src/phase_data/forstSTS__vapour_curve.txt"
old_phase_path = "src/phase_data/duniteN__vapour_curve.txt"

new_phase_df = pd.read_fwf(new_phase_path, skiprows=1,
                           names=["temperature", "density_sol_liq", "density_vap", "pressure",
                                  "entropy_sol_liq", "entropy_vap"])
old_phase_df = pd.read_fwf(old_phase_path, skiprows=1,
                           names=["temperature", "density_sol_liq", "density_vap", "pressure",
                                  "entropy_sol_liq", "entropy_vap"])


def get_all_sims(angle, high=True):
    fformat = "{}_{}_{}"
    tformat = "{}{}{}"
    names = []
    titles = []
    for runs in ["new", "old"]:
        n = "n"
        if runs == "old":
            n = "o"
        for cd in cutoff_densities:
            output_name = fformat.format(cd, angle, runs)
            title_name = tformat.format(cd, angle, n)
            titles.append(title_name)
            names.append(output_name)
    if high:
        output_name = fformat.format(5, angle, "new") + "_high"
        names.append(output_name)
        title_name = tformat.format(5, angle, "n") + "-high"
        titles.append(title_name)
    return names, titles

def __build_scene(iteration, times, dfs, sims, titles, impact_angles, target_coms, impactor_coms, mom_vectors,
                  to_save_path, r_vectors, v_vectors):
    num_new = len([i for i in sims if "new" in i])
    num_old = len([i for i in sims if "old" in i])
    num_rows = max([num_new, num_old])
    plt.style.use("dark_background")
    fig, axs = plt.subplots(num_rows, 2, figsize=(16, 32), sharex='all',
                            gridspec_kw={"hspace": 0.10, "wspace": 0.12})
    fig.patch.set_facecolor('xkcd:black')
    axs = axs.flatten()
    for ax in axs:
        ax.set_xlim(-hsquare_scale, hsquare_scale)
        ax.set_ylim(-vsquare_scale, 0)
    index_new, index_old = 0, 1
    for s, t in zip(sims, titles):
        df, impact_angle, target_com, impactor_com, time, mom_vector, r_vector, v_vector = dfs[t][-1], impact_angles[t][-1], \
                                                                       target_coms[t][-1], impactor_coms[t][-1], \
                                                                       times[t][-1], mom_vectors[t][-1], \
                                                                        r_vectors[t][-1], v_vectors[t][-1]
        df = df[df['z'] < 0]
        target_silicate = df[df['tag'] == 0]
        target_iron = df[df['tag'] == 1]
        impactor_silicate = df[df['tag'] == 2]
        impactor_iron = df[df['tag'] == 3]
        t1 = ["Target Silicate", "Target Iron", "Impactor Silicate", "Impactor Iron"]
        t2 = [target_silicate, target_iron, impactor_silicate, impactor_iron]

        to_index = index_new
        if "old" in s:
            to_index = index_old
        for a, b, in zip(t1, t2):
            axs[to_index].scatter(
                b['x'], b['y'], s=1, marker='.', label=a
            )

        axs[to_index].scatter(
            impactor_com[0], impactor_com[1], s=60, c='pink', marker="*", label="Impactor COM"
        )
        axs[to_index].scatter(
            target_com[0], target_com[1], s=60, c='green', marker="*", label="Target COM"
        )
        axs[to_index].plot(
            [target_com[0], impactor_com[0]], [target_com[1], impactor_com[1]], linewidth=2.0, color="aqua"
        )
        axs[to_index].plot(
            [target_com[0], target_com[0]], [target_com[1], impactor_com[1]], linewidth=2.0, color="aqua"
        )
        axs[to_index].plot(
            [impactor_com[0], target_com[0]], [impactor_com[1], impactor_com[1]], linewidth=2.0, color="aqua"
        )

        r_vector = np.array(impactor_com) - np.array(target_com)
        axs[to_index].quiver(impactor_com[0], impactor_com[1], mom_vector[0], mom_vector[1], color='green', label="Spec. Mom. Vector")
        axs[to_index].quiver(target_com[0], target_com[1], r_vector[0], r_vector[1], color='yellow', label="Radial Vector")


        axs[to_index].set_title("{} {} hrs. ({} - {} deg.)".format(t, time, iteration, round(impact_angle, 2)))
        if "old" in s:
            index_old += 2
        else:
            index_new += 2
    axs[0].legend(loc='upper left')
    plt.savefig(to_save_path + "/{}.png".format(iteration), format='png', dpi=200)


