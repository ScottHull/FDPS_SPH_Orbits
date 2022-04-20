import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 14, })
# plt.style.use("dark_background")
plt.style.use('seaborn-colorblind')

iteration = 1800
angles = ['b073', 'b075']
cutoff_densities = [5, 500, 1000, 2000]
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
            if cd == 5 and high and runs == "new":
                output_name = fformat.format(cd, angle, runs) + "_high"
                names.append(output_name)
                title_name = tformat.format(cd, angle, n) + "-high"
                titles.append(title_name)
    return names, titles


fig, axs = plt.subplots(2, 2, figsize=(16, 9), gridspec_kw={"hspace": 0.24, "wspace": 0.20})
axs = axs.flatten()
b073n_index, b073o_index, b075n_index, b075o_index = 0, 1, 2, 3
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

for ax in axs:
    ax.grid(alpha=0.4)

for ax in [axs[0], axs[2]]:
    critical_point = max(new_phase_df['temperature'])
    ax.plot(
        new_phase_df['entropy_sol_liq'],
        new_phase_df['temperature'],
        linewidth=2.0,
        color='black'
    )
    ax.set_ylabel("Temperature (K)")
for ax in [axs[1], axs[3]]:
    ax.plot(
        old_phase_df['entropy_sol_liq'],
        old_phase_df['temperature'],
        linewidth=2.0,
        color='black'
    )
    ax.fill_between(
        x=old_phase_df['entropy_sol_liq'],
        y1=old_phase_df['temperature'],
        color='red',
        alpha=0.2,
        label="Mixed"
    )
    ax.fill_between(
        x=old_phase_df['entropy_sol_liq'],
        y1=old_phase_df['temperature'],
        color='red',
        alpha=0.2,
        label="Mixed"
    )
for ax in [axs[2], axs[3]]:
    ax.set_xlabel("Entropy")


def plot_phase_diagrams():
    for angle in angles:
        high = True
        if angle == "b075":
            high = False
        sims, titles = get_all_sims(angle=True, high=high)
        for s, t in zip(sims, titles):
            to_path = base_path + "{}/circularized_{}/{}.csv".format(s, s, iteration)
            df = pd.read_csv(to_path)
            disk = df[df['label'] == "DISK"]
            disk_filtered = disk[disk['circ_entropy_delta'] < 5000]
            temp, entropy = disk_filtered['temperature'], disk_filtered['entropy'] + disk_filtered['circ_entropy_delta']
            cd = int(s.split("_")[0])
            color = colors[cutoff_densities.index(cd)]
            phase_curve = new_phase_df
            if "old" in s:
                phase_curve = old_phase_df
            to_index = 0
            if "b073" in s and "new" in s:
                to_index = 0
            elif "b073" in s and "old" in s:
                to_index = 1
            elif "b075" in s and "new" in s:
                to_index = 2
            else:
                to_index = 3
            axs[to_index].scatter(
                temp, entropy, s=5, color=color, label=t
            )
