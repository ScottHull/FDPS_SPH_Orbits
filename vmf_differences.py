import os
import pandas as pd
import matplotlib.pyplot as plt

from src.interpolation import NearestNeighbor1D

# use seaborn colorblind
plt.style.use('seaborn-colorblind')
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

runs = [
    '5_b073_new', '5_b073_old', '5_b073_new_high', '500_b073_new', '500_b073_old', '1000_b073_new', '1000_b073_old',
    '2000_b073_new', '2000_b073_old',
    '5_b075_new', '5_b075_old', '500_b075_new', '500_b075_old',
    '1000_b075_new', '1000_b075_old', '2000_b075_new', '2000_b075_old', '2000_b075_old_low'
]
cutoff_densities = [5, 500, 1000, 2000]
max_iteration = 1800
base_path = "/home/theia/scotthull/Paper1_SPH/gi/"
new_phase_path = "src/phase_data/forstSTS__vapour_curve.txt"
old_phase_path = "src/phase_data/duniteN__vapour_curve.txt"
new_phase_df = pd.read_fwf(new_phase_path, skiprows=1,
                           names=["temperature", "density_sol_liq", "density_vap", "pressure",
                                  "entropy_sol_liq", "entropy_vap"])
old_phase_df = pd.read_fwf(old_phase_path, skiprows=1,
                           names=["temperature", "density_sol_liq", "density_vap", "pressure",
                                  "entropy_sol_liq", "entropy_vap"])


def get_phase_data_for_particle(s, t, nearest_neighbor, phase_df):
    supercritical = max(phase_df['temperature'])
    # get the location of the phase curve that corresponds to the supercritical point
    nearest_temperature_index = nearest_neighbor.neighbor_index(given_val=t,
                                                                array=list(phase_df['temperature']))
    entropy_liq = phase_df['entropy_sol_liq'][nearest_temperature_index]
    entropy_vap = phase_df['entropy_vap'][nearest_temperature_index]
    if t >= supercritical:
        return None
    if s < entropy_liq or s > entropy_vap:
        return None
    else:
        return (s - entropy_liq) / (entropy_vap - entropy_liq)


fig = plt.figure(figsize=(12, 16))
ax_new = fig.add_subplot(2, 1, 1)
ax_old = fig.add_subplot(2, 1, 2)
for ax in [ax_new, ax_old]:
    ax.set_xlabel("Entropy")
    ax.set_ylabel("s - s_liq / (s_vap - s_liq)")
    ax.grid(True)
for cd in cutoff_densities:
    c = colors[cutoff_densities.index(cd)]
    ax_old.scatter(
        [], [], c=c, label="{} kg/m^3".format(cd)
    )
for run in runs:
    cd = int(run.split('_')[0])
    c = colors[cutoff_densities.index(cd)]
    if "high" in run or "low" in run:
        c = len(colors)
    phase_df = new_phase_df if 'new' in run else old_phase_df
    ax = ax_new if 'new' in run else ax_old
    circ_path = base_path + f"{run}/circularized_{run}/{max_iteration}.csv"
    circ_df = pd.read_csv(circ_path)
    circ_df = circ_df[circ_df['circ_entropy_delta'] < 5000]
    circ_df = circ_df[circ_df['label'] == "DISK"]
    circ_df = circ_df[circ_df['tag'] % 2 == 0]
    circ_df['entropy_w_circ'] = circ_df['entropy'] + circ_df['circ_entropy_delta']
    d = [(s, t, get_phase_data_for_particle(s, t, NearestNeighbor1D, phase_df)) for s, t in zip(circ_df['entropy_w_circ'], circ_df['temperature'])]
    d = [x for x in d if x[2] is not None]

    ax.scatter(
        [x[0] for x in d],
        [x[3] for x in d],
        s=2,
        color=c
    )

plt.savefig("vmf_phase_curve_differences.png")
