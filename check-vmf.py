import os
import pandas as pd
import matplotlib.pyplot as plt

from src.interpolation import NearestNeighbor1D

# use seaborn colorblind
plt.style.use('seaborn-colorblind')
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

# read in the data
data_path = r"C:\Users\Scott\OneDrive\Desktop\500b073S_1800.csv"
df = pd.read_csv(data_path)
disk_df = df[df['label'] == "DISK"]
disk_df = disk_df[disk_df['tag'] % 2 == 0]
disk_df['entropy_w_circ'] = disk_df['entropy'] + disk_df['circ_entropy_delta']
disk_df = disk_df[disk_df['circ_entropy_delta'] < 5000]

# read in the phase curve data
new_phase_path = "src/phase_data/forstSTS__vapour_curve.txt"
phase_df = pd.read_fwf(new_phase_path, skiprows=1,
                       names=["temperature", "density_sol_liq", "density_vap", "pressure",
                              "entropy_sol_liq", "entropy_vap"])


def shade_plot(s, ax):
    cp = max(phase_df['temperature'])

    ax.plot(
        phase_df['entropy_sol_liq'],
        phase_df['temperature'],
        linewidth=2.0,
        color='black'
    )
    ax.plot(
        phase_df['entropy_vap'],
        phase_df['temperature'],
        linewidth=2.0,
        color='black'
    )
    ax.fill_between(
        x=phase_df['entropy_vap'],
        y1=phase_df['temperature'],
        y2=cp,
        color=colors[-1],
        alpha=0.2,
        # label="100% Vapor"
    )
    ax.fill_between(
        x=phase_df['entropy_sol_liq'],
        y1=phase_df['temperature'],
        y2=cp,
        color=colors[-2],
        alpha=0.2,
        # label="100% Liquid"
    )
    ax.fill_between(
        x=sorted(list(phase_df['entropy_sol_liq']) + list(phase_df['entropy_vap'])),
        y1=cp,
        y2=1e10,
        color=colors[-3],
        alpha=0.2,
        # label="Supercritical"
    )
    ax.fill_between(
        x=phase_df['entropy_sol_liq'],
        y1=phase_df['temperature'],
        color=colors[-4],
        edgecolor="none",
        alpha=0.2,
        # label="Mixed"
    )
    ax.fill_between(
        x=phase_df['entropy_vap'],
        y1=phase_df['temperature'],
        color=colors[-4],
        edgecolor="none",
        alpha=0.2,
    )

    ax.set_xlabel("Entropy (J/K)")
    ax.set_ylabel("Temperature (K)")


def get_phase_data_for_particle(s, t, nearest_neighbor):
    nearest_temperature_index = nearest_neighbor.neighbor_index(given_val=t,
                                                                array=list(phase_df['temperature']))
    entropy_liq = phase_df['entropy_sol_liq'][nearest_temperature_index]
    entropy_vap = phase_df['entropy_vap'][nearest_temperature_index]
    return entropy_liq, entropy_vap


def calc_vmf(data_df, phase_df):
    nearest_neighbor = NearestNeighbor1D()
    supercritical = max(phase_df['temperature'])
    # get the location of the phase curve that corresponds to the supercritical point
    supercritical_idx = phase_df[phase_df['temperature'] == supercritical].index[0]
    # get the entropy values for the supercritical point
    supercritical_entropy = phase_df.iloc[supercritical_idx]['entropy_sol_liq']
    vmf = 0.0
    num_particles = 0
    particles = {'liq': [], 'vap': [], 'mix': [], 'supercritical': []}
    for s, t in zip(disk_df['entropy_w_circ'], data_df['temperature']):
        nearest_temperature_index = nearest_neighbor.neighbor_index(given_val=t,
                                                                    array=list(phase_df['temperature']))
        entropy_liq = phase_df['entropy_sol_liq'][nearest_temperature_index]
        entropy_vap = phase_df['entropy_vap'][nearest_temperature_index]
        if t >= supercritical:
            vmf += 1.0
            particles['supercritical'].append((s, t, entropy_liq, entropy_vap))
        elif s < entropy_liq:
            vmf += 0.0
            particles['liq'].append((s, t, entropy_liq, entropy_vap))
        elif entropy_liq <= s <= entropy_vap:
            vmf += (s - entropy_liq) / (entropy_vap - entropy_liq)
            particles['mix'].append((s, t, entropy_liq, entropy_vap))
        elif s > entropy_vap:
            vmf += 1.0
            particles['vap'].append((s, t, entropy_liq, entropy_vap))
        num_particles += 1
    return vmf / num_particles, (supercritical_entropy, supercritical), particles


fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
shade_plot(s="new", ax=ax)
ax.grid(alpha=0.4)
ax.set_xlim(1800, 12000)
ax.set_ylim(0, 12500)
# ax.scatter(
#     disk_df['entropy'], disk_df['temperature'], s=1
# )
# ax.scatter(
#     disk_df['entropy_w_circ'], disk_df['temperature'], s=1
# )
# calculate VMF and annotate on upper right of plot
vmf = calc_vmf(disk_df, phase_df)
ax.scatter(
    vmf[1][0], vmf[1][1], s=100, color='black', marker='x', label="Supercritical Point"
)
ax.scatter(
    [p[0] for p in vmf[2]['liq']], [p[1] for p in vmf[2]['liq']], s=8, label="Liquid"
)
ax.scatter(
    [p[0] for p in vmf[2]['vap']], [p[1] for p in vmf[2]['vap']], s=8, label="Vapor"
)
ax.scatter(
    [p[0] for p in vmf[2]['mix']], [p[1] for p in vmf[2]['mix']], s=8, label="Mixed"
)
ax.scatter(
    [p[0] for p in vmf[2]['supercritical']], [p[1] for p in vmf[2]['supercritical']], s=8, label="Supercritical"
)
# pick 3 random particles from the data set and get their phase data
random_particles = disk_df[disk_df['entropy_w_circ'] > 5000].sample(n=3)
for s, t in zip(random_particles['entropy_w_circ'], random_particles['temperature']):
    entropy_liq, entropy_vap = get_phase_data_for_particle(s, t, nearest_neighbor=NearestNeighbor1D())
    # plot a line from the particle to the phase curve
    ax.plot(
        [s, entropy_liq],
        [t, t],
        color='blue',
        linewidth=1.0,
    )
    ax.plot(
        [s, entropy_vap],
        [t, t],
        color='green',
        linewidth=1.0,
    )
    # scatter the particle
    ax.scatter(
        s, t, s=100, color='magenta', marker='o'
    )

ax.annotate(
    f"VMF = {vmf[0] * 100.0:.2f}",
    xy=(0.95, 0.95),
    xycoords='axes fraction',
    horizontalalignment='right',
    verticalalignment='top',
    fontsize=16
)
ax.legend(loc='upper left')
# make the dots in the legend larger
for handle in ax.get_legend().legendHandles:
    handle._sizes = [100]
plt.show()
